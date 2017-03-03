# Converts a file path to Cygwin compliant paths
# Author : Sylvain Mareschal <mareschal@ovsa.fr>

cygwinPath = function(path) {
	if(grepl("^([A-Za-z]):", path)) {
		return(sub("^([A-Za-z]):", "/cygdrive/\\1", gsub("\\\\", "/", path)))
	} else {
		return(path)
	}
}


# Fuzzy matching of various length probes in multiple chromosomes using BLAT
# http://genome.ucsc.edu/goldenPath/help/blatSpec.html
# Author : Sylvain Mareschal <mareschal@ovsa.fr>

localize = function(
		probeFile,
		chromFiles,
		chromPattern = "^(.+)\\.[^\\.]+$",
		blatArgs = character(0),
		rawOutput = FALSE,
		noMulti = TRUE,
		noOverlap = TRUE,
		noPartial = TRUE,
		verbose = 2
		)
	{
	# Installation check
	blatExe <- sprintf("%s/blat/blat.exe", find.package("cghRA", quiet=TRUE, verbose=FALSE))
	if(!file.exists(blatExe)) stop("BLAT exe not found, see blatInstall()")
	
	# Inter-chromosome storage
	if(isTRUE(rawOutput)) {
		hits = NULL
	} else {
		name = character(0)
		chrom = character(0)
		strand = character(0)
		start = integer(0)
		end = integer(0)
		insertions = integer(0)
		deletions = integer(0)
		mismatches = integer(0)
		freeEnds = integer(0)
	}
	
	if(verbose == 1) message("Processing chromosome ", appendLF=FALSE)
	
	for(chromFile in chromFiles) {
		
		chromName = sub(chromPattern, "\\1", basename(chromFile))
		if(verbose == 1) message(chromName, appendLF=FALSE)
		if(verbose >= 2) message("[Processing chromosome \"", chromName, "\"]")
		
		# Temporary BLAT output file
		outFile <- tempfile()
		
		# BLAT call
		output <- suppressWarnings(
			system2(
				command = blatExe,
				args = c(blatArgs, shQuote(cygwinPath(chromFile)), shQuote(cygwinPath(probeFile)), shQuote(cygwinPath(outFile))),
				stderr = TRUE,
				stdout = TRUE
			)
		)
		
		# BLAT error handling
		if(!file.exists(outFile) || file.info(outFile)['size'] == 0) {
			stop("[BLAT] ", paste(output, collapse="\n"))
		}
		
		# File check
		header = scan(outFile, "", sep="\n", n=1, quiet=TRUE)
		if(header != "psLayout version 3") {
			stop("Unhandled output type (psLayout version 3 expected)")
		}
		
		# PSL 3 parsing
		tab = utils::read.table(
			file = outFile,
			sep = "\t",
			dec = ".",
			header = TRUE,
			stringsAsFactors = FALSE,
			col.names = c("matches", "misMatches", "repMatches", "nCount", "qNumInsert", "qBaseInsert", "tNumInsert", "tBaseInsert", "strand", "qName", "qSize", "qStart", "qEnd", "tName", "tSize", "tStart", "tEnd", "blockCount", "blockSizes", "qStarts", "tStarts"),
			colClasses = c("integer", "integer", "integer", "integer", "integer", "integer", "integer", "integer", "character", "character", "integer", "integer", "integer", "character", "integer", "integer", "integer", "integer", "character", "character", "character"),
			skip = 5
		)
		
		# Data appending
		if(isTRUE(rawOutput)) {
			hits = rbind(hits, tab)
		} else {
			name = c(name, tab$qName)
			chrom = c(chrom, rep(chromName, nrow(tab)))
			strand = c(strand, tab$strand)
			start = c(start, tab$tStart+1L)
			end = c(end, tab$tEnd)
			insertions = c(insertions, tab$qBaseInsert)
			deletions = c(deletions, tab$tBaseInsert)
			mismatches = c(mismatches, tab$misMatches)
			freeEnds = c(freeEnds, tab$qStart + (tab$qSize - tab$qEnd))
		}
		
		if(verbose == 1) message(", ", appendLF=FALSE)
		if(verbose >= 2) message(paste(output, collapse="\n"), "\n")
		
	}
	
	if(verbose == 1) message("done")
	
	# Hit table
	if(!isTRUE(rawOutput)) {
		
		if(verbose == 1) message("Post-processing")
		if(verbose >= 2) message("[Post-processing]")
		
		# Table building
		hits = data.frame(name, chrom, strand, start, end, insertions, deletions, mismatches, freeEnds, stringsAsFactors=FALSE)
		
		# Count
		if(verbose >= 2) {
			probeCount <- sum(grepl("^>", scan(probeFile, what="", sep="\n", quote=NULL, na.strings=NULL, quiet=TRUE, comment.char="")))
			hitCount <- length(unique(hits$name))
			message(probeCount - hitCount, " probes not found")
		}
		
		# Unique probe filtering
		if(isTRUE(noMulti)) {
			exclude = unique(hits[ duplicated(hits$name) , "name" ])
			if(verbose >= 2) message(length(exclude), " multi-hitting probes filtered out")
			if(length(exclude) > 0) {
				hits = hits[ -which(hits$name %in% exclude) ,]
			}
		}
		
		# Overlap filtering
		if(isTRUE(noOverlap)) {
			obj <- track.table(hits[, c("name","chrom","start","end","strand") ], warn=FALSE)
			exclude <- obj$extract(obj$cross(obj, type="count") > 1, "name")
			if(verbose >= 2) message(length(exclude), " overlapping probes filtered out")
			if(length(exclude) > 0) {
				hits = hits[ -which(hits$name %in% exclude) ,]
			}
		}
		
		# Partial filtering
		if(isTRUE(noPartial)) {
			exclude = unique(hits[ hits$insertions > 0 | hits$deletions > 0 | hits$mismatches > 0 | hits$freeEnds > 0 , "name" ])
			if(verbose >= 2) message(length(exclude), " partially matching probes filtered out")
			if(length(exclude) > 0) {
				hits = hits[ -which(hits$name %in% exclude) ,]
			}
		}
		
		# Ordering
		hits = hits[ order(hits$chrom, hits$start, hits$name) ,]
		row.names(hits) = NULL
	}
	
	return(hits)
}
