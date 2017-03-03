# Produces a cghRA.regions object from a simplified karyotype formula
# Author : Sylvain Mareschal <mareschal@ovsa.fr>

parseKaryo <- function(formula, bandTrack, name=as.character(NA), design=NULL, alteratedOnly=TRUE) {
	# Split clones
	clones <- strsplit(formula, split=" ; ")[[1]]
	
	# Initialize bands
	if(!bandTrack$isArmed()) stop("'bandTrack' is supposed to be armed")
	bands <- bandTrack$copy()
	starts <- with(bands$extract(), tapply(start, chrom, head, 1))
	ends <- with(bands$extract(), tapply(end, chrom, tail, 1))
	
	# Add chromosome boundaries to bands
	arm <- c("p", "p", "q", "q")
	boundary <- c("ter", "cen", "cen", "ter")
	coord <- list(starts, ends, starts, ends)
	for(i in 1:4) {
		x <- coord[[i]][ grep(sprintf("%s$", arm[i]), names(coord[[i]])) ]
		bands$addDataFrame(
			data.frame(
				name = sprintf("%s%s", names(x), boundary[i]),
				chrom = factor(names(x), bands$chromosomes()),
				strand = factor("+",c("-","+")),
				start = as.integer(x),
				end = as.integer(x),
				stringsAsFactors = FALSE
			)
		)
	}
	
	# Erase arms
	bands$eraseArms()
	
	genomes <- list()
	weights <- integer(0)
	normals <- integer(0)
	ploidies <- integer(0)
	for(clone in clones) {
		# Meta extraction
		regex <- "^([0-9]+)(?:-([0-9]+))?<([0-9]+)n>((?: *, *[0-9]+\\([0-9pqctXY-]+\\))*)(?: *\\[(cp)?([0-9]+)\\])?$"
		if(!grepl(regex, clone)) stop("Meta parse error in clone \"", clone, "\"")
		minChrom <- as.integer(sub(regex, "\\1", clone))
		maxChrom <- as.integer(sub(regex, "\\2", clone))
		ploidy <- as.integer(sub(regex, "\\3", clone))
		cnv <- sub(regex, "\\4", clone)
		cp <- sub(regex, "\\5", clone) == "cp"
		nmito <- as.integer(sub(regex, "\\6", clone))
		
		# Default profile (ploidy)
		sizes <- with(bands$extract(), tapply(end, chrom, max))
		profile <- data.frame(
			chrom = names(sizes),
			start = rep(1L, length(sizes)),
			end = as.integer(sizes),
			copies = as.numeric(ploidy),
			stringsAsFactors = FALSE
		)
		
		# Local alteration collection
		locals <- NULL
		regex1 <- "^([0-9XY]+)\\((([0-9XY]+)(?:[pq][0-9\\.]*)?)\\)$"
		regex2 <- "^([0-9XY]+)\\(([0-9XY]+)([pq][0-9\\.ct]+)-([pq][0-9\\.ct]+)\\)$"
		elements <- strsplit(cnv, split=" *, *")[[1]][-1]
		for(element in elements) {
			# Data extraction
			if(grepl(regex1, element)) {
				# Single band or chromosome
				copies <- as.numeric(sub(regex1, "\\1", element))
				from <- to <- sub(regex1, "\\2", element)
				chrom <- sub(regex1, "\\3", element)
			} else if(grepl(regex2, element)) {
				# Localized band(s)
				copies <- as.numeric(sub(regex2, "\\1", element))
				chrom <- sub(regex2, "\\2", element)
				from <- sub(regex2, sprintf("%s\\3", chrom), element)
				to <- sub(regex2, sprintf("%s\\4", chrom), element)
			} else stop("Parse error in \"", element, "\"")
			
			# Band selection expression
			patterns <- sprintf("^%s", c(from, to))
			patterns[ !grepl("[pq]", patterns) ] <- sprintf("%s[pq]", patterns[ !grepl("[pq]", patterns) ])
			expr <- parse(text=sprintf("grepl(\"%s\", name) | grepl(\"%s\", name)", patterns[1], patterns[2]))
			
			# Selection coordinates
			region <- bands$extract(expr)
			start <- min(region$start)
			end <- max(region$end)
			
			# Add to local CNV table
			locals <- rbind(locals, data.frame(chrom=chrom, start=start, end=end, copies=copies, stringsAsFactors=FALSE))
		}
		
		if(is.null(locals)) {
			# No local alteration
			para <- profile
			
			# Normal formula, record for later filtering
			if(ploidy == 2L) normals <- c(normals, length(genomes)+1L)
		} else {
			# Merge regions
			para <- parallelize(list(profile=profile, locals=locals), value="copies")
	
			# Merge copies
			copies <- para$locals
			copies[ is.na(copies) ] <- para[ is.na(copies) , "profile" ]
			para$copies <- copies
			para$profile <- para$locals <- NULL
		}
		
		# Add to clone list
		genomes[[ length(genomes)+1 ]] <- para
		weights[ length(weights)+1 ] <- nmito
		ploidies[ length(ploidies)+1 ] <- ploidy
	}
	
	# Normal clone filtering
	if(isTRUE(alteratedOnly)) {
		if(length(normals) > 0 && length(normals) < length(genomes)) {
			# There are some normal clones to filter out, and it will leave some alterated ones
			genomes <- genomes[ -normals ]
			weights <- weights[ -normals ]
			ploidies <- ploidies[ -normals ]
		} # else there is no normal clone, or there are only normal clones that can be merged safely
	}
	
	if(length(genomes) > 1) {
		# Merge clones
		para <- parallelize(genomes, value="copies")
		para$copies <- apply(para[,-(1:3)], 1, stats::weighted.mean, w=weights)
		para <- para[,c("chrom","start","end","copies")]
	} else {
		# Single clone
		para <- genomes[[1]]
	}
	
	# Convert to track.table
	para$strand <- factor("+", levels=c("-","+"))
	para$chrom <- factor(para$chrom, levels=bands$chromosomes())
	para <- track.table(para, .name=name, .organism=bands$organism, .assembly=bands$assembly, .makeNames=TRUE)
	
	# cghRA.regions conversion
	if(is(design, "cghRA.design")) {
		# Add cghRA.copies columns
		para$cross(design$eraseArms(TRUE), colname="probes", type="count")
		para$addColumn(rep(as.double(NA), para$rowCount), "logRatio")

		# Convert to cghRA.copies
		object <- new("cghRA.copies")
		object$import(para)

		# Switch return
		para <- object
	}
	
	# Clone summary
	clones <- weights
	names(clones) <- sprintf("%in", ploidies)
	
	return(list(clones=clones, copies=para))
}

