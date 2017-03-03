# Produces a cghRA.design object from an Agilent TDT file
# Author : Sylvain Mareschal <mareschal@ovsa.fr>

Agilent.design = function(
		file,
		name = NULL,
		organism = as.character(NA),
		assembly = as.character(NA),
		chromosomes	= NULL,
		...
	) {
	# File parsing
	dat <- utils::read.table(
		file = file,
		sep = "\t",
		dec = ".",
		header = TRUE,
		comment.char = "",
		quote = NULL,
		stringsAsFactors = FALSE
	)
	
	# Coordinates
	pattern <- "^chr([0-9XxYy]+):([0-9]+)\\-([0-9]+)$"
	nameCol <- dat$ID
	nameCol[ grepl("^NA\\.[0-9]+$", nameCol) ] <- NA
	chrom <- strand <- start <- end <- rep(NA, nrow(dat))
	valid <- grepl(dat$ChromosomalLocation, pattern=pattern)
	chrom[valid] <- gsub(pattern, "\\1", dat$ChromosomalLocation[valid])
	start[valid] <- as.integer(gsub(pattern, "\\2", dat$ChromosomalLocation[valid]))
	end[valid] <- as.integer(gsub(pattern, "\\3", dat$ChromosomalLocation[valid]))
	strand[valid] <- "+"
	
	# Factorial chrom
	if(is.null(chromosomes)) {
		chrom <- factor(chrom)
	} else {
		chrom <- factor(chrom, levels=chromosomes)
	}
	
	# Factorial strand
	strand <- factor(strand, levels=c("-","+"))
	
	# Default design name
	if(is.null(name)) {
		name <- sprintf("Agilent %s x %s", max(dat$Row), max(dat$Column))
	}
	
	# New object
	design <- cghRA.design(
		name = dat$ID,
		chrom = chrom,
		strand = strand,
		start = start,
		end = end,
		id = dat$RefNumber,
		control = dat$ControlType,
		row = dat$Row,
		col = dat$Column,
		quality = dat$PerformanceScore,
		.name = name,
		.organism = organism,
		.assembly = assembly
	)
	
	return(design)
}
