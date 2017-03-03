# Produces a cghRA.design object from a custom TSV file
# Author : Sylvain Mareschal <mareschal@ovsa.fr>

custom.design = function(
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
		quote = "\"",
		na.strings = c("NA", ""),
		stringsAsFactors = FALSE
	)
	
	# Mandatory columns
	missing <- setdiff(c("chrom", "start", "end"), colnames(dat))
	if(length(missing) > 0L) stop("Mandatory columns are missing : ", paste(missing, collapse=", "))
	
	# Factorial strand
	if("strand" %in% colnames(dat)) { dat$strand <- factor(dat$strand, levels=c("-","+"))
	} else                          { dat$strand <- factor(NA, levels=c("-","+"))
	}
	
	# Use provided chromosome levels
    if(is.null(chromosomes)) { dat$chrom <- factor(dat$chrom)
    } else                   { dat$chrom <- factor(dat$chrom, levels=chromosomes)
    }
    
	# Merging ID
	if(!"id" %in% colnames(dat)) {
		dat$id <- 1:nrow(dat)
		warning("No 'id' column provided, assuming that design and probe files will be ordered in the same way")
	}
	
	# Without name, use ID as name
	if(!"name" %in% colnames(dat)) {
		dat$name <- as.character(dat$id)
		warning("No 'name' column provided (not recommended), using 'id' as probe name")
	}
	
	# Use file name as default
	if(is.null(name)) {
		name <- basename(file)
	}
	
	# New object
	design <- cghRA.design(
		dat,
		.name = name,
		.organism = organism,
		.assembly = assembly
	)
	
	return(design)
}
