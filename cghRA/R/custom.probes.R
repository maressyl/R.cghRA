# Produces a cghRA.probes object from a custom TSV file
# Author : Sylvain Mareschal <mareschal@ovsa.fr>

custom.probes = function(
	file,
	columns = NULL,
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
	missing <- setdiff(c("logRatio"), colnames(dat))
	if(length(missing) > 0L) stop("Mandatory columns are missing : ", paste(missing, collapse=", "))
	
	# Merging ID
	if(!"id" %in% colnames(dat)) {
		dat$id <- 1:nrow(dat)
		warning("No 'id' column provided, assuming that design and probe files will be ordered in the same way")
	}
	
	# Rename flag columns
	if(length(columns) > 0L) {
		for(i in 1:length(columns)) colnames(dat)[ colnames(dat) == columns[i] ] <- names(columns)[i]
	}
	
	# New object
	probes <- cghRA.probes(dat, .name=basename(file))
	
	return(probes)
}
