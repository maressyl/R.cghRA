# Various genomic context based values computed on probes
# Author : Sylvain Mareschal <mareschal@ovsa.fr>

bias <- function(
		chromFiles,
		probeChrom,
		probeStarts,
		probeEnds,
		chromPattern = "^(.+)\\.[^\\.]+$",
		fragSites = c(AluI="AG|CT", RsaI="GT|AC"),
		digits = 6,
		verbose = 1
		)
	{
	# Probe count consistency
	if(length(probeChrom) == 0)                           stop("'probeChrom' must not be empty")
	if(length(probeStarts) != length(probeChrom))         stop("'probeChrom' and 'probeStarts' must have same lengths")
	if(length(probeEnds) != length(probeChrom))           stop("'probeChrom' and 'probeEnds' must have same lengths")
	
	# fragSites check
	if(!is.character(fragSites) || any(is.na(fragSites))) stop("'fragSites' must be a non NA character vector")
	if(any(!grepl("^[ACGT]+\\|[ACGT]+$", fragSites)))     stop("Invalid 'fragSites' format")
	
	# chromFiles check
	chromFiles = as.character(chromFiles)
	for(i in 1:length(chromFiles)) {
		if(!file.exists(chromFiles[i]))                   stop(sprintf("'chromFiles[%d]' does not exist", i))
	}
	
	# Output allocation
	output = matrix(
		data = as.double(NA),
		ncol = 5,
		nrow = length(probeChrom),
		dimnames = list(
			NULL,
			c("wGC150", "wGC500", "wGCprobe", "wGCfrag", "wFragSize")
		)
	)
	
	if(verbose == 1) message("Processing chromosome ", appendLF=FALSE)
	
	for(chromFile in chromFiles) {
		# Chromosome name
		chromName = sub(chromPattern, "\\1", basename(chromFile))
		if(verbose == 1) message(chromName, appendLF=FALSE)
		
		# Probes mapped on the current chromosome
		selection = (!is.na(probeChrom) & probeChrom == chromName)
		probeCount = sum(selection)
		if(probeCount > 0) {
			# Output allocation
			wGC150 <- double(probeCount)
			wGC500 <- double(probeCount)
			wGCprobe <- double(probeCount)
			wGCfrag <- double(probeCount)
			wFragSize <- integer(probeCount)
			
			# C level processing
			results = .C("R_WACA", PACKAGE="cghRA", NAOK=FALSE, DUP=TRUE,
				chromFile,
				fragSites,
				length(fragSites),
				probeCount,
				probeStarts[selection],
				probeEnds[selection],
				wGC150,
				wGC500,
				wGCprobe,
				wGCfrag,
				wFragSize
			)
			
			# Matrix filling
			output[selection, "wGC150"] = round(results[[7]], digits)
			output[selection, "wGC500"] = round(results[[8]], digits)
			output[selection, "wGCprobe"] = round(results[[9]], digits)
			output[selection, "wGCfrag"] = round(results[[10]], digits)
			output[selection, "wFragSize"] = results[[11]]
		}
		
		if(verbose == 1) message(", ", appendLF=FALSE)
	}
	
	if(verbose == 1) message("done")
	
	return(output)
}
