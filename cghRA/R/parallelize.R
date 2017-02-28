# Summarize log-ratios or copy-numbers of segments in a data.frame with a column for each sample
# Author : Sylvain Mareschal <mareschal@ovsa.fr>

parallelize = function(
		segTables,
		value = "logRatio",
		digits = 3,
		quiet = FALSE,
		chroms = NULL
		)
	{
	# Arg checks
	if(is.data.frame(segTables)) {
		# Checks
		if(any(!c("sample", "chrom", "start", "end", value) %in% names(segTables))) stop(sprintf("'segTables' must have 'sample', 'chrom', 'start', 'end' and [value] columns", i))
		
		# Sample merging
		allSeg <- segTables
		groupList <- unique(allSeg$sample)
	} else if(is.list(segTables)) {
		# Checks
		for(i in 1:length(segTables)) {
			if(!is.data.frame(segTables[[i]]))                                     stop(sprintf("'segTables[[%i]]' must be a data.frame", i))
			if(any(!c("chrom", "start", "end", value) %in% names(segTables[[i]]))) stop(sprintf("'segTables[[%i]]' must have 'chrom', 'start', 'end' and [value] columns", i))
		}
		
		# Sample merging
		allSeg <- NULL
		if(is.null(names(segTables))) { groupList <- 1:length(segTables)
		} else                        { groupList <- names(segTables)
		}
		for(i in 1:length(segTables)) {
			segTables[[i]]$sample <- groupList[i]
			allSeg <- rbind(allSeg, segTables[[i]])
		}
	} else {
		stop("'segTables' must be a list of data.frames or a single data.frame")
	}
	
	# Rounding
	if(!is.na(digits)) {
		allSeg[,value] = round(allSeg[,value], digits=digits)
	}
	
	regions = NULL
	if(!isTRUE(quiet)) message("Parallelizing chromosome ", appendLF=FALSE)
	if(is.null(chroms)) chroms <- unique(allSeg$chrom)
	for(chrom in chroms) {
		if(!isTRUE(quiet)) message(chrom, ", ", appendLF=FALSE)
		# Merging
		regions = rbind(
			regions,
			data.frame(
				chrom = chrom,
				segPile(
					starts = allSeg[allSeg$chrom == chrom, "start"],
					ends = allSeg[allSeg$chrom == chrom, "end"],
					values = allSeg[allSeg$chrom == chrom, value],
					groups = allSeg[allSeg$chrom == chrom, "sample"],
					groupList = groupList
				),
				stringsAsFactors = FALSE,
				check.names = FALSE
			)
		)
	}
	if(!isTRUE(quiet)) message("done")
	
	return(regions)
}
