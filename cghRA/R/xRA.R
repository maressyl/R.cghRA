# Short/Long Recurrent Abnormality computing from CGH penetrance (Lenz et al. PNAS 2008, PMID:18765795)
# Author : Sylvain Mareschal <mareschal@ovsa.fr>

SRA = function(...) {
	return(
		xRA(
			...,
			lengthMax = 25e6,
			lengthMin = NA,
			gaps.width = 500e3,
			gaps.ratio = 1
		)
	)
}

LRA = function(...) {
	return(
		xRA(
			...,
			lengthMax = NA,
			lengthMin = 15e6,
			gaps.width = 10e6,
			gaps.ratio = 1.5
		)
	)
}

xRA = function(
		segTables,
		value = "copies",
		states = list(deletion=c(-Inf, -0.5), gain=c(0.5, Inf)),
		sampleMin = 2,
		quiet = FALSE,
		lengthMax,
		lengthMin,
		gaps.width,
		gaps.ratio
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
		groupList <- names(segTables)
		for(i in 1:length(segTables)) {
			segTables[[i]]$sample <- i
			allSeg <- rbind(allSeg, segTables[[i]])
		}
	} else {
		stop("'segTables' must be a list of data.frames or a single data.frame")
	}
	
	# Sample min
	if(sampleMin < 1) {
		sampleMin <- sampleMin*length(groupList)
	}
	
	# Segment lengths
	allSeg$length <- allSeg$end - allSeg$start
	
	# For each state
	output <- list()
	for(stateName in names(states)) {
		if(!isTRUE(quiet)) message(sprintf("Processing state \"%s\" :", stateName))
		
		# Get a fresh segment list
		tab <- allSeg
		
		# Filter on current state
		if(!isTRUE(quiet)) message("   Filtering ... ", appendLF=FALSE)
		if(length(states[[stateName]]) == 2)        { filter <- tab[,value] >= states[[stateName]][1] & tab[,value] < states[[stateName]][2]
		} else if(length(states[[stateName]]) == 1) { filter <- tab[,value] == states[[stateName]]
		} else                                      { stop("Each element in 'states' must contain two boundaries or a single value")
		}
		tab <- tab[ !is.na(filter) & filter ,]
		if(!isTRUE(quiet)) message("done")
		
		# Fast exit without suitable segment
		if(nrow(tab) == 0) {
			output[[ stateName ]] <- data.frame(
				chrom = allSeg$chrom[0],
				value = allSeg[[value]][0],
				overlap.start = integer(0),
				overlap.end = integer(0),
				start = integer(0),
				end = integer(0),
				extended.start = integer(0),
				extended.end = integer(0),
				stringsAsFactors = FALSE
			)
			next
		}
		
		# Neighbor merging
		if(!is.na(gaps.width) && !is.na(gaps.ratio)) {
			r <- 2L
			if(!isTRUE(quiet)) message("   Merging ... ", appendLF=FALSE)
			while(r <= nrow(tab)) {
				if(tab[r-1,"sample"] == tab[r,"sample"] && tab[r-1,"chrom"] == tab[r,"chrom"]) {
					gapLength <- tab[r,"start"] - tab[r-1,"end"]
					if(gapLength < gaps.width && gapLength < tab[r-1,"length"]*gaps.ratio && gapLength < tab[r,"length"]*gaps.ratio) {
						# Merge
						tab[r-1,"end"] <- tab[r,"end"]
						tab[r-1,"length"] <- tab[r-1,"end"] - tab[r-1,"start"]
						tab[r-1,"copies"] <- NA
						tab <- tab[-r,]
					} else {
						r <- r + 1L
					}
				} else {
					r <- r + 1L
				}
			}
			if(!isTRUE(quiet)) message("done")
		}
		
		# Segment size filtering
		if(!is.na(lengthMax)) {
			tab <- tab[ tab$length <= lengthMax ,]
		}
		if(!is.na(lengthMin)) {
			tab <- tab[ tab$length >= lengthMin ,]
		}
		
		# Skip copies
		tab$segment <- 1:nrow(tab)
		rownames(tab) <- NULL
		
		# Parallelization
		if(!isTRUE(quiet)) message("   Parallelization ... ", appendLF=FALSE)
		para <- parallelize(tab, "segment", quiet=TRUE)
		paraMap <- para[,1:3]
		paraSeg <- para[,4:ncol(para)]
		paraPen <- rowSums(!is.na(paraSeg))
		if(!isTRUE(quiet)) message("done")
		
		# Loop over the genome
		RA <- NULL
		if(!isTRUE(quiet)) message("   Identifying RA : ", appendLF=FALSE)
		repeat{
			# New location
			maxIndex <- which.max(paraPen)
			maxValue <- paraPen[maxIndex]
			if(maxValue < sampleMin) break
			
			# Indexes in tab of segments in the "overlapping group"
			overlappingGroup <- as.integer(paraSeg[maxIndex,])
			overlappingGroup <- overlappingGroup[ !is.na(overlappingGroup) ]
			
			# Penetrance with only "overlapping group" segments
			ogPara <- parallelize(tab[overlappingGroup,], "segment", quiet=TRUE)
			ogPene <- rowSums(!is.na(ogPara[4:ncol(ogPara)]))
			
			# Extended regions
			sra23 <- ogPara[ ogPene > maxValue * 2/3 ,]
			sra13 <- ogPara[ ogPene > maxValue * 1/3 ,]
			
			# Add to results
			RA <- rbind(
				RA,
				data.frame(
					chrom = paraMap[maxIndex,"chrom"],
					value = maxValue/length(groupList),
					overlap.start = paraMap[maxIndex,"start"],
					overlap.end = paraMap[maxIndex,"end"],
					start = min(sra23$start),
					end = max(sra23$end),
					extended.start = min(sra13$start),
					extended.end = max(sra13$end),
					stringsAsFactors = FALSE
				)
			)
			
			# Mask in penetrance
			paraPen[ paraMap$chrom == paraMap[maxIndex,"chrom"] & paraMap$start >= min(sra13$start) & paraMap$end <= max(sra13$end) ] <- NA
			if(!isTRUE(quiet)) message(sum(!is.na(paraPen) & paraPen >= sampleMin), ", ", appendLF=FALSE)
		}
		if(!isTRUE(quiet)) message("done")
		
		# Finitions
		rownames(RA) <- NULL
		output[[ stateName ]] <- RA
	}
	
	return(output)
}
