# Computes the proportion of studied samples with a specific copy state for each position of the genome
# Author : Sylvain Mareschal <mareschal@ovsa.fr>

penetrance = function(
		segParallel,
		states = list(deletion=c(-Inf, -0.5), gain=c(0.5, Inf)),
		na = c("fill", "keep", "false"),
		mergeOnValue = FALSE,
		bool = FALSE,
		quiet = FALSE
		)
	{
	# Arg checks
	na <- match.arg(na)
	if(!is.data.frame(segParallel))                                     stop("'segParallel' must be a data.frame")
	if(!identical(names(segParallel)[1:3], c("chrom", "start", "end"))) stop("'segParallel' first columns must be 'chrom', 'start' and 'end'")
	if(ncol(segParallel) <= 3)                                          stop("'segParallel' must have a distinct column for each sample")
	
	# State assignation
	output = list()
	for(stateName in names(states)) {
		if(!isTRUE(quiet)) message("Computing penetrance for state ", stateName, ": ", appendLF=FALSE)
		
		if(!isTRUE(quiet)) message("calling, ", appendLF=FALSE)
		if(length(states[[stateName]]) == 2) {
			# Interval
			isInState <- segParallel[, 4:ncol(segParallel)] >= states[[stateName]][1] & segParallel[, 4:ncol(segParallel)] < states[[stateName]][2]
		} else if(length(states[[stateName]]) == 1) {
			# Single value
			isInState <- segParallel[, 4:ncol(segParallel)] == states[[stateName]]
		} else {
			stop("Each element in 'states' must contain two boundaries or a single value")
		}
		
		# Fill NA (gaps and above-limit segments)
		if(na == "fill") {
			if(!isTRUE(quiet)) message("NA filling, ", appendLF=FALSE)
			
			# Chromosome boundaries
			startingLines <- sort(tapply(X=1:nrow(segParallel), INDEX=segParallel$chrom, FUN=min))
			endingLines <- c(utils::tail(startingLines - 1L, -1L), nrow(segParallel))
						
			# Mask full NA chromosomes
			unMasked <- matrix(TRUE, ncol=ncol(isInState), nrow=nrow(isInState), dimnames=dimnames(isInState))
			for(i in 1:length(startingLines)) {
				emptySamples <- which(apply(is.na(isInState[ startingLines[i] : endingLines[i] , , drop=FALSE ]), 2, all))
				unMasked[ startingLines[i] : endingLines[i] , emptySamples ] <- FALSE
			}
			
			# Chromosome starts
			for(l in startingLines) {
				samplesToFill <- which(is.na(isInState[l,]) & unMasked[l,])
				for(k in samplesToFill) {
					to <- l + 1L
					while(is.na(isInState[to,k])) to <- to + 1L
					isInState[ l:to , k ] <- isInState[ to , k ]
				}
			}
			
			# Chromosome ends
			for(l in endingLines) {
				samplesToFill <- which(is.na(isInState[l,]) & unMasked[l,])
				for(k in samplesToFill) {
					to <- l - 1L
					while(is.na(isInState[to,k])) to <- to - 1L
					isInState[ l:to , k ] <- isInState[ to , k ]
				}
			}
			
			# Internal gaps
			beforeGap <- afterGap <- which(is.na(isInState) & unMasked)
			while(any(lgc <- is.na(isInState[ beforeGap ]))) beforeGap[lgc] <- beforeGap[lgc] - 1L
			while(any(lgc <- is.na(isInState[ afterGap ])))  afterGap[lgc] <- afterGap[lgc] + 1L
			gaps <- cbind(before=beforeGap, after=afterGap)
			gaps <- unique(gaps)
			
			# Gaps to fill (NA between identical states, not necessary identical values)
			gaps <- gaps[ isInState[ gaps[,"before"] ] == isInState[ gaps[,"after"] ] ,]
			
			# Gap filling
			if(nrow(gaps) > 0) for(i in 1:nrow(gaps)) isInState[ (gaps[i,"before"] + 1L) : (gaps[i,"after"] - 1L) ] <- isInState[ gaps[i,"before"] ]
		} else if(na == "false") {
			if(!isTRUE(quiet)) message("NA filling, ", appendLF=FALSE)
			
			# NA always considered as "not in state"
			isInState[ is.na(isInState) ] <- FALSE
		}
		
		if(isTRUE(bool)) {
			# Logical table
			pen <- data.frame(
				segParallel[, 1:3],
				isInState,
				stringsAsFactors = FALSE,
				check.names = FALSE
			)
		} else {
			# Penetrance computation
			if(na == "false") { value <- rowSums(isInState) / ncol(isInState)
			} else            { value <- rowSums(isInState, na.rm=TRUE) / rowSums(!is.na(isInState))
			}
			
			# Penetrance table
			pen <- data.frame(
				segParallel[, 1:3],
				value = value,
				stringsAsFactors = FALSE
			)
			
			# Add column to distinct same penetrance with different samples
			if(!isTRUE(mergeOnValue)) {
				if(!isTRUE(quiet)) message("pre-merging, ", appendLF=FALSE)
				samples <- 1:ncol(isInState)
				samplesIn <- function(lgc) {
					paste(which(lgc), collapse=",")
				}
				pen$samplesIn <- apply(isInState, 1, samplesIn)
			}
			
			# Similar region grouping
			if(!isTRUE(quiet)) message("merging, ", appendLF=FALSE)
			pen <- segMerge(pen)
			
			# Remove dummy column
			if(!isTRUE(mergeOnValue)) pen$samplesIn <- NULL
		}
		
		# Add to output list
		output[[stateName]] <- pen
		
		if(!isTRUE(quiet)) message("done")
	}
	
	return(output)
}
