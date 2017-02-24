# Compute CNV likelihood for each sample segment
# - sample.map       : sample segments mapped to the CGH-array design used, as returned by map2design()
# - dgv.map          : constitutionnal CNVs mapped to the CGH-array design used, as returned by map2design()
# - hangingThreshold : segments to score must cover at least this proportion of union(CNV, segment) for a CNV to be retained
# - minGain          : CNVs must add at least this value to the path's score to be retained
# - maxPaths         : maximal path to be computed for each segment, considering that best paths are computed first (actually true for single-element ones) and final score focus on them
# - trace            : whether to return data on path built (list with 'traces' and 'scores' components) or only the final scores (vector)
# - expand           : whether to return a vector with one row for each element of 'sample.map' or expand back to the mapped sample track
#                    (some segments may be missing or merged in the map). Notice traces are not expanded.
# Author : Sylvain Mareschal <mareschal@ovsa.fr>

cnvScore <- function(sample.map, dgv.map, hangingThreshold=0.8, minGain=0.1, maxPaths=NA, trace=FALSE, expand=TRUE, quiet=TRUE) {
	# Checks
	if(!inherits(sample.map, "segmentMap")) stop("'sample.map' must inherit 'segmentMap' class")
	if(!inherits(dgv.map, "segmentMap"))    stop("'dgv.map' must inherit 'segmentMap' class")
	
	# Final segment score storage
	out <- double(nrow(sample.map$map))
	names(out) <- rownames(sample.map$map)
	
	# Path trace
	if(isTRUE(trace)) allPaths <- NULL
	
	# For each segment (including normal ones, as can be a CNV amp in a del)
	for(i in 1:nrow(sample.map$map)) {
		if(!isTRUE(quiet) && i %% 100 == 0) message("Segment ", i, " / ", nrow(sample.map$map))
		
		# Overlapping CNVs to be used
		overlaps <- dgv.map$map[ dgv.map$map[,"end"] >= sample.map$map[i,"start"] & dgv.map$map[,"start"] <= sample.map$map[i,"end"] , , drop=FALSE ]
		
		# Filter out hanging CNVs
		candidateLength <- sample.map$map[i,"end"] - sample.map$map[i,"start"] + 1L
		unionLength <- pmax(sample.map$map[i,"end"], overlaps[,"end"]) - pmin(sample.map$map[i,"start"], overlaps[,"start"]) + 1L
		overlaps <- overlaps[ candidateLength / unionLength >= hangingThreshold , , drop=FALSE ]
		
		if(nrow(overlaps) > 0) {
			# Fast exit
			intersectionLength <- pmin(sample.map$map[i,"end"], overlaps[,"end"]) - pmax(sample.map$map[i,"start"], overlaps[,"start"]) + 1L
			unionLength <- pmax(sample.map$map[i,"end"], overlaps[,"end"]) - pmin(sample.map$map[i,"start"], overlaps[,"start"]) + 1L
			if(max(intersectionLength / unionLength) >= minGain) {
				# Path scores for the current segment
				scores <- double(0)
				
				# Build as many tilling paths as possible
				j <- 0L
				repeat {
					# Tilling path count
					j <- j + 1L
					if(!is.na(maxPaths) && j >= maxPaths) break
					
					# Reset current tilling path
					elements <- 0L                                                                   # CNV count in the current tilling path
					isCovered <- logical(sample.map$map[i,"end"] - sample.map$map[i,"start"] + 1L)   # Is each sample probe covered by the path or not
					unionStart <- sample.map$map[i,"start"]                                          # Start of union between segment and path (possibly overruning segment boundaries)
					unionEnd <- sample.map$map[i,"end"]                                              # End of union between segment and path (possibly overruning segment boundaries)
					jaccard <- 0                                                                     # Jaccard index of the current tilling path
					score <- 0                                                                       # Score of the current tilling path
					
					# Update trace
					if(isTRUE(trace)) retained <- integer(0)
					
					# Store selected overlap indexes to try path duplication (valid until some row is deleted !)
					dupIndexes <- integer(0)
					
					# Single tilling path with remaining CNVs
					while(nrow(overlaps) > 0) {
						# 'isCovered' versions for each overlap candidate
						from <- pmax(overlaps[,"start"] - sample.map$map[i,"start"] + 1L, 1L)
						to <- pmin(overlaps[,"end"] - sample.map$map[i,"start"] + 1L, length(isCovered))
						inter.new <- integer(nrow(overlaps))
						for(k in 1:nrow(overlaps)) {
							tmp <- isCovered
							tmp[ from[k]:to[k] ] <- TRUE
							inter.new[k] <- sum(tmp)
						}
						
						# Union boundaries versions for each overlap candidate
						unionStart.new <- pmin(unionStart, overlaps[,"start"])
						unionEnd.new <- pmax(unionEnd, overlaps[,"end"])
						
						# Jaccard and score versions for each overlap candidate
						jaccard.new <- inter.new / (unionEnd.new - unionStart.new + 1L)
						score.new <- jaccard.new ^ (elements + 1L)
						
						# Stop if no significant gain can be obtained with remaining CNVs
						if(max(score.new) - score < minGain) break
						
						# CNV choosen to enhance current path
						best <- which.max(jaccard.new)
						
						# Update coverage
						isCovered[ from[best]:to[best] ] <- TRUE
						unionStart <- unionStart.new[best]
						unionEnd <- unionEnd.new[best]
						jaccard <- jaccard.new[best]
						score <- score.new[best]
						
						# Add to retained CNVs
						elements <- elements + 1L
						if(isTRUE(trace)) retained <- c(retained, match(rownames(overlaps)[best], rownames(dgv.map$map)))
						
						# Remove from usable overlaps
						if(overlaps[best,"count"] > 1) {
							# Pick one from an identical segment set
							overlaps[best,"count"] <- overlaps[best,"count"] - 1L
							dupIndexes <- c(dupIndexes, best)
						} else                         {
							# Last segment with these boudaries, no possible duplication
							overlaps <- overlaps[-best,,drop=FALSE]
							dupIndexes <- c(dupIndexes, NA)
						}
					}
					
					# No more significant tilling path
					if(elements == 0L) break
					
					# Store final path score
					scores[j] <- score
					
					# Duplicate path if used segments are still available
					if(all(!is.na(dupIndexes))) {
						# Duplication amount
						pathDuplicates <- min(overlaps[ dupIndexes , "count" ])
						
						# Add paths
						scores[ j+(1:pathDuplicates) ] <- score
						j <- j + pathDuplicates
						
						# Remove used overlaps (some remain)
						overlaps[ dupIndexes , "count" ] <- overlaps[ dupIndexes , "count" ] - pathDuplicates
						overlaps <- overlaps[ overlaps[,"count"] > 0 , , drop=FALSE ]
					} else {
						# No duplication
						pathDuplicates <- 0L
					}
					
					# Store path results
					if(isTRUE(trace)) {
						allPaths <- rbind(
							allPaths,
							data.frame(                                              # One row per distinct path built for a given segment in a given sample
								seg = rownames(sample.map$map)[i],                   # List of sample row numbers for assayed segment
								seg.score = as.double(NA),                           # Final CNV score for the segment, all paths comprised
								path.count = pathDuplicates + 1L,                    # How many times the CNV path described was built
								path.jaccard = jaccard,                              # Jaccard index between the assayed segment and the CNV path described
								path.cnvCount = elements,                            # How many CNVs are included in the CNV path described
								path.cnvList = paste(retained, collapse=", "),       # Indexes in 'dgv.map' of the CNVs retained in the CNV path described
								path.score = scores[j],                              # 'path.jaccard' corrected for the amount of CNVs included in the CNV path described
								stringsAsFactors = FALSE
							)
						)
					}
				}
				
				# Combine paths into a single score
				scores <- sort(scores, decreasing=TRUE)
				out[i] <- sum(scores^(1+1:length(scores)))
				if(isTRUE(trace)) allPaths[ allPaths$seg == rownames(sample.map$map)[i] , "seg.score" ] <- out[i]
			} else {
				# No significant CNV, fast exit
				out[i] <- 0
			}
		} else {
			# No overlapping CNV, fast exit
			out[i] <- 0
		}
	}
	
	# Expand back from sample map to sample track (overlapping and invisible segments)
	if(isTRUE(expand)) {
		x <- rep(as.double(NA), sample.map$trackSize)
		ranges <- strsplit(names(out), split=",")
		indexes <- as.integer(unlist(ranges))
		values <- rep(out, sapply(ranges, length))
		x[ indexes ] <- values
		out <- x
	}
	
	if(isTRUE(trace)) { return(list(scores=out, traces=allPaths))
	} else            { return(out)
	}
}

