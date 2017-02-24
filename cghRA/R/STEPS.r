# "Selective Trends Evidenced by Penetrance Surges" computing from CGH penetrance
# Author : Sylvain Mareschal <mareschal@ovsa.fr>

STEPS <- function(
		segPenetrance,
		dpen = 2,
		vpen = 0.8,
		gpen = 0.3,
		threshold = NA,
		nested = c("merge", "flag", "none"),
		digits = 3,
		chromEnd = FALSE,
		quiet = FALSE
		)
	{
	# Arg checks
	nested <- match.arg(nested)
	if(!is.data.frame(segPenetrance)) stop("'segPenetrance' must be a data.frame")
	if(any(! c("chrom", "start", "end", "value") %in% names(segPenetrance))) {
		stop("'segPenetrance' columns must have 'chrom', 'start', 'end' and 'value' columns")
	}
	
	# Ordering
	segPenetrance <- segPenetrance[ order(segPenetrance$chrom, segPenetrance$start) ,]
	row.names(segPenetrance) <- NULL
	
	# Storage
	str.leftScore  <- rep(as.double(NA), nrow(segPenetrance))
	str.rightScore <- rep(as.double(NA), nrow(segPenetrance))
	str.leftIndex  <- rep(as.integer(NA), nrow(segPenetrance))
	str.rightIndex <- rep(as.integer(NA), nrow(segPenetrance))
	str.score      <- rep(as.double(NA), nrow(segPenetrance))
	
	if(!isTRUE(quiet)) message("Processing chromosome ", appendLF=FALSE)
	
	# Score computation
	values <- segPenetrance$value * 100
	widths <- segPenetrance$end - segPenetrance$start
	ranges <- split(1:nrow(segPenetrance), segPenetrance$chrom)
	for(chrom in names(ranges)) {
		if(!isTRUE(quiet)) message(chrom, ", ", appendLF=FALSE)
		k <- ranges[[chrom]]
		
		# Compute score for each penetrance segment (excluding boundaries)
		if(length(k) >= 3) for(mcr in k[ 2:(length(k)-1) ]) {
			# Unilateral score storage
			earlyBreak <- FALSE
			uniScore <- rep(as.double(NA), nrow(segPenetrance))
			uniScore[mcr] <- 0
			
			# Right score computation
			first <- TRUE
			for(i in (mcr+1L):k[length(k)]) {
				dval <- values[i-1L] - values[i]
				if(dval > 0) {
					# Going down, increase score
					uniScore[i] <- uniScore[i-1L] + dval
				} else if(first) {
					# First move is up, not a MCR
					earlyBreak <- TRUE
					break
				} else {
					# Going up, penalize
					uniScore[i] <- uniScore[i-1L] + dval * dpen
				}
				first <- FALSE
			}
			if(earlyBreak) next
			
			# Left score computation
			first <- TRUE
			for(i in (mcr-1L):k[1L]) {
				dval <- values[i+1L] - values[i]
				if(dval > 0)     {
					# Going down, increase score
					uniScore[i] <- uniScore[i+1L] + dval
				} else if(first) {
					# First move is up, not a MCR
					earlyBreak <- TRUE
					break
				} else {
					# Going up, penalize
					uniScore[i] <- uniScore[i+1L] + dval * dpen
				}
				first <- FALSE
			}
			if(earlyBreak) next
			
			# Find best combination
			uniScore[ uniScore < 0 ] <- 0
			starts <- segPenetrance$start
			ends <- segPenetrance$end
			
			### DEBUG : Store score
			### mtx <- matrix(as.double(NA), nrow=length(uniScore), ncol=length(uniScore))
			
			# Consider any block end on the left/right of the MCR
			leftCandidates  <- (mcr-1L):k[1L]
			rightCandidates <- (mcr+1L):k[length(k)]
			
			# Retain only uniScore local maxima (penetrance local minima)
			leftCandidates  <- leftCandidates[ localMax(uniScore)[ leftCandidates ] ]
			rightCandidates <- rightCandidates[ localMax(uniScore)[ rightCandidates ] ]
			leftCandidates <- leftCandidates[ !is.na(leftCandidates) ]
			rightCandidates <- rightCandidates[ !is.na(rightCandidates) ]
			
			# Ignore local maximum filter for chromosome boundaries
			leftCandidates  <- c(leftCandidates, k[1L])
			rightCandidates <- c(rightCandidates, k[length(k)])
			
			# Retain only candidates with lower penetrance than the MCR
			leftCandidates  <- leftCandidates[ values[leftCandidates] < values[mcr] ]
			rightCandidates  <- rightCandidates[ values[rightCandidates] < values[mcr] ]
			
			# Assess candidates
			smax <- 0
			lmax <- rmax <- as.numeric(NA)
			for(right in rightCandidates) {
				for(left in leftCandidates) {
					score <- finalScore(
						leftScore = uniScore[left],
						rightScore = uniScore[right],
						leftWidth = starts[mcr] - starts[left],
						rightWidth = ends[right] - ends[mcr],
						vpen=vpen, gpen=gpen
					)
					
					### DEBUG : Store score
					### mtx[left,right] <- score
					
					# Retain best score encountered
					if(score > smax || (score == 0 && is.na(lmax))) {
						lmax <- left
						rmax <- right
						smax <- score
					}
				}
			}
			
			### DEBUG : Plot start / end combinations
			### layout(matrix(c(3,4,1,2), ncol=2), widths=c(1, 2), heights=c(2, 1))
			### 
			### par(mar=c(1,1,1,1))
			### image(mtx, axes=FALSE)
			### 
			### par(mar=c(4,1,1,1))
			### plot(x=1:length(uniScore), y=uniScore, xlab="Left end", ylab="Score")
			### abline(v=lmax)
			### 
			### par(mar=c(1,4,1,1))
			### plot(y=1:length(uniScore), x=uniScore, ylab="Right end", xlab="Score")
			### abline(h=rmax)
			### 
			### plot(x=NA, y=NA, xlim=0:1, ylim=0:1, xlab="", ylab="", xaxt="n", yaxt="n", bty="n")
			
			# Store result
			str.leftScore[mcr] <- uniScore[lmax]
			str.rightScore[mcr] <- uniScore[rmax]
			str.leftIndex[mcr] <- lmax
			str.rightIndex[mcr] <- rmax
			str.score[mcr] <- smax
		}
	}
	
	if(!isTRUE(quiet)) message("done")
	if(!isTRUE(quiet)) message("Post-processing ... ", appendLF=FALSE)
	
	# Data merging
	segPenetrance$score <- str.score
	segPenetrance$start.baseline <- segPenetrance[ str.leftIndex , "start" ]
	segPenetrance$leftScore <- str.leftScore
	segPenetrance$end.baseline <- segPenetrance[ str.rightIndex , "end" ]
	segPenetrance$rightScore <- str.rightScore
	segPenetrance$index <- 1:nrow(segPenetrance)
	
	# Filter out regions with higher neighbors
	mcrTable <- segPenetrance[ !is.na(segPenetrance$score) ,]
	
	# Score filtering
	if(!is.na(threshold)) mcrTable <- mcrTable[ mcrTable$score >= threshold ,]
	
	# Explicite start / end
	colnames(mcrTable)[ colnames(mcrTable) == "start" ] <- "start.MCR"
	colnames(mcrTable)[ colnames(mcrTable) == "end" ] <- "end.MCR"
	colnames(mcrTable)[ colnames(mcrTable) == "value" ] <- "penetrance.MCR"
	
	# Nested filter
	if(nested != "none" && nrow(mcrTable) > 1L) {
		# Allocate nests
		n <- 1L
		nest <- 1L
		for(i in 2:nrow(mcrTable)) {
			if(mcrTable[i,"chrom"] != mcrTable[i-1,"chrom"] || (mcrTable[i,"start.baseline"] >= mcrTable[i-1,"end.MCR"] && mcrTable[i-1,"end.baseline"] <= mcrTable[i,"start.MCR"])) n <- n + 1L
			nest[i] <- n
		}
		
		if(nested == "flag") {
			# Return nest ID
			mcrTable$nest <- nest
		} else if(nested == "merge") {
			# Keep the best from each nest		
			nestMembers <- tapply(X=rownames(mcrTable), INDEX=nest, FUN=c)
			nestBest <- tapply(X=mcrTable$score, INDEX=nest, FUN=which.max)
			sel <- mapply(FUN="[", nestMembers, nestBest)
		
			# Apply filtering
			mcrTable <- mcrTable[ sel ,]
		}
	}
	
	# Ordering
	mcrTable <- mcrTable[ order(mcrTable$score, decreasing=TRUE) ,]
	row.names(mcrTable) <- NULL
	
	if(!isTRUE(quiet)) message("done")
	
	return(mcrTable)
}


# Returns a logical vector, identifying local maxima
# Author : sylvain.mareschal@etu.univ-rouen.fr
localMax <- function(x) {
	n <- length(x)
	if(n > 2) { return(c(FALSE, x[1:(n-2)] <= x[2:(n-1)] & x[2:(n-1)] >= x[3:n], FALSE))
	} else    { return(rep(FALSE, n))
	}
}


# Combine left and right scores into a final MCR score
# Author : sylvain.mareschal@etu.univ-rouen.fr
finalScore <- compiler::cmpfun(
	function(leftScore, rightScore, leftWidth, rightWidth, vpen, gpen) {
		# Base score
		score <- leftScore + rightScore
		
		# Weight for vertical equilibrium (ratio of penetrance values at ends)
		if(leftScore < rightScore) { score <- score * ((1+leftScore) / (1+rightScore))^vpen
		} else                     { score <- score * ((1+rightScore) / (1+leftScore))^vpen
		}
		
		# Weight for horizontal equilibrium (ratio of genomic widths)
		if(leftWidth < rightWidth) { score <- score * ((1+leftWidth) / (1+rightWidth))^gpen
		} else                     { score <- score * ((1+rightWidth) / (1+leftWidth))^gpen
		}
		
		return(score)
	},
	options = list("optimize"=3)
)

