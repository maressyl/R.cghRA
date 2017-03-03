# Copy number modelization based on recurring log-ratios across segments
# Author : Sylvain Mareschal <mareschal@ovsa.fr>

model.auto = function(
		segLogRatios,
		segChroms,
		segLengths = rep(1L, length(segLogRatios)),
		from = 0.02,
		to = 0.5,
		by = 0.001,
		precision = 512,
		maxPeaks = 8,
		minWidth = 0.15,
		maxWidth = 0.9,
		minDensity = 0.001,
		peakFrom = -2,
		peakTo = 1.3,
		ploidy = 0,
		discreet = FALSE,
		method = c("stm", "sdd", "ptm"),
		exclude = c("X", "Y", "Xp", "Xq", "Yp", "Yq")
		)
	{
	# Length checks
	if(length(segLogRatios) == 0 || !is.numeric(segLogRatios))                stop("'segLogRatios' must be a non empty numeric vector")
	if(length(segLengths) != length(segLogRatios) || !is.integer(segLengths)) stop("'segLengths' must be a non empty integer vector")
	
	# Arg check
	method <- match.arg(method)
	minPeaks <- ifelse(method == "sdd", 3, 2)
	
	# Log-ratio related Copy Numbers (LCN)
	segLCN <- LCN(segLogRatios, exact=TRUE)
	segWeights <- segLengths[ ! segChroms %in% exclude ]
	segWeights <- segWeights / sum(segWeights)
	
	# Bandwidth testing
	scores <- matrix(
		data = as.numeric(NA),
		ncol = 10L,
		nrow = 1L + as.integer(ceiling((to - from) / by)),
		dimnames = list(
			NULL,
			c("bw", "peaks", "peakFrom", "peakTo", "center", "width", "ploidy", "sdd", "ptm", "stm")
		)
	)
	b <- 1
	bw <- from
	repeat {
		# Test bandwidth
		testDensity <- stats::density(
			x = segLCN[ ! segChroms %in% exclude ],
			bw = bw,
			kernel = "gaussian",
			weights = segWeights,
			n = precision
		)
		
		# Local maxima
		y <- testDensity$y / max(testDensity$y)
		i <- 2:(length(y)-1)
		iMax <- 1 + which(y[i-1] <= y[i] & y[i+1] <= y[i] & y[i+1] != y[i-1] & y[i] > minDensity)
		xMax <- testDensity$x[ iMax ]
		yMax <- testDensity$y[ iMax ]
		
		# Peak filtering
		if(!is.na(peakFrom)) {
			fromLCN <- LCN(peakFrom, exact=TRUE)
			iMax <- iMax[xMax >= fromLCN]
			xMax <- xMax[xMax >= fromLCN]
			yMax <- yMax[xMax >= fromLCN]
		}
		if(!is.na(peakTo)) {
			toLCN <- LCN(peakTo, exact=TRUE)
			iMax <- iMax[xMax <= toLCN]
			xMax <- xMax[xMax <= toLCN]
			yMax <- yMax[xMax <= toLCN]
		}
		
		# Bandwidth
		scores[b, "bw"] <- bw
		scores[b, "peaks"] <- length(iMax)
		scores[b, "peakFrom"] <- peakFrom
		scores[b, "peakTo"] <- peakTo
		
		if(scores[b, "peaks"] >= minPeaks) {
			# Model
			scores[b, "center"] <- xMax[ which.max(yMax) ]
			scores[b, "width"] <- stats::median(diff(xMax))
			scores[b, "ploidy"] <- ploidy
			
			# Standard Deviation of Diffs
			scores[b, "sdd"] <- stats::sd(diff(xMax))
			
			# Peaks to Model
			ptm <- copies(xMax, from="LCN", center=scores[b, "center"], width=scores[b, "width"], ploidy=scores[b, "ploidy"], exact=TRUE)
			ptm <- mean(abs(ptm - round(ptm)))
			scores[b, "ptm"] <- ptm
			
			# Segments to Model
			stm <- copies(segLCN[ ! segChroms %in% exclude ], from="LCN", center=scores[b, "center"], width=scores[b, "width"], ploidy=scores[b, "ploidy"], exact=TRUE)
			stm <- stats::weighted.mean(abs(stm - round(stm)), segLengths[ ! segChroms %in% exclude ])
			scores[b, "stm"] <- stm
		}
		
		# Repeat until bandwidth is too large
		bw = bw + by
		if(bw > to) {  ### if(scores[b, "peaks"] < minPeaks || bw > to) {
			break
		} else {
			b <- b +1
		}
	}
	
	# Selection
	scores <- scores[ !is.na(scores[,method]) & scores[,"peaks"] < maxPeaks & scores[,"width"] >= minWidth & scores[,"width"] <= maxWidth , , drop=FALSE ]
	if (nrow(scores) > 0) {
		return(scores[ which.min(scores[,method]) , , drop=TRUE ])
	} else {
		if(isTRUE(discreet)) {
			output <- rbind(scores[0,], NA)[,,drop=TRUE]
			output["peakFrom"] <- peakFrom
			output["peakTo"] <- peakTo
			return(output)
		} else {
			stop("No satisfying 'bandwidth' provided")
		}
	}
}
