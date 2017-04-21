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
		minPeaks = 1,
		maxPeaks = 8,
		minWidth = 0.15,
		maxWidth = 1.1,
		defWidth = 1,
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
	if(defWidth > maxWidth || defWidth < minWidth) warning("'defWidth' is outside the interval defined by 'minWidth' and 'maxWidth'")
	
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
	b <- 1L
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
		iMax <- 1L + which(y[i-1] <= y[i] & y[i+1] <= y[i] & y[i+1] != y[i-1] & y[i] > minDensity)
		xMax <- testDensity$x[ iMax ]
		yMax <- testDensity$y[ iMax ]
		
		# Peak filtering
		if(!is.na(peakFrom)) {
			fromLCN <- LCN(peakFrom, exact=TRUE)
			iMax <- iMax[xMax >= fromLCN]
			yMax <- yMax[xMax >= fromLCN]
			xMax <- xMax[xMax >= fromLCN]
		}
		if(!is.na(peakTo)) {
			toLCN <- LCN(peakTo, exact=TRUE)
			iMax <- iMax[xMax <= toLCN]
			yMax <- yMax[xMax <= toLCN]
			xMax <- xMax[xMax <= toLCN]
		}
		
		# Bandwidth
		scores[b, "bw"] <- bw
		
		# Peaks
		nPeaks <- length(iMax)
		scores[b, "peaks"] <- nPeaks
		scores[b, "peakFrom"] <- peakFrom
		scores[b, "peakTo"] <- peakTo
		
		# Model
		if(nPeaks > 0L) scores[b, "center"] <- xMax[ which.max(yMax) ]
		if(nPeaks > 1L) { scores[b, "width"] <- stats::median(diff(xMax))
		} else          { scores[b, "width"] <- defWidth
		}
		scores[b, "ploidy"] <- ploidy
		
		# Standard Deviation of Diffs
		if(nPeaks > 2L) scores[b, "sdd"] <- stats::sd(diff(xMax))
		
		if(nPeaks > 1L) {
			# Peaks to Model
			ptm <- copies(xMax, from="LCN", center=scores[b, "center"], width=scores[b, "width"], ploidy=scores[b, "ploidy"], exact=TRUE)
			ptm <- mean(abs(ptm - round(ptm)))
			scores[b, "ptm"] <- ptm
		}
		
		if(nPeaks > 0L) {
			# Segments to Model
			stm <- copies(segLCN[ ! segChroms %in% exclude ], from="LCN", center=scores[b, "center"], width=scores[b, "width"], ploidy=scores[b, "ploidy"], exact=TRUE)
			stm <- stats::weighted.mean(abs(stm - round(stm)), segLengths[ ! segChroms %in% exclude ])
			scores[b, "stm"] <- stm
		}
		
		# Repeat until bandwidth is too large
		bw <- bw + by
		if(bw > to) {
			break
		} else {
			b <- b + 1L
		}
	}
	
	# Selection
	filter <- which(!is.na(scores[,method]) & scores[,"peaks"] <= maxPeaks & scores[,"peaks"] >= minPeaks & scores[,"width"] >= minWidth & scores[,"width"] <= maxWidth)
	if(any(filter)) {
		# Best matching model
		output <- scores[ filter , , drop=FALSE ]
		output <- output[ which.min(output[,method]) , , drop=TRUE ]
	} else {
		# No matching model
		if(isTRUE(discreet)) {
			# Silently return a NA model
			output <- rbind(scores[0,], NA)[,,drop=TRUE]
			output["peakFrom"] <- peakFrom
			output["peakTo"] <- peakTo
			output["ploidy"] <- ploidy
		} else {
			stop("No satisfying 'bandwidth' provided")
		}
	}
	
	return(output)
}

