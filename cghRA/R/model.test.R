# Plot a copy-number model as produced by model.auto()
# Author : Sylvain Mareschal <mareschal@ovsa.fr>

model.test = function(
		segLogRatios,
		segChroms,
		segLengths = rep(1L, length(segLogRatios)),
		model = NA,
		center = model['center'],
		width = model['width'],
		ploidy = model['ploidy'],
		bw = model['bw'],
		minDensity = 0.001,
		peakFrom = model['peakFrom'],
		peakTo = model['peakTo'],
		graph = TRUE,
		parameters = TRUE,
		returnPar = FALSE,
		xlim = c(0, 5),
		ylim = c(0, max(segLengths)),
		xlab = "Segment copy number",
		ylab = "Segment length",
		cex.seg = 0.4,
		cex.leg = 0.7,
		cex.l2r = 0.7,
		exclude = c("X", "Y", "Xp", "Xq", "Yp", "Yq"),
		title = NULL,
		panel = FALSE,
		klim = NULL,
		l2r.points = c(35, 10),
		...
		)
	{
	# Length checks
	if(length(segLogRatios) == 0 || !is.numeric(segLogRatios))                stop("'segLogRatios' must be a non empty numeric vector")
	if(length(segLengths) != length(segLogRatios) || !is.integer(segLengths)) stop("'segLengths' must be a non empty integer vector")
	
	# Default ploidy
	if(is.na(ploidy)) ploidy <- 0
	
	# Log-ratio related Copy Numbers (LCN)
	segLCN <- LCN(segLogRatios, exact=TRUE)
	
	# Filter
	f <- ! segChroms %in% exclude & !is.na(segLCN)
	if(!any(f)) stop("Unable to plot density with only NAs")
	segLogRatios <- segLogRatios[f]
	segLengths <- segLengths[f]
	segChroms <- segChroms[f]
	segLCN <- segLCN[f]
	
	# Weights
	segWeights <- segLengths / sum(segLengths)
	
	if(!is.na(bw)) {
		# Segment density
		segDensity <- stats::density(x=segLCN, bw=bw, kernel="gaussian", weights=segWeights)
		
		# Local maxima
		y <- segDensity$y / max(segDensity$y)
		i <- 2:(length(y)-1)
		iMax <- 1 + which(y[i-1] <= y[i] & y[i+1] <= y[i] & y[i+1] != y[i-1] & y[i] > minDensity)
		xMax <- segDensity$x[ iMax ]
		
		# Peak filtering
		if(!is.na(peakFrom)) {
			fromLCN <- LCN(peakFrom, exact=TRUE)
			iMax <- iMax[xMax >= fromLCN]
			xMax <- xMax[xMax >= fromLCN]
		}
		if(!is.na(peakTo)) {
			toLCN <- LCN(peakTo, exact=TRUE)
			iMax <- iMax[xMax <= toLCN]
			xMax <- xMax[xMax <= toLCN]
		}
		
		# Peak amount
		peaks <- length(iMax)
	} else {
		peaks <- NA
	}
	
	# Model statistics
	if(!is.na(bw) && !is.na(center) && !is.na(width)) {
		# Standard Deviation of Diffs
		sdd <- stats::sd(diff(xMax))
		
		# Peaks to Model
		ptm <- copies(xMax, from="LCN", center=center, width=width, ploidy=ploidy, exact=TRUE)
		ptm <- mean(abs(ptm - round(ptm)))
		
		# Segments to Model
		stm <- copies(segLCN[ ! segChroms %in% exclude ], from="LCN", center=center, width=width, ploidy=ploidy, exact=TRUE)
		stm <- stats::weighted.mean(abs(stm - round(stm)), segLengths[ ! segChroms %in% exclude ], na.rm=TRUE)
	} else {
		sdd <- NA
		ptm <- NA
		stm <- NA
	}
	
	if(isTRUE(graph)) {
		# xlim in modelized copies
		if(is.numeric(klim) && length(klim) == 2 && all(!is.na(klim))) {
			# kilm is provided and valid
			if(!is.na(ploidy) && !is.na(width) && !is.na(center)) { xlim <- (klim - ploidy) * width + center
			} else                                                { stop("Unable to use 'klim' without a copy-number model")
			}
		}
		
		# Background for ploting (relative copies)
		if(isTRUE(panel)) { graphics::plot(x=NA, y=NA, xlim=rev(ylim), ylim=xlim, xlab="", xaxt="n", ylab="", yaxt="n", ...)
		} else            { graphics::plot(x=NA, y=NA, xlim=xlim,      ylim=ylim, xlab="", xaxt="n", ylab="Segment length", ...)
		}
		savePar <- graphics::par()
		
		# Peak range
		if(!is.na(peakFrom) & !isTRUE(panel)) {
			graphics::rect(
				xleft = graphics::par("usr")[1],
				xright = LCN(peakFrom, exact=TRUE),
				ybottom = graphics::par("usr")[3],
				ytop = graphics::par("usr")[4],
				col = "#DDDDDD",
				border = "#DDDDDD"
			)
		}
		if(!is.na(peakTo) & !isTRUE(panel)) {
			graphics::rect(
				xleft = LCN(peakTo, exact=TRUE),
				xright = graphics::par("usr")[2],
				ybottom = graphics::par("usr")[3],
				ytop = graphics::par("usr")[4],
				col = "#DDDDDD",
				border = "#DDDDDD"
			)
		}
		
		# LogRatios axis
		at <- pretty(xlim, n=l2r.points[panel+1L])
		at <- at[ at > 0L ]
		labels <- round(log(at/2, 2), 2)
		if(!isTRUE(panel)) {
			graphics::mtext(
				side = ifelse(isTRUE(panel), 2, 3),
				text = "LogRatios",
				line = 3
			)
		}
		graphics::axis(
			tck = 1,
			col = "#CCCCCC",
			side = ifelse(isTRUE(panel), 2, 3),
			at = at,
			labels = FALSE
		)
		graphics::mtext(
			side = ifelse(isTRUE(panel), 2, 3),
			at = at,
			cex = cex.l2r,
			text = ifelse(labels == 0 & isTRUE(panel), "logRatio", labels),
			las = ifelse(isTRUE(panel), 1, 3),
			line = ifelse(isTRUE(panel), 0.15, 0.5)
		)
		
		if(!is.na(center) && !is.na(width) && !isTRUE(panel)) {
			# Copies axis
			labels <- ploidy + (-3:3)
			at <- (labels - ploidy) * width + center
			if(ploidy == 0) labels <- c("-3", "-2", "-1", "n", "+1", "+2", "+3")
			graphics::axis(
				side = 1,
				at = at,
				labels = labels
			)
			graphics::title(xlab="Segment copy number")
		}
		
		# Real peaks
		if(!is.na(bw)) {
			if(isTRUE(panel)) { graphics::abline(h=xMax, col="#CC0000", lty="solid")
			} else            { graphics::abline(v=xMax, col="#CC0000", lty="solid")
			}
		}
		
		# Plot density
		if(!is.na(bw)) {
			graphics::par(new=TRUE)
			if(isTRUE(panel)) {
				graphics::plot(x=NA, y=NA, xlim=range(-segDensity$y), ylim=xlim, xlab="", ylab="", xaxt="n", yaxt="n", ...)
				graphics::polygon(x=-segDensity$y, y=segDensity$x, col="#AAAAAA", border="#AAAAAA")
			} else {
				graphics::plot(x=segDensity$x, y=segDensity$y, type="l", xlim=xlim, col="#AAAAAA", xlab="", ylab="", xaxt="n", yaxt="n", ...)
			}
		}
		
		# Plot segments
		graphics::par(new=TRUE)
		if(isTRUE(panel)) {
			pch <- 20L
			graphics::plot(x=segLengths, y=segLCN, xlim=rev(ylim), ylim=xlim, pch=pch, cex=cex.seg, xlab="", ylab="", xaxt="n", yaxt="n", ...)
		} else {
			pch <- rep(3L, length(segLCN))
			pch[ segChroms %in% exclude ] <- 1L
			graphics::plot(x=segLCN, y=segLengths, xlim=xlim,      ylim=ylim, pch=pch, cex=cex.seg, xlab="", ylab="", xaxt="n", yaxt="n", ...)
		}
		
		if(parameters) {
			# Parameters legend
			if(isTRUE(panel)) {
				graphics::text(x=graphics::par("usr")[1], y=graphics::par("usr")[3], labels=sprintf(" %.3f~%.3f\n STM = %.5f", center, width, stm), adj=c(0, -0.25), cex=cex.leg)
			} else {
				graphics::legend(
					x = "topleft",
					legend = c(
						sprintf("bw = %.3f   ", bw),
						sprintf("peaks = %i   ", peaks),
						sprintf("peakFrom = %.3f   ", peakFrom),
						sprintf("peakTo = %.3f   ", peakTo),
						sprintf("center = %.3f   ", center),
						sprintf("width = %.3f   ", width),
						sprintf("ploidy = %s   ", ifelse(ploidy == 0, "relative", as.character(ploidy))),
						sprintf("sdd = %.5f   ", sdd),
						sprintf("ptm = %.5f   ", ptm),
						sprintf("stm = %.5f   ", stm)
					),
					bg = "#EEEEEE",
					inset = 0.01,
					cex = cex.leg,
					title = title
				)
			}
		}
	} else {
		savePar <- NULL
	}
	
	# Return
	if(isTRUE(returnPar)) {
		invisible(savePar)
	} else {
		return(
			c(
				bw = as.double(bw),
				peaks = as.double(peaks),
				peakFrom = as.double(peakFrom),
				peakTo = as.double(peakTo),
				center = as.double(center),
				width = as.double(width),
				ploidy = as.double(ploidy),
				sdd = as.double(sdd),
				ptm = as.double(ptm),
				stm = as.double(stm)
			)
		)
	}
}

