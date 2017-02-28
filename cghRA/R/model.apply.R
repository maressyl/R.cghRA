# Transform log-ratios into copy-numbers using a copy-number model produced by model.auto()
# Author : Sylvain Mareschal <mareschal@ovsa.fr>

model.apply = function(
		segStarts,
		segEnds,
		segChroms,
		segLogRatios,
		segLengths,
		model = NA,
		center = model['center'],
		width = model['width'],
		ploidy = model['ploidy'],
		exact = FALSE,
		merge = TRUE
		)
	{
	# Checks
	if(length(segStarts) == 0 || !is.integer(segStarts))       stop("'segStarts' must be a non empty integer vector")
	if(length(segEnds) == 0 || !is.integer(segEnds))           stop("'segEnds' must be a non empty integer vector")
	if(length(segChroms) == 0 || !is.atomic(segChroms))        stop("'segChroms' must be a non empty atomic vector")
	if(length(segLogRatios) == 0 || !is.numeric(segLogRatios)) stop("'segLogRatios' must be a non empty numeric vector")
	if(length(segLengths) == 0 || !is.integer(segLengths))     stop("'segLengths' must be a non empty integer vector")
	if(is.na(center) || is.na(width) || is.na(ploidy))         stop("'center', 'width' and 'ploidy' must be provided directly or via 'model'")
	
	# Copies
	segCopies <- copies(x=segLogRatios, center=center, width=width, ploidy=ploidy, exact=exact, from="logRatios")
	if(!isTRUE(exact) && isTRUE(merge)) {
		# Ordering
		segOrder <- order(segChroms, segStarts)
		segStarts <- segStarts[ segOrder ]
		segEnds <- segEnds [ segOrder ]
		segChroms <- segChroms[ segOrder ]
		segLogRatios <- segLogRatios[ segOrder ]
		segLengths <- segLengths[ segOrder ]
		
		# Segments to merge
		g <- 1
		segGroups <- integer(length(segChroms))
		segGroups[1] <- g
		for(i in 2:length(segChroms)) {
			if(segChroms[i-1] != segChroms[i] || segEnds[i-1] != segStarts[i] || segCopies[i-1] != segCopies[i]) g <- g + 1
			segGroups[i] <- g
		}
		
		# Merging
		segStarts <- tapply(X=segStarts, INDEX=segGroups, FUN=min)
		segEnds <- tapply(X=segEnds, INDEX=segGroups, FUN=max)
		segChroms <- tapply(X=segChroms, INDEX=segGroups, FUN=unique, simplify=FALSE)    # Returns an array of mode 'list', to preserve factors
		segCopies <- tapply(X=segCopies, INDEX=segGroups, FUN=unique)
		segLengths <- tapply(X=segLengths, INDEX=segGroups, FUN=sum)
		
		# Arrays to vectors
		attributes(segStarts) <- list()
		attributes(segEnds) <- list()
		attributes(segChroms) <- list()
		segChroms <- unlist(segChroms)                                                  # Unlist to a vector
		attributes(segCopies) <- list()
		attributes(segLengths) <- list()
		
		# Theoretic logRatios
		if(ploidy == 0) segLogRatios <- rep(as.numeric(NA), length(segCopies))
		else            segLogRatios <- copies(x=segCopies, ploidy=ploidy, from="copies")
	}
	
	return(
		data.frame(
			segStarts = segStarts,
			segEnds = segEnds,
			segChroms = segChroms,
			segLogRatios = segLogRatios,
			segCopies = segCopies,
			segLengths = segLengths,
			stringsAsFactors = FALSE
		)
	)
}
