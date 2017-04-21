# Reference class for the modelized regions element of cghRA objects
# Author : Sylvain Mareschal <mareschal@ovsa.fr>

# R5 sub-class definition
setRefClass(
	Class = "cghRA.copies",
	contains = "cghRA.regions",
	methods = list(

check = function(warn=TRUE) {
"Raises an error if the object is not valid, else returns TRUE"
	
	# refTable, track.table, cghRA.regions
	callSuper(warn)
	
	# Copies column
	if(! "copies" %in% .self$getColNames()) stop("Must contain a 'copies' column")
	if(.self$types("copies") != "numeric")  stop("'copies' column must be of class 'numeric'")
	
	return(TRUE)
},

defaultParams = function(...) {
"Returns class-specific defaults for graphical parameters. Inheriting class should overload it to define their own defaults.
- ...   : may be used by inheriting methods, especially for inter-dependant parameters."
	
	params <- callSuper(...)
	
	if(.self$modelized()) {
		# cghRA model
		params$column <- "copies"
		params$ylim <- .self$model['ploidy'] + c(-3, 3)
		params$origin <- .self$model['ploidy']
	} else {
		# May be relevant if issued by an external copy model
		params$column <- "copies"
		params$ylim <- c(-3, 3)
		params$origin <- 0
	}
	
	return(params)
},

permute = function(design, check=TRUE, replace=FALSE, storeHash=FALSE, keepTer=TRUE, quiet=TRUE, seed=NULL) {
"Relocate altered segments randomly in the genome preserving their sizes (in probe counts) and copy number.
- design      : the cghRA.design object corresponding to this array.
- check       : single logical value, whether to check permutated data consistency or not.
- replace     : single logical value, whether to update the object itself or return a data.frame with permuted coordinates.
- storeHash   : single logical value, whether to store the coordinate hash in 'design' to enhance speed or not.
- keepTer     : single logical value, whether to force starting/ending segments of chromosome to terminal chromosome positions or not.
- quiet       : single logical value, whether to display how many possibilities were left to relocate a segment or not.
- seed        : single integer value to control randomness (see set.seed manual for further details).
"
	
	# Reproducible randomness
	if(!is.null(seed)) set.seed(seed)
	
	# Hash for design coordinates
	if("startHash" %in% design$getColNames()) {
		startHash <- design$extract(,"startHash")
	} else {
		startHash <- sprintf("%s:%i", design$extract(,"chrom"), design$extract(,"start"))
		if(isTRUE(storeHash)) design$addColumn(startHash, "startHash")
	}
	
	# Permuted segment storage
	out <- NULL
	
	# Last probe of each free ranges (starting with whole chromosomes)
	freeRangeEnds <- design$index
	
	# Probes occupied by already permuted segments
	occupied <- integer(0)
	
	# Convert genomic coordinates to probe indexes
	startCoords <- sprintf("%s:%i", .self$extract(,"chrom"), .self$extract(,"start"))
	endCoords <- sprintf("%s:%i", .self$extract(,"chrom"), .self$extract(,"end"))
	startIndexes <- match(startCoords, startHash)
	endIndexes <- match(endCoords, startHash)
	if(any(is.na(startIndexes))) stop("Some segment 'start' don't correspond to any probe 'start' in 'design'")
	if(any(is.na(endIndexes)))   stop("Some segment 'end' don't correspond to any probe 'start' in 'design'")
	
	# Is segment starting or ending a chromosome ?
	isTerminal <- !duplicated(.self$extract(,"chrom")) | !duplicated(.self$extract(,"chrom"), fromLast=TRUE)
	
	# Compute size (in probes)
	sizes <- endIndexes - startIndexes + 1L
	
	# Process altered segments from the widest to the smallest
	segIndexes <- which(.self$extract(,"copies") != 0L)
	segIndexes <- segIndexes[ order(sizes[ segIndexes ], decreasing=TRUE) ]
	
	# Permuted segment storage
	out <- .self$extract(segIndexes)
	out$.oldSize <- sizes[segIndexes]
	out$.startIndex <- as.integer(NA)
	out$.endIndex <- as.integer(NA)
	
	# Processed altered segment one by one
	for(i in 1:length(segIndexes)) {
		segIndex <- segIndexes[i]
		
		if(isTRUE(keepTer) && isTerminal[segIndex]) {
			# Start by allowing no probe as new potential starting points
			ptl <- rep(FALSE, tail(design$index, 1L))
			
			# Allow to stick to chromosome start (chromosomes too short are filtered out later)
			ptl[ c(1L, head(design$index, -1)+1L) ] <- TRUE
			
			# Allow to stick to chromosome end (must stay on current chromosome)
			cdd <- design$index - sizes[segIndex] + 1L
			cdd <- cdd[ cdd >= c(1L, head(design$index, -1)) ]
			ptl[ cdd ] <- TRUE
		} else {
			# Start by considering all mapped probes as new potential starting points
			ptl <- rep(TRUE, tail(design$index, 1L))
		}
		
		# Filter out starting points too close to a chromosome end or a permuted segment start
		for(k in 1:length(freeRangeEnds)) {
			mask <- freeRangeEnds[k] - (0L : (sizes[segIndex] - 2L))
			mask <- mask[ mask > 0L ]
			ptl[ mask ] <- FALSE
		}
		
		# Filter out probes occupied by previous permuted segments
		ptl[ occupied ] <- FALSE
		
		# Choose a new starting point
		n <- sum(ptl)
		if(n == 0L) stop("No more room to place a segment")
		if(!isTRUE(quiet)) {
			message(
				"Choosing 1 new position out of ",
				sum(ptl),
				" for a ",
				sizes[segIndex],
				" probe long ",
				ifelse(isTRUE(keepTer) && isTerminal[segIndex], "terminal ", ""),
				"segment"
			)
		}
		newStartIndex <- sample(which(ptl), 1L)
		newEndIndex <- newStartIndex + sizes[segIndex] - 1L
		
		# Store new location index
		out$.startIndex[i] <- newStartIndex
		out$.endIndex[i] <- newEndIndex
		
		# Remind permuted segment to prevent overlap
		freeRangeEnds <- c(freeRangeEnds, newStartIndex - 1L)
		occupied <- c(occupied, newStartIndex:newEndIndex)
	}
	
	# Sort permuted segments
	out <- out[ order(out$.startIndex) ,]
	
	# Add gap-filling WT segments between all permuted segments
	segEnds <- c(design$index, out$.startIndex - 1L)
	segStarts <- c(1L, head(design$index + 1L, -1L), out$.endIndex + 1L)
	wt <- out[NA,][1:length(segStarts),]
	wt$.startIndex <- sort(segStarts)
	wt$.endIndex <- sort(segEnds)
	wt$.oldSize <- as.integer(NA)
	wt$copies <- 0L
	
	# Exclude empty filling segments (occurs when 2 permuted segments end next to each others)
	wt <- wt[ wt$.startIndex <= wt$.endIndex ,]
	
	# Merge altered and WT segments
	out <- rbind(out, wt)
	out <- out[ order(out$.startIndex) ,]
	
	# Replace genomic coordinates
	out$chrom <- design$extract(out$.startIndex, "chrom")
	out$start <- design$extract(out$.startIndex, "start")
	out$end <- design$extract(out$.endIndex, "start")
	
	# Recompute probe counts (higher than original ones as no probe is now excluded)
	out$probes <- out$.endIndex - out$.startIndex + 1L
	
	# Cosmetics
	rownames(out) <- NULL
	
	# Checks
	if(isTRUE(check)) {
		if(any(tail(out$.startIndex, -1) != head(out$.endIndex, -1) + 1L))                          stop("Index coordinate inconsistency")
		if(any(any(out$.endIndex - out$.startIndex + 1L != out$probes)))                            stop("Permuted segment length inconsistency")
		if(any(out$.startIndex > as.integer(tail(design$index, 1))))                                stop("Out of range start index")
		if(any(out$.endIndex > as.integer(tail(design$index, 1))))                                  stop("Out of range end index")
		if(any(design$extract(out$.startIndex, "chrom") != design$extract(out$.endIndex, "chrom"))) stop("Start and end probes on dictinct chromosomes")
		if(any(! design$index %in% out$.endIndex))                                                  stop("Some chromosome is not covered until the end")
		if(any(!c(1L, head(design$index, -1L) + 1L) %in% out$.startIndex))                          stop("Some chromosome is not covered from the start")
		if(any(!is.na(out$.oldSize) && out$probes != out$.oldSize))                                 stop("Some segments lose their probe sizes")
	}
	
	# Remove indexes
	out$.startIndex <- NULL
	out$.endIndex <- NULL
	out$.oldSize <- NULL
	
	if(isTRUE(replace)) {
		# Replace content
		.self$erase()
		.self$addDataFrame(out)
	} else {
		# Return new content
		return(out)
	}
},

show = function(include=FALSE, fieldWidth=10) {
"Interactive printing
- include   : single logical value, if TRUE class name will not be printed."
	
	# Class name
	if(!isTRUE(include)) {
		cat("\n  \"cghRA.copies\" reference class object\n")
	} else {
		cat("\n  Extends \"cghRA.copies\"\n")
	}
	
	# Inherited show()
	callSuper(include=TRUE, fieldWidth=fieldWidth)
}

	)
)

# Constructor
cghRA.copies <- function(..., warn=TRUE) {
	# Inheritance
	object <- new("cghRA.copies")
	object$import(cghRA.regions(..., warn=warn))
	
	# Check
	object$check(warn=warn)
	
	return(object)
}
