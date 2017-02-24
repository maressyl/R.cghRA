# Reference class for genomic segments
# Author : Sylvain Mareschal <mareschal@ovsa.fr>

# R5 sub-class definition
setRefClass(
	Class = "cghRA.regions",
	fields = list(
		model = "numeric",
		segmentCall = "call",
		modelizeCall = "call"
	),
	contains = "track.table",
	methods = list(

check = function(warn=TRUE) {
"Raises an error if the object is not valid, else returns TRUE"
	
	# refTable, track.table
	callSuper(warn)
	
	# Probes column
	if(! "probes" %in% .self$getColNames())   stop("Must contain a 'probes' column")
	if(.self$types("probes") != "integer")    stop("'probes' column must be of class 'integer'")
	
	# LogRatio column
	if(! "logRatio" %in% .self$getColNames()) stop("Must contain a 'logRatio' column")
	if(.self$types("logRatio") != "numeric")  stop("'logRatio' column must be of class 'numeric'")
	
	return(TRUE)
},

defaultParams = function(...) {
"Returns class-specific defaults for graphical parameters. Inheriting class should overload it to define their own defaults.
- ...   : may be used by inheriting methods, especially for inter-dependant parameters."
	
	params <- callSuper(...)
	
	params$drawFun <- "draw.hist"
	params$column <- "logRatio"
	params$origin <- 0
	params$ylim <- c(-3, 3)
	params$yaxt <- "s"
	params$colorVal <- NA
	params$colorFun <- function() {
		output = rep(as.character(NA), nrow(slice))
		output[ slice[[column]] == origin ] = "#888888"
		output[ slice[[column]] > origin ] = "#8888FF"
		output[ slice[[column]] < origin ] = "#FF8888"
		return(output)
	}
	
	return(params)
},

fillGaps = function(...) {
"Apply the fillGaps() function to extend regions in order to fill inter-segment gaps."
	
	# Processing
	newRegions <- cghRA::fillGaps(
		segTable = .self$extract(),
		...
	)
	
	# Region field replacing
	.self$erase()
	.self$addDataFrame(newRegions)
	
	# Check
	.self$check()
},

initialize = function(model=NA_real_, segmentCall=new("call"), modelizeCall=new("call"), ...) {
	callSuper(...)
	initFields(model=model, segmentCall=segmentCall, modelizeCall=modelizeCall)
},


karyotype = function(bandTrack, value="logRatio", thresholds=c(-0.3, 0.3), precision=2) {
"Returns a karyotype formula of altered regions.
- bandTrack    : a track.table object, as returned by track.UCSC_bands().
- value        : column to use to select altered regions.
- thresholds   : length 2 numeric vector defining altered values.
- precision    : single integer value from 1 to 4, amount of digits to consider in banding."
	
	# Merge arms
	regions <- .self$eraseArms(temp=TRUE)
	
	# Altered regions
	alter <- rep(as.character(NA), regions$rowCount)
	alter[ regions$extract(,value) < thresholds[1] ] <- "del" 
	alter[ regions$extract(,value) > thresholds[2] ] <- "dup"
	
	# Cytoband boundaries
	mergedBands <- bandTrack$eraseArms(temp=TRUE)
	loci <- strsplit(regions$cross(mergedBands, type="name", maxElements=NA), split=", ")
	
	# Subset on altered regions
	loci <- loci[ !is.na(alter) ]
	alter <- alter[ !is.na(alter) ]
	
	# Locus boundaries
	loci <- cbind(sapply(loci, head, 1), sapply(loci, tail, 1))
	
	# Band precision
	loci <- switch(
		as.character(as.integer(precision)),
		"1" = sub("[0-9](\\.[0-9]+)?$", "0", loci),
		"2" = sub("\\.[0-9]+$", "", loci),
		"3" = sub("(?<=\\.[0-9])[0-9]$", "", loci, perl=TRUE),
		"4" = loci,
		stop("Invalid precision (from 1 to 4)")
	)
	
	# Reshape loci
	loci <- apply(loci, 1, unique)
	for(i in 1:length(loci)) loci[[i]] <- sub(sprintf("^%s", names(loci)[i]), "", loci[[i]])
	loci <- sapply(loci, paste, collapse=";")
	
	# Formula build
	karyo <- paste(sprintf("%s(%s)(%s)", alter, names(loci), loci), collapse=", ")
	return(karyo)
},

show = function(include=FALSE, fieldWidth=10) {
"Interactive printing
- include   : single logical value, if TRUE class name will not be printed."
	
	# Class name
	if(!isTRUE(include)) {
		cat("\n  \"cghRA.regions\" reference class object\n")
	} else {
		cat("\n  Extends \"cghRA.regions\"\n")
	}
	
	# Fields
	cat(sprintf("  %-*s : <%s>\n", fieldWidth, "model", .self$modelized()))
	
	# Inherited show()
	callSuper(include=TRUE, fieldWidth=fieldWidth)
},

status = function(chrom, start, end, value="copies", na=c("fill", "keep"), fuzzy=FALSE, states=list(deletion=c(-Inf, -0.5), gain=c(0.5, Inf))) {
"Returns the copy states in various windows, mimicing penetrance behavior.
- chrom    : character vector, chromosome location of the regions to query.
- start    : integer vector, starting position on the chromosome for the regions to query.
- end      : integer vector, ending position on the chromosome for the regions to query.
- value    : single character value, name of the column to use for state assignation.
- na       : single character value, see penetrance() help page for details ('false' is not handled).
- fuzzy    : single logical value, whether to assign the state when some sub-regions are out or not.
- states   : list of states, see penetrance help page for details."
	
	na <- match.arg(na)
	chrom <- as.character(chrom)
	output <- matrix(
		data = as.logical(NA),
		nrow = length(chrom),
		ncol = length(states),
		dimnames = list(NULL, names(states))
	)
	
	for(i in 1:length(chrom)) {
		# Copies in the region
		x <- slice(chrom[i], start[i], end[i])[[ value ]]
		
		# NA : try to fill
		if(length(x) == 0L) {
			if(na == "fill") {
				chromSeg <- slice(chrom[i], 0, max(extract(,"end")))
				boundaries <- as.vector(rbind(chromSeg$start, chromSeg$end))
				if(start[i] > max(boundaries)) {
					# After chromosome end
					x <- chromSeg[ nrow(chromSeg), value ]
				} else if(end[i] < min(boundaries)) {
					# Before chromosome start
					x <- chromSeg[ 1L, value ]
				} else {
					# Between segments
					extendedStart <- boundaries[ tail(which(boundaries < start[i]), 1) ]
					extendedEnd <- boundaries[ head(which(boundaries > end[i]), 1) ]
					x <- slice(chrom[i], extendedStart, extendedEnd)[[ value ]]
				}
			} else {
				# Keep NA
				output[i,] <- NA; next;
			}
		}
		
		for(stateName in names(states)) {
			# State attribution
			if(length(states[[stateName]]) == 2)      s <- x >= states[[stateName]][1] & x < states[[stateName]][2]
			else if(length(states[[stateName]]) == 1) s <- x == states[[stateName]]
			else                                      stop("Each element in 'states' must contain two boundaries or a single value")
			
			# Consensus
			if(isTRUE(fuzzy))               output[i,stateName] <- any(s)
			else if(length(unique(s)) == 1) output[i,stateName] <- unique(s)
			else                            output[i,stateName] <- NA
		}
	}
	
	return(output)
},

modelized = function() {
"Does the object embed a complete model or not"
	
	needed <- c("bw", "peaks", "peakFrom", "peakTo", "center", "width", "ploidy")
	return(all(needed %in% names(model)) && all(!is.na(model[ needed ])))
},

model.apply = function(...) {
"Call the model.apply() function to produce a cghRA.copies object with predicted copy number for each region."
	
	# Copy calling
	seg <- cghRA::model.apply(
		segStarts = .self$extract(,"start"),
		segEnds = .self$extract(,"end"),
		segChroms = .self$extract(,"chrom"),
		segLogRatios = .self$extract(,"logRatio"),
		segLengths = .self$extract(,"probes"),
		model = model,
		...
	)
	
	# cghRA.copies object
	output <- cghRA.copies(
		chrom = seg$segChroms,
		strand = factor(rep(NA, nrow(seg)), levels=c("-","+")),
		start = seg$segStarts,
		end = seg$segEnds,
		probes = seg$segLengths,
		logRatio = seg$segLogRatios,
		copies = seg$segCopies,
		.name = .self$name,
		.organism = .self$organism,
		.assembly = .self$assembly,
		.model = .self$model,
		.parameters = .self$parameters,
		.makeNames = TRUE
	)
	
	return(output)
},

model.auto = function(save=TRUE, ...) {
"Call the model.auto() function to automatically fit a copy-number prediction model.
- save   : single logical value, whether to save the model or only return it"
	
	# Apply
	tmpModel <- cghRA::model.auto(
		segLogRatios = .self$extract(,"logRatio"),
		segChroms = .self$extract(,"chrom"),
		segLengths = .self$extract(,"probes"),
		...
	)
	
	# Return
	if(isTRUE(save)) {
		model <<- tmpModel
		invisible(tmpModel)
	} else {
		return(tmpModel)
	}			
},

model.test = function(...) {
"Call the model.test() function to plot the current copy-number model."
	
	# Apply and return
	return(
		cghRA::model.test(
			segLogRatios = .self$extract(,"logRatio"),
			segChroms = .self$extract(,"chrom"),
			segLengths = .self$extract(,"probes"),
			model = .self$model,
			title = .self$name,
			...
		)
	)
},
		
proportions = function(chrom=chromosomes(), value="copies", states=list(deletion=c(-Inf, -0.5), gain=c(0.5, Inf)), mode=c("bp", "probes")) {
"Returns the proportion of the chromosomes in given states (in bp involved).
- chrom    : character vector, chromosome location of the regions to query. Consider track.table$eraseArms() to focus on chromosome arms.
- value    : single character value, name of the column to use for state assignation.
- states   : list of states, see penetrance help page for details."
	
	# Arg check
	mode <- match.arg(mode)
	
	# Output matrix
	output <- matrix(
		data = as.double(NA),
		nrow = length(chrom),
		ncol = length(states),
		dimnames = list(chrom, names(states))
	)
	
	for(k in chrom) {
		# Get regions
		expr <- parse(text=sprintf("chrom == \"%s\"", k))
		regions <- extract(expr)
		
		if(nrow(regions) > 0) {				
			for(stateName in names(states)) {
				# State attribution
				if(length(states[[stateName]]) == 2)      isInState <- regions[[value]] >= states[[stateName]][1] & regions[[value]] < states[[stateName]][2]
				else if(length(states[[stateName]]) == 1) isInState <- regions[[value]] == states[[stateName]]
				else                                      stop("Each element in 'states' must contain two boundaries or a single value")
				
				# Computation
				if(mode == "bp")            { output[ k , stateName ] <- sum(regions[isInState,"end"] - regions[isInState,"start"]) / (max(regions$end) - min(regions$start))
				} else if(mode == "probes") { output[ k , stateName ] <- sum(regions[isInState,"probes"]) / sum(regions$probes)
				}
			}
		}
	}
	
	return(output)
}

	)
)

# Constructor
cghRA.regions <- function(..., .model, warn=TRUE) {
	# Inheritance
	object <- new("cghRA.regions")
	if(!missing(.model)) object$model <- .model
	object$import(track.table(..., warn=warn))
	
	# Check
	object$check(warn=warn)
	
	return(object)
}

