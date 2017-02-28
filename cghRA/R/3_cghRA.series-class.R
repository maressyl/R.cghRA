# Reference class for series of cghRA.regions objects
# Author : Sylvain Mareschal <mareschal@ovsa.fr>

# R5 class definition
setRefClass(
	Class = "cghRA.series",
	fields = list(
		name = "character",
		arrays = "list"
	),
	methods = list(

add = function(object) {
"Add an object to the series"
	
	if(!is(object, "cghRA.regions"))           stop("'object' must be a 'cghRA.regions' object")
	if(object$name %in% .self$getArrayNames()) stop("'object' name (\"", object$name, "\") is already used in this series")
	
	arrays[[ object$name ]] <<- object
},

applyMethod = function(.method, ..., .simplify=TRUE, .quiet=TRUE) {
"Calls a method on each array of the series
- .method     : single character value, the method to be called.
- ...         : arguments to be passed to the method.
- .simplify   : same behavior as sapply() 'simplify' argument.
- .quiet      : single logical value, whether to print iterations or not."
	
	# Modified from Martin Morgan's answer on Stack Overflow
	# http://stackoverflow.com/questions/5841339/using-notation-for-reference-class-methods
	
	# Collect results
	returns <- list()
	for(i in 1:length(arrays)) {
		if(!isTRUE(.quiet)) message("Applying ", .method, "() to ", names(arrays)[i])
		returns[[i]] <- eval(substitute(arrays[[i]]$FUN(...), list(FUN=.method, ...)))
	}
	
	# Reshape results (inspired from sapply())
	if(length(returns) > 0) names(returns) <- names(arrays)
	if(!identical(.simplify, FALSE) && length(returns)) {
		returns <- simplify2array(returns, higher=(.simplify == "array"))
	}
	
	invisible(returns)
},

check = function(warn=TRUE) {
"Raises an error if the object is not valid, else returns TRUE"
	
	# Fields
	if(length(name) != 1)             stop("'name' must be a single character value")
	
	# Arrays
	for(arr in arrays) {
		if(!is(arr, "cghRA.regions")) stop("'arrays' must be a list of 'cghRA.regions' objects")
		arr$check(warn=warn)
	}
	
	# Warnings
	if(isTRUE(warn)) {
		if(is.na(name))               warning("'name' should not be NA")
	}
	
	return(TRUE)
},

get = function(arrayName) {
"Returns an element from the series"
	
	if(!arrayName %in% .self$getArrayNames()) stop("'arrayName' can not be found in the series")
	return(.self$arrays[[arrayName]])
},

getArrayNames = function() {
"Returns a vector of array names"
	
	return(names(.self$arrays))
},

initialize = function(name=NA_character_, arrays=list(), ...) {
	initFields(name=name, arrays=arrays)
},

last = function() {
"Refers to the last array added in the series"
	
	if(length(arrays) == 0) stop("The 'arrays' list is empty")
	return(tail(arrays, 1)[[1]])
},

LRA = function(value = "logRatio", tracks = TRUE, ...) {
"Apply the LRA() function to list Long Recurrent Abnormalities (Lenz et al, PNAS 2008).
- value    : single character value, the name of the column to use as copy number estimate ('copies' or 'logRatio').
- tracks   : single logical value, whether to convert output to track.table class or not."
	
	# All segments
	allSeg <- .self$pool(tracks=FALSE, value=value, states=NULL)
	
	# Computation
	output <- cghRA::LRA(
		segTables = allSeg,
		value = value,
		...
	)
	
	# Conversion
	if(isTRUE(tracks)) {
		if(length(arrays) == 0) stop("Unable to guess annotation from series content (empty series)")
		
		lraTracks <- list()
		for(state in names(output)) {
			# Track conversion
			output[[ state ]]$strand <- rep("+", nrow(output[[ state ]]))
			lraTracks[[ state ]] <- track.table(
				output[[ state ]],
				.name = sprintf("%s LRA (%s)", name, state),
				.organism = arrays[[1]]$organism,
				.assembly = arrays[[1]]$assembly,
				.chromosomes = arrays[[1]]$chromosomes(),
				.makeNames = TRUE,
				warn = FALSE
			)
			
			# Graphical parameters
			lraTracks[[ state ]]$setParam("drawFun", "draw.steps")
			lraTracks[[ state ]]$setParam("startColumns", c("extended.start", "start", "overlap.start"))
			lraTracks[[ state ]]$setParam("endColumns", c("extended.end", "end", "overlap.end"))
			lraTracks[[ state ]]$setParam("colorVal", NA)
			lraTracks[[ state ]]$setParam("border", "#000000")
			lraTracks[[ state ]]$setParam(
				"colorFun",
				function() {
					grey(1 - slice$value)
				}
			)
		}
		
		# Output
		return(lraTracks)
	} else {
		# Raw output
		return(output)
	}
},

parallelize = function(value = "logRatio", quiet = FALSE, tracks = TRUE, ...) {
"Apply the parallelize() function to build a summary matrix of the series.
- tracks   : single logical value, whether to convert output to track.table class or not."
	
	# All segments
	allSeg <- pool(tracks=FALSE, value=value, states=NULL, quiet=quiet)
	
	# Parallelize
	para <- cghRA::parallelize(
		segTables = allSeg,
		value = value,
		...
	)
	
	# Conversion
	if(isTRUE(tracks)) {
		if(!isTRUE(quiet)) message("Converting to track.table ... ", appendLF=FALSE)
		if(length(arrays) == 0) stop("Unable to guess annotation from series content (empty series)")
		
		# Track conversion
		para$strand <- rep("+", nrow(para))
		paraTrack <- track.table(
			para,
			.name = sprintf("%s parallelization (%s)", name, value),
			.organism = arrays[[1]]$organism,
			.assembly = arrays[[1]]$assembly,
			.chromosomes = arrays[[1]]$chromosomes(),
			.makeNames = TRUE,
			warn = FALSE
		)
		
		if(!isTRUE(quiet)) message("done")
		
		# Output
		return(paraTrack)
	} else {
		# Raw output
		return(para)
	}
},

penetrance = function(tracks = TRUE, ...) {
"Apply the penetrance() function to compute the proportion of altered samples for each genomic position.
- tracks   : single logical value, whether to convert output to track.table class or not."
	
	# Argument dispatch
	allArgs = list(...)
	if(any(names(allArgs) == "")) stop("All arguments must be named for dispatch")
	peneArgs = allArgs[ intersect(names(allArgs), names(formals(cghRA::penetrance))) ]
	paraArgs = allArgs[ intersect(names(allArgs), names(formals(cghRA::parallelize))) ]
	
	# Parallelization
	paraArgs$tracks <- FALSE
	peneArgs$segParallel <- do.call(
		what = .self$parallelize,
		args = paraArgs
	)
	
	# Penetrance
	pene <- do.call(
		what = cghRA::penetrance,
		args = peneArgs
	)
	
	# Conversion
	if(isTRUE(tracks)) {
		if(length(arrays) == 0) stop("Unable to guess annotation from series content (empty series)")
		
		peneTracks <- list()
		for(state in names(pene)) {
			# Track conversion
			pene[[ state ]]$strand <- rep("+", nrow(pene[[ state ]]))
			peneTracks[[ state ]] <- track.table(
				pene[[ state ]],
				.name = sprintf("%s penetrance (%s)", name, state),
				.organism = arrays[[1]]$organism,
				.assembly = arrays[[1]]$assembly,
				.chromosomes = arrays[[1]]$chromosomes(),
				.makeNames = TRUE,
				warn = FALSE
			)
			
			# Drawing parameters
			peneTracks[[ state ]]$setParam("drawFun", "draw.hist")
			if(state == "gain")     {
				peneTracks[[ state ]]$setParam("colorVal", "#8888FF")
				peneTracks[[ state ]]$setParam("border", "#8888FF")
			}
			if(state == "deletion") {
				peneTracks[[ state ]]$setParam("colorVal", "#FF8888")
				peneTracks[[ state ]]$setParam("border", "#FF8888")
			}
			if(state == "loss") {
				peneTracks[[ state ]]$setParam("colorVal", "#880000")
				peneTracks[[ state ]]$setParam("border", "#880000")
			}
		}
		
		# Output
		return(peneTracks)
	} else {
		# Raw output
		return(pene)
	}
},

pool = function(tracks = TRUE, value = "copies", group = TRUE, states = list(deletion=c(-Inf, -0.5), gain=c(0.5, Inf)), others=NULL, quiet=FALSE) {
"Collect and pool all alterated segments from the various samples of the series.
- tracks   : single logical value, whether to convert output to track.table class or not.
- value    : column on which apply a filtering.
- group    : single logical value, whether to visually group segments per samples or not (valid only for tracks=TRUE).
- states   : list of states, see penetrance help page for details. If 'states' is not empty, segments without state will be filtered out.
- others   : character vector, names of other columns to keep.
- quiet    : single logical value, whether to throw diagnosis messages or not."
	
	# All segments
	allSeg <- NULL
	if(!isTRUE(quiet)) message("Collecting segments ... ", appendLF=FALSE)
	for(arrayName in .self$getArrayNames()) {
		allSeg <- rbind(
			allSeg,
			cbind(
				sample = arrayName,
				.self$get(arrayName)$extract(, c("chrom", "start", "end", value, others))
			)
		)
	}
	if(!isTRUE(quiet)) message("done")
	
	# State calling
	if(!isTRUE(quiet) && length(states) > 0) message("Calling state ", appendLF=FALSE)
	for(stateName in names(states)) {
		if(!isTRUE(quiet)) message(stateName, ", ", appendLF=FALSE)
		if(length(states[[stateName]]) == 2) {
			# Interval
			isInState <- allSeg[[value]] >= states[[stateName]][1] & allSeg[[value]] < states[[stateName]][2]
		} else if(length(states[[stateName]]) == 1) {
			# Single value
			isInState <- allSeg[[value]] == states[[stateName]]
		} else {
			stop("Each element in 'states' must contain two boundaries or a single value")
		}
		
		# Attribution
		allSeg[ isInState , "type" ] <- stateName
	}
	if(!isTRUE(quiet) && length(states) > 0) message("done")
	
	# Filtering
	if(length(states) > 0) allSeg <- allSeg[ !is.na(allSeg$type) ,]
	
	# Conversion
	if(isTRUE(tracks)) {
		if(!isTRUE(quiet)) message("Converting to track.CNV ... ", appendLF=FALSE)
		if(length(arrays) == 0) stop("Unable to guess annotation from series content (empty series)")
		
		# Names built on sample names
		newNames <- make.unique(as.character(allSeg$sample), sep="#")
		newNames[ !grepl("#[0-9]+$", newNames) ] <- sprintf("%s#0", newNames[ !grepl("#[0-9]+$", newNames) ])
		allSeg$name <- newNames
		
		# Track conversion
		allSeg$strand <- "+"
		poolTrack <- track.CNV(
			allSeg,
			.name = sprintf("%s pool (%s)", name, value),
			.organism = arrays[[1]]$organism,
			.assembly = arrays[[1]]$assembly,
			.chromosomes = arrays[[1]]$chromosomes(),
			.makeNames = FALSE,
			warn = FALSE
		)
		
		# Drawing parameters tunning
		poolTrack$setParam("height", 1)
		
		if(isTRUE(group)) {
			# Merging
			poolTrack$buildGroupPosition(groupBy="sample")
			poolTrack$buildGroupSize(groupBy="sample")
			poolTrack$setParam("groupBy", "sample")
			poolTrack$setParam("border", "#000000")
		}
						
		if(!isTRUE(quiet)) message("done")
		
		# Output
		return(poolTrack)
	} else {
		# Raw output
		return(allSeg)
	}
},

show = function(include=FALSE, fieldWidth=10) {
"Interactive printing
- include   : single logical value, if TRUE class name will not be printed."
	
	# Class name
	if(!isTRUE(include)) {
		cat("\n  \"cghRA.series\" reference class object\n")
	} else {
		cat("\n  Extends \"cghRA.series\"\n")
	}
	
	# Fields
	cat(sprintf("  %-*s : %s\n", fieldWidth, "name", name))
	
	cat("\n")
	# Regions summary
	indexList <- names(arrays)
	if(is.null(indexList)) indexList <- rep("", length(arrays))
	nameList <- character(0)
	countList <- character(0)
	if(length(arrays) > 0) for(i in 1:length(arrays)) {
		nameList[i] <- arrays[[i]]$name
		countList[i] <- arrays[[i]]$rowCount
	}
	print(data.frame(index=indexList, name=nameList, segments=countList))
},

SRA = function(value = "logRatio", tracks = TRUE, ...) {
"Apply the SRA() function to list Short Recurrent Abnormalities (Lenz et al, PNAS 2008).
- value    : single character value, the name of the column to use as copy number estimate ('copies' or 'logRatio').
- tracks   : single logical value, whether to convert output to track.table class or not."
	
	# All segments
	allSeg <- .self$pool(tracks=FALSE, value=value, states=NULL)
	
	# Computation
	output <- cghRA::SRA(
		segTables = allSeg,
		value = value,
		...
	)
	
	# Conversion
	if(isTRUE(tracks)) {
		if(length(arrays) == 0) stop("Unable to guess annotation from series content (empty series)")
		
		sraTracks <- list()
		for(state in names(output)) {
			# Track conversion
			output[[ state ]]$strand <- rep("+", nrow(output[[ state ]]))
			sraTracks[[ state ]] <- track.table(
				output[[ state ]],
				.name = sprintf("%s SRA (%s)", name, state),
				.organism = arrays[[1]]$organism,
				.assembly = arrays[[1]]$assembly,
				.chromosomes = arrays[[1]]$chromosomes(),
				.makeNames = TRUE,
				warn = FALSE
			)
			
			# Graphical parameters
			sraTracks[[ state ]]$setParam("drawFun", "draw.steps")
			sraTracks[[ state ]]$setParam("startColumns", c("extended.start", "start", "overlap.start"))
			sraTracks[[ state ]]$setParam("endColumns", c("extended.end", "end", "overlap.end"))
			sraTracks[[ state ]]$setParam("colorVal", NA)
			sraTracks[[ state ]]$setParam("border", "#000000")
			sraTracks[[ state ]]$setParam(
				"colorFun",
				function() {
					grey(1 - slice$value)
				}
			)
		}
		
		# Output
		return(sraTracks)
	} else {
		# Raw output
		return(output)
	}
},

STEPS = function(tracks = TRUE, ...) {
"Apply the STEPS() function to prioritize commonly altered regions.
- tracks   : single logical value, whether to convert output to track.table class or not."
	
	# Argument dispatch
	allArgs = list(...)
	if(any(names(allArgs) == "")) stop("All arguments must be named for dispatch")
	stepsArgs = allArgs[ intersect(names(allArgs), names(formals(cghRA::STEPS))) ]
	
	# Penetrance
	segPenetrance <- .self$penetrance(tracks=FALSE, ...)
	
	# Computation for each state
	output = list()
	for(state in names(segPenetrance)) {
		stepsArgs$segPenetrance <- segPenetrance[[state]]
		output[[state]] <- do.call(
			what = cghRA::STEPS,
			args = stepsArgs
		)
	}
	
	# Conversion
	if(isTRUE(tracks)) {
		if(length(arrays) == 0) stop("Unable to guess annotation from series content (empty series)")
		
		stepsTracks <- list()
		for(state in names(output)) {
			# Track conversion
			output[[ state ]]$strand <- rep("+", nrow(output[[ state ]]))
			output[[ state ]]$start <- output[[ state ]]$start.MCR
			output[[ state ]]$end <- output[[ state ]]$end.MCR
			stepsTracks[[ state ]] <- track.table(
				output[[ state ]],
				.name = sprintf("%s STEPS (%s)", name, state),
				.organism = arrays[[1]]$organism,
				.assembly = arrays[[1]]$assembly,
				.chromosomes = arrays[[1]]$chromosomes(),
				.makeNames = TRUE,
				warn = FALSE
			)
			
			# Graphical parameters
			stepsTracks[[ state ]]$setParam("drawFun", "draw.steps")
			stepsTracks[[ state ]]$setParam("startColumns", c("start.baseline", "start.MCR"))
			stepsTracks[[ state ]]$setParam("endColumns", c("end.baseline", "end.MCR"))
			stepsTracks[[ state ]]$setParam("colorVal", NA)
			stepsTracks[[ state ]]$setParam("border", "#000000")
			stepsTracks[[ state ]]$setParam(
				"colorFun",
				function() {
					grey(pmax(1 - slice$score/10, 0))
				}
			)
		}
		
		# Output
		return(stepsTracks)
	} else {
		# Raw output
		return(output)
	}
}

	)
)

# Constructor
cghRA.series <- function(..., .name, warn=TRUE) {
	# Importation
	object <- new("cghRA.series")
	if(!missing(.name)) object$name <- .name
	
	# Addition
	ignored <- 0L
	elements <- list(...)
	if(length(elements) == 1 && is.list(elements[[1]])) elements <- elements[[1]]
	for(element in elements) {
		if(is(element, "cghRA.regions")) {
			# Correct object
			object$add(element)
		} else if(is.character(element)) {
			# Vector of RDT file names
			for(e in element) {
				if(grepl("\\.rdt$", e, ignore.case=TRUE) && file.exists(e)) {
					tmp <- readRDT(e)
					if(is(tmp, "cghRA.regions")) { object$add(tmp)
					} else                       { ignored <- ignored + 1L; warning("\"", e, "\" is an existing RDT file, but it does not contain a cghRA.regions object")
					}
				} else                            { ignored <- ignored + 1L; warning("\"", e, "\" is not an existing RDT file")
				}
			}
		} else {
			# Other type, not recognized
			ignored <- ignored + 1L			
		}
	}
	
	# Ignored warning
	if(ignored > 0L) warning(ignored, " argument elements were ignored")
	
	# Check
	object$check(warn=warn)
	
	return(object)
}
