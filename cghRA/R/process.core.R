# Sub-function for process(), processing a single array. Not intended to be called by the user.
# Author : Sylvain Mareschal <mareschal@ovsa.fr>

process.core = function(
		input,
		inputName,
		output,
		steps = c("parse", "mask", "replicates", "waca", "export", "spatial", "segment", "fill", "fixLast", "modelize", "export", "fittest", "export", "applyModel", "export"),
		...
		)
	{
	# Header
	T0 <- proc.time()['elapsed']
	message(
		format(Sys.time(), "[%H:%M:%S] "),
		inputName,
		" : ",
		appendLF = TRUE
	)
	
	# Argument list
	args <- list(
		input = input,
		...	
	)
	
	first <- TRUE
	for(step in steps) {
		# Get function
		fun <- sprintf("process.%s", step)
		if(!exists(fun)) stop("No \"", fun, "\" function found, unknown step")
		fun <- get(fun)
		
		# Apply function
		message(step, ", ", appendLF=FALSE)
		args$input <- do.call(what=fun, args=args)
	}
	
	# Output of the last step
	if(isTRUE(output)) { out <- args$input
	} else             { out <- TRUE
	}
	
	# Cleaning
	message("clean, ", appendLF=FALSE)
	rm(args)
	gc()
	gc()
	
	# Footer
	T1 <- proc.time()['elapsed']
	d <- as.integer(T1 - T0);
	message(sprintf("done (%02i:%02i:%02i)\n", d %/% 3600, d %% 3600 %/% 60, d %% 60))
	
	invisible(out)
}



### PROCESS STEPS ###

process.applyModel <- function(input, ...) {
	# Multiple inputs
	multiple <- is.list(input)
	if(!multiple) input <- list(input)
	
	out <- vector(mode="list", length=length(input))
	for(i in 1:length(input)) {
		# Checks
		if(!is(input[[i]], "cghRA.regions")) stop("'input' must be cghRA.regions objects (#", i, ")")
		
		# Apply model
		out[[i]] <- input[[i]]$model.apply()
	}
	
	if(multiple) { return(out)
	} else       { return(out[[1]])
	}
}

process.cnvScore <- function(input, design, dgv.map, cnvScoreCol="cnvScore", ...) {
	# Checks
	if(!is.character(design) || length(design) != 1 || is.na(design) || !file.exists(design))     stop("'design' must be an existing file name")
	if(!is.character(dgv.map) || length(dgv.map) != 1 || is.na(dgv.map) || !file.exists(dgv.map)) stop("'dgv.map' must be an existing file name")
	if(!grepl("\\.rdt$", design, ignore.case=TRUE))                                               stop("'design' must be a \".rdt\" file")
	if(!grepl("\\.rds$", dgv.map, ignore.case=TRUE))                                              stop("'dgv.map' must be a \".rds\" file")
	
	# Import components
	dgv.map <- readRDS(dgv.map)
	design <- readRDT(design)
	
	# Multiple inputs
	multiple <- is.list(input)
	if(!multiple) input <- list(input)
	
	# Update 'input' directly
	for(i in 1:length(input)) {
		# Checks
		if(!is(input[[i]], "track.table")) stop("'input' must consist of track.table inheriting objects (#", i, ")")
		
		# Remap segments to current design
		obj.map <- map2design(input[[i]], design, quiet=TRUE)
		
		# Compute CNV score
		score <- cnvScore(obj.map, dgv.map, expand=TRUE, quiet=TRUE)
		
		# Store / Replace in table
		if(cnvScoreCol %in% input[[i]]$getColNames()) input[[i]]$delColumns(cnvScoreCol)
		input[[i]]$addColumn(score, cnvScoreCol)
	}
	
	if(multiple) { return(input)
	} else       { return(input[[1]])
	}
}

process.export <- function(input, outDirectory, ...) {
	# Multiple inputs
	multiple <- is.list(input)
	if(!multiple) input <- list(input)
	
	# Collect scores
	stm <- double(length(input))
	for(i in 1:length(input)) {
		# Array export
		if(is(input[[i]], "cghRA.array")) { object <- input[[i]]$probes
		} else                            { object <- input[[i]]
		}
		
		# Checks
		if(!is(object, "refTable")) stop("'input' must consist of refTable inheriting objects (#", i, ")")
		
		# Multi-export
		if(multiple) { suffix <- sprintf("#%i", i)
		} else       { suffix <- ""
		}
		
		# File name	
		if(grepl("^cghRA\\.", class(object))) { fileName <- sprintf("%s/%s%s.%s.rdt", outDirectory, object$name, suffix, sub("^cghRA\\.", "", class(object)))
		} else                                { fileName <- sprintf("%s/%s%s.rdt", outDirectory, object$name, suffix)
		}
		
		# Export
		saveRDT(object, file=fileName)
	}
	
	if(multiple) { return(input)
	} else       { return(input[[1]])
	}
}

process.fill <- function(input, ...) {
	# Multiple inputs
	multiple <- is.list(input)
	if(!multiple) input <- list(input)
	
	for(i in 1:length(input)) {
		# Checks
		if(!is(input[[i]], "cghRA.regions")) stop("'input' must consist of cghRA.regions inheriting objects (#", i, ")")
	
		# Process
		input[[i]]$fillGaps()
	}
	
	if(multiple) { return(input)
	} else       { return(input[[1]])
	}
}

process.filter <- function(input, filter=NULL, ...) {
	# Checks
	if(!is.expression(filter)) stop("'filter' must be an R expression")
	
	# Multiple inputs
	multiple <- is.list(input)
	if(!multiple) input <- list(input)
	
	out <- vector(mode="list", length=length(input))
	for(i in 1:length(input)) {
		# Checks
		if(!is(input[[i]], "refTable")) stop("'input' must consist of refTable inheriting objects (#", i, ")")
		
		# Extract
		out[[i]] <- input[[i]]$extract(filter, asObject=TRUE)
	}
	
	if(multiple) { return(out)
	} else       { return(out[[1]])
	}
}

process.fittest <- function(input, ...) {
	# Multiple inputs
	multiple <- is.list(input)
	if(!multiple) input <- list(input)
	
	# Collect scores
	stm <- double(length(input))
	for(i in 1:length(input)) {
		# Checks
		if(!is(input[[i]], "cghRA.regions")) stop("'input' must consist of cghRA.regions inheriting objects (#", i, ")")
		
		# Extract STM score
		stm[i] <- input[[i]]$model['stm']
	}
	
	# Return best model
	best <- which.min(stm)
	if(length(best) == 1) { out <- input[[ best ]]
	} else                { stop("No fittest segmentation found")
	}
	
	return(out)
}

process.fixLast <- function(input, design, ...) {
	# Checks
	if(!is.character(design) || length(design) != 1 || is.na(design) || !file.exists(design)) stop("'design' must be an existing file name")
	if(!grepl("\\.rdt$", design, ignore.case=TRUE))                                           stop("'design' must be a \".rdt\" file")
	
	# Import design
	design <- readRDT(design)
	
	# Multiple inputs
	multiple <- is.list(input)
	if(!multiple) input <- list(input)
	
	for(i in 1:length(input)) {
		# Checks
		if(!is(input[[i]], "cghRA.regions")) stop("'input' must consist of cghRA.regions inheriting objects (#", i, ")")
	
		# Process
		input[[i]]$fixLast(design)
	}
	
	if(multiple) { return(input)
	} else       { return(input[[1]])
	}
}

process.mask <- function(input, ...) {
	# Checks
	if(!is(input, "cghRA.array")) stop("'input' must be a cghRA.array inheriting object")
	
	# Process
	input$maskByFlag(
		flags = "^flag_",
		pattern = TRUE,
		multiple = "any"
	)
	
	return(input)
}

process.modelize <- function(input, modelizeArgs=process.default("modelizeArgs"), ...) {
	# Multiple inputs
	multiple <- is.list(input)
	if(!multiple) input <- list(input)
	
	out <- vector(mode="list", length=length(input))
	for(i in 1:length(input)) {
		# Checks
		if(!is(input[[i]], "cghRA.regions")) stop("'input' must consist of cghRA.regions inheriting objects (#", i, ")")
		
		# Build and check call
		regions <- input[[i]]
		modelizeCall <- try(parse(text=sprintf("regions$model.auto(discreet=TRUE, %s)", modelizeArgs)), silent=TRUE)
		if(is(modelizeCall, "try-error")) stop("Parse error in cghRA.copies parameters (#", i, ")")
		
		# Modelization call
		eval(modelizeCall)
		regions$modelizeCall <- modelizeCall[[1]]
		
		# Store
		out[[i]] <- regions
	}
	
	if(multiple) { return(out)
	} else       { return(out[[1]])
	}
}

process.parse <- function(input, design, probeParser=Agilent.probes, probeArgs=list(), ...) {
	# Checks
	if(!is.character(input) || length(input) != 1 || is.na(input) || !file.exists(input))     stop("'input' must be an existing file name")
	if(!is.character(design) || length(design) != 1 || is.na(design) || !file.exists(design)) stop("'design' must be an existing file name")
	if(!grepl("\\.rdt$", design, ignore.case=TRUE))                                           stop("'design' must be a \".rdt\" file")
	
	# Parse
	probeArgs$file <- input
	probes <- do.call(what=probeParser, args=probeArgs)
	
	# Parse design
	designObject <- readRDT(design)
	
	# Assemble
	arr <- cghRA.array(
		.design = designObject,
		.probes = probes,
		warn = FALSE
	)
	
	return(arr)
}

process.probes <- function(input, design, ...) {
	# Checks
	if(!is.character(input) || length(input) != 1 || is.na(input) || !file.exists(input))     stop("'input' must be an existing file name")
	if(!is.character(design) || length(design) != 1 || is.na(design) || !file.exists(design)) stop("'design' must be an existing file name")
	if(!grepl("\\.rdt$", input, ignore.case=TRUE))                                            stop("'input' must be a \".rdt\" file")
	if(!grepl("\\.rdt$", design, ignore.case=TRUE))                                           stop("'design' must be a \".rdt\" file")
	
	# Import components
	probes <- readRDT(input)
	design <- readRDT(design)
	
	# Assemble
	arr <- cghRA.array(
		.design = design,
		.probes = probes,
		warn = FALSE
	)
	
	return(arr)
}

process.regions <- function(input, ...) {
	# Checks
	if(!is.character(input) || any(is.na(input)) || any(!file.exists(input))) stop("'input' must be one or many existing file name(s)")
	if(any(!grepl("\\.rdt$", input, ignore.case=TRUE)))                       stop("'input' must be one or many \".rdt\" file(s)")
	
	# Import probe file
	out <- vector(mode="list", length=length(input))
	for(i in 1:length(input)) out[[i]] <- readRDT(input[i])
	
	if(length(out) > 1) { return(out)
	} else              { return(out[[1]])
	}
}

process.replicates <- function(input, replicateFun=stats::median, ...) {
	# Checks
	if(!is(input, "cghRA.array")) stop("'input' must be a cghRA.array inheriting object")
	
	# Process
	input$replicates(
		fun = replicateFun,
		na.rm = TRUE
	)
	
	return(input)
}

process.segment <- function(input, segmentArgs=process.default("segmentArgs"), ...) {
	# Checks
	if(!is(input, "cghRA.array")) stop("'input' must be a cghRA.array inheriting object")
	
	# Process
	array <- input
	out <- vector(mode="list", length=length(segmentArgs))
	for(i in 1:length(segmentArgs)) {
		# Build and check call
		segmentCall <- try(parse(text=sprintf("regions <- array$DNAcopy(verbose=0, %s)", segmentArgs[i])), silent=TRUE)
		if(is(segmentCall, "try-error")) stop(sprintf("Parse error in DNAcopy profile #%i", i))
		
		# Apply call
		eval(segmentCall)
		regions$segmentCall <- segmentCall[[1]][[3]]
		
		# Store
		out[[i]] <- regions
	}
	
	if(length(out) > 1) { return(out)
	} else              { return(out[[1]])
	}
}

process.spatial <- function(input, outDirectory, ...) {
	# Checks
	if(!is(input, "cghRA.array")) stop("'input' must be a cghRA.array inheriting object")
	
	# Process
	input$spatial(filename=sprintf("%s/%s.spatial.png", outDirectory, input$name))
	
	return(input)
}

process.waca <- function(input, ...) {
	# Checks
	if(!is(input, "cghRA.array")) stop("'input' must be a cghRA.array inheriting object")
	
	# Process
	input$WACA()
	
	return(input)
}

