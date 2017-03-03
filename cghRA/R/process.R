# Command line array processing
# Author : Sylvain Mareschal <mareschal@ovsa.fr>

process = function(
		inputs,
		logFile = "process.log",
		cluster = NA,
		...
		)
	{
	if(!identical(cluster, FALSE)) {
		if(!requireNamespace("parallel")) {
			warning("Install the 'parallel' package to enhance computation speed by using multiple CPUs")
			cluster <- FALSE
		} else if(is.na(cluster)) {
			# Default cluster definition
			cluster <- list(spec=parallel::detectCores())
		}
	}
	
	# Log file
	if(is.na(logFile)) { logFile <- ""
	} else             { unlink(logFile)
	}
	
	# Log header
	cat(sprintf("%s\n\n", Sys.time()), file=logFile, append=FALSE)
	cat(sprintf("R %s.%s (%s)\n", R.version$major, R.version$minor, R.version$platform), file=logFile, append=TRUE)
	for(package in c("Rgb", "cghRA", "DNAcopy")) {
		if(length(find.package(package, quiet=TRUE)) > 0) cat(sprintf("%-10s version %s\n", package, utils::packageVersion(package)), file=logFile, append=TRUE)
	}
	cat("\n\n", file=logFile, append=TRUE)
	
	# Input names
	if(!is.null(names(inputs))) {
		# Provided
		inputNames <- names(inputs)
	} else if(is.list(inputs)) {
		# List
		if(all(sapply(inputs, is.atomic))) {
			# Consider the first elements
			inputNames <- sapply(inputs, "[[", 1)
			if(is.character(inputNames)) inputNames <- basename(inputNames)
			inputNames <- sprintf("%s [1:%i]", inputNames, sapply(inputs, length))
		} else {
			# Mixed or non-atomic list
			inputNames <- sprintf("#%i", 1:length(inputs))
		}
	} else if(is.character(inputs)) {
		# Probably a vector of file names
		inputNames <- basename(inputs)
	} else if(is.atomic(inputs)) {
		# As is
		inputNames <- inputs
	} else {
		# Single non-atomic object ?
		inputNames <- sprintf("#%i", 1:length(inputs))
	}
	
	if(identical(cluster, FALSE)) {
		# Linear local processing
		for(i in 1:length(inputs)) {
			process.log(
				input = inputs[[i]],
				inputName = inputNames[i],
				logFile = logFile,
				...
			)
		}
	} else {
		# Cluster processing
		cluster <- do.call(what=parallel::makeCluster, args=cluster)
		on.exit(parallel::stopCluster(cluster), add=TRUE)
		parallel::clusterMap(
			cl = cluster,
			fun = process.log,
			input = inputs,
			inputName = inputNames,
			MoreArgs = list(
				logFile = logFile,
				...
			),
			RECYCLE = FALSE,
			SIMPLIFY = FALSE,
			USE.NAMES = FALSE,
			.scheduling = "dynamic"
		)
	}
	
	# Finalize
	lastMessage <- format(Sys.time(), "[%H:%M:%S] All done\n")
	cat(lastMessage, file=logFile, append=TRUE)
	
	invisible(TRUE)
}
