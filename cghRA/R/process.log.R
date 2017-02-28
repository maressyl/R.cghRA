# Sub-function for process(), processing a single array. Not intended to be called by the user.
# Author : Sylvain Mareschal <mareschal@ovsa.fr>

process.log = function(
		...,
		logFile
		)
	{
	# Extensive requiring, for sub-processes
	require(cghRA, quietly=TRUE)
	
	if(!is.na(logFile) && file.exists(logFile)) {
		# Initialize
		logLine <- ""
		opt <- options(warn=1)
		on.exit(options(opt))
	
		# Catch error, warnings and messages
		handle(
			expr = {
				process.core(...)
			},
			messageHandler = function(m) {
				logLine <<- sprintf("%s%s", logLine, conditionMessage(m))
			},
			warningHandler = function(w) {
				logLine <<- sprintf("%s | WARNING : %s | ", logLine, conditionMessage(w))
			},
			errorHandler = function(e) {
				logLine <<- sprintf("%s | ERROR : %s\n\n", logLine, conditionMessage(e))
			}
		)
	
		# Write log
		cat(logLine, file=logFile, append=TRUE)
	} else {
		# No logging
		process.core(...)
	}
	
	invisible(TRUE)
}
