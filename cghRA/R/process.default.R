# Default values for process() and tk.process()
# Author : Sylvain Mareschal <mareschal@ovsa.fr>

process.default <- function(argName, profileName) {
	# No argument
	if(missing(argName)) {
		return(
			list(
				segmentArgs = c("fast", "fastest", "CBS default", "accurate"),
				modelizeArgs = c("default")
			)
		)
	}
	
	# Get default
	if(argName == "segmentArgs") {
		if(missing(profileName) || profileName == "fast") {
			# Fast segmentation (default)
			return(
				c(
					"alpha=0.0100, nperm=25/0.0100, undo.splits='sdundo', undo.SD=1",
					"alpha=0.0100, nperm=25/0.0100, undo.splits='sdundo', undo.SD=2"
				)
			)
		} else if(profileName == "fastest") {
			# Fastest segmentation
			return("alpha=0.05, nperm=20/0.05, undo.splits='sdundo', undo.SD=1")
		} else if(profileName == "CBS default") {
			# Default values
			return("")
		}else if(profileName == "accurate") {
			# Most accurate segmentation
			return(
				c(
					"alpha=0.1000, nperm=25/0.1000, undo.splits='sdundo', undo.SD=1",
					"alpha=0.1000, nperm=25/0.1000, undo.splits='sdundo', undo.SD=2",
					"alpha=0.0500, nperm=25/0.0500, undo.splits='sdundo', undo.SD=1",
					"alpha=0.0500, nperm=25/0.0500, undo.splits='sdundo', undo.SD=2",
					"alpha=0.0100, nperm=25/0.0100, undo.splits='sdundo', undo.SD=1",
					"alpha=0.0100, nperm=25/0.0100, undo.splits='sdundo', undo.SD=2",
					"alpha=0.0010, nperm=25/0.0010, undo.splits='sdundo', undo.SD=1",
					"alpha=0.0010, nperm=25/0.0010, undo.splits='sdundo', undo.SD=2",
					"alpha=0.0001, nperm=25/0.0001, undo.splits='sdundo', undo.SD=1",
					"alpha=0.0001, nperm=25/0.0001, undo.splits='sdundo', undo.SD=2"
				)
			)
		} else stop("Unknown 'profileName' for \"", argName, "\" argument")
	} else if(argName == "modelizeArgs") {
		if(missing(profileName) || profileName == "default") {
			# Standard modelization
			return("")
		} else stop("Unknown 'profileName' for \"", argName, "\" argument")
	}
}
