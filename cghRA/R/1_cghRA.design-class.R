# Reference class for the design element of cghRA objects
# Author : Sylvain Mareschal <mareschal@ovsa.fr>

# R5 sub-class definition
setRefClass(
	Class = "cghRA.design",
	contains = "track.table",
	methods = list(

bias = function(...) {
"Computes the Waves aCGH Correction Algorithm (Lepretre et al. 2009) bias for the current design.
- ...   : arguments to be passed to the bias() function (except from 'probeChrom', 'probeStarts' and 'probeEnds')."
	
	# Bias computation
	biasResults = cghRA::bias(
		probeChrom = .self$extract(,"chrom"),
		probeStarts = .self$extract(,"start"),
		probeEnds = .self$extract(,"end"),
		...
	)
	
	# Adding bias columns in design
	for(k in dimnames(biasResults)[[2]]) {
		.self$addColumn(biasResults[,k], k)
	}
},

check = function(warn=TRUE) {
"Raises an error if the object is not valid, else returns TRUE"
	
	# refTable, track.table
	callSuper(warn)
	
	# Columns
	if(!"id" %in% .self$getColNames()) stop("Must contain a 'id' column")
	if(.self$types("id") != "integer") stop("'id' column must be of class 'integer'")
	
	return(TRUE)
},


defaultParams = function(...) {
"Returns class-specific defaults for graphical parameters. Inheriting class should overload it to define their own defaults.
- ...   : may be used by inheriting methods, especially for inter-dependant parameters."
	
	params <- callSuper(...)
	
	params$maxElements <- 5000L
	params$label <- FALSE
	
	return(params)
},

initialize = function(...) {
	callSuper(...)
	.self
},

remap = function(...) {
"Recomputes the coordinates of the probes from the probes and genome sequences. Forces 'chrom' to factor, keeping levels if available.
- ...   : arguments to be passed to the localize() function."
	
	# Localization
	probeLoc = cghRA::localize(...)
	
	# Ordering probeLoc according to the design
	probeLoc = probeLoc[ factor(x=.self$extract(,"name"), levels=probeLoc$name, ordered=TRUE) ,]
	
	# Object updating
	.self$fill(, "chrom", factor(probeLoc$chrom, levels=.self$chromosomes()))
	.self$fill(, "strand", factor(probeLoc$strand, levels=c("-","+")))
	.self$fill(, "start", probeLoc$start)
	.self$fill(, "end", probeLoc$end)
	
	# Object reordering
	rowOrder(c("chrom", "start"))
},

show = function(include=FALSE, fieldWidth=10) {
"Interactive printing
- include   : single logical value, if TRUE class name will not be printed."
	
	# Class name
	if(!isTRUE(include)) {
		cat("\n  \"cghRA.design\" reference class object\n")
	} else {
		cat("\n  Extends \"cghRA.design\"\n")
	}
	
	# Inherited show()
	callSuper(include=TRUE, fieldWidth=fieldWidth)
}

	)
)

# Constructor
cghRA.design <- function(..., warn=TRUE) {
	# Inheritance
	object <- new("cghRA.design")
	object$import(track.table(..., warn=FALSE))
	
	# Check
	object$check(warn=warn)
	
	return(object)
}
