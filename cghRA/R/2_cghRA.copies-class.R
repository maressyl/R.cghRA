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
