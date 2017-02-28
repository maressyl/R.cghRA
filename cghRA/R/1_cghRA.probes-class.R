# Reference class for the probes element of cghRA.array objects
# Author : Sylvain Mareschal <mareschal@ovsa.fr>

# R5 sub-class definition
setRefClass(
	Class = "cghRA.probes",
	contains = "refTable",
	fields = list(
		name = "character"
	),
	methods = list(

check = function(warn=TRUE) {
"Raises an error if the object is not valid, else returns TRUE"
	
	# refTable
	callSuper(warn)
	
	# Specificities
	if(isTRUE(rowNamed))   stop("Must not be rowNamed")
	
	# Fields
	if(length(name) != 1)  stop("'name' must be a single character value")
	
	# Warnings
	if(isTRUE(warn)) {
		if(is.na(name))   warning("'name' should not be NA")
		if(rowCount == 0) warning("The table is empty")
	}
	
	# Columns
	if(!"id" %in% .self$getColNames())       stop("Must contain a 'id' column")
	if(!"logRatio" %in% .self$getColNames()) stop("Must contain a 'logRatio' column")
	if(.self$types("id") != "integer")       stop("'id' column must be of class 'integer'")
	if(.self$types("logRatio") != "numeric") stop("'logRatio' column must be of class 'numeric'")
	
	return(TRUE)
},

initialize = function(name=NA_character_, ...) {
	callSuper(...)
	initFields(name=name)
},

show = function(include=FALSE, fieldWidth=10) {
"Interactive printing
- include   : single logical value, if TRUE class name will not be printed."
	
	# Class name
	if(!isTRUE(include)) { cat("\n  \"cghRA.probes\" reference class object\n")
	} else               { cat("\n  Extends \"cghRA.probes\"\n")
	}
	
	# Fields
	cat(sprintf("  %-*s : %s\n", fieldWidth, "name", name[1]))
	
	# Inherited show()
	callSuper(include=TRUE, fieldWidth=fieldWidth)
}		

	)
)

# Constructor
cghRA.probes <- function(..., .name, warn=TRUE) {
	# Inheritance
	object <- new("cghRA.probes")
	if(!missing(.name)) object$name <- .name
	object$import(refTable(..., warn=warn))
	
	# Check
	object$check(warn=warn)
	
	return(object)
}

