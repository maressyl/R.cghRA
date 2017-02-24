# Reference class for a segment list located by probe IDs rather than genomic position
# Author : Sylvain Mareschal <mareschal@ovsa.fr>

# R5 sub-class definition
setRefClass(
	Class = "segmentMap",
	fields = list(
		map = "matrix",
		trackName = "character",
		trackSize = "integer",
		designName = "character",
		designSize = "integer"
	),
	methods = list(

check = function(warn=TRUE) {
"Raises an error if the object is not valid, else returns TRUE"
	
	# Fields
	if(storage.mode(map) != "integer")                        stop("'map' is supposed to be an integer matrix")
	if(!identical(colnames(map), c("start", "end", "count"))) stop("'map' columns must be 'start', 'end' and 'count'")
	if(any(is.na(map)))                                       stop("'map' must not contain NA values")
	if(any(map < 0L))                                         stop("'map' must not contain negative values")
	if(any(!grepl("^[0-9,]+$", rownames(map))))               stop("'map' row names must be comma-separated integer values")
	if(length(trackName) != 1L)                               stop("'trackName' must contain a single value")
	if(length(trackSize) != 1L)                               stop("'trackSize' must contain a single value")
	if(length(designName) != 1L)                              stop("'designName' must contain a single value")
	if(length(designSize) != 1L)                              stop("'designSize' must contain a single value")
	if(is.na(trackSize))                                      stop("'trackSize' must not be NA")
			
	# Warnings
	if(isTRUE(warn)) {
		if(is.na(trackName))                                  warning("'trackName' should not be NA")
		if(is.na(designName))                                 warning("'designName' should not be NA")
		if(is.na(designSize))                                 warning("'designSize' should not be NA")
	}
	
	return(TRUE)
},

initialize = function(
		map = matrix(as.integer(NA), nrow=0L, ncol=3L, dimnames=list(NULL, c("start", "end", "count"))),
		trackName = as.character(NA),
		trackSize = as.integer(NA),
		designName = as.character(NA),
		designSize = as.integer(NA),
		...
	) {
	# Initialize an empty table
	.self$initFields(
		map = map,
		trackName = trackName,
		trackSize = trackSize,
		designName = designName,
		designSize = designSize
	)
	.self
},

show = function(include=FALSE, fieldWidth=10) {
"Interactive printing
- include   : single logical value, if TRUE class name will not be printed."
	
	# Class name
	if(!isTRUE(include)) {
		cat("\n  \"segmentMap\" reference class object\n")
	} else {
		cat("\n  Extends \"segmentMap\"\n")
	}
}

	)
)

