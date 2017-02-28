# Fill gaps between segments, extending on the right side
# Author : Sylvain Mareschal <mareschal@ovsa.fr>

fillGaps = function(
		segTable,
		isOrdered = FALSE
		)
	{
	# Arg checks
	if(!is.data.frame(segTable))                              stop("'segTable' must be a data.frame")
	if(any(!c("chrom", "start", "end") %in% names(segTable))) stop("'segTable' must have 'chrom', 'start' and 'end' columns")
	
	# Ordering
	if(!isTRUE(isOrdered)) {
		segTable = segTable[ order(segTable$chrom, segTable$start) ,]
	}
	
	# Last index of each chromosome
	breaks = c(
		0,
		tapply(
			X = 1:nrow(segTable),
			INDEX = segTable$chrom,
			FUN = max
		)
	)   
	
	# All indexes without first ones of each chromosome
	starts = setdiff(1:nrow(segTable), breaks+1)
	
	# All indexes without last ones of each chromosome
	stops = setdiff(1:nrow(segTable), breaks)
	
	# New end : the following start
	if (length(starts) > 0) segTable[ stops , "end" ] = segTable[ starts , "start" ]
	
	return(segTable)
}
