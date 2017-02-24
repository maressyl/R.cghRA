# Builds a data.frame with a column for each 'groups' with according 'values' in them, and a distinct row for each section
# Typically used to merge logRatio (values) of segments (starts : ends) from distinct samples (groups)
# Author : Sylvain Mareschal <mareschal@ovsa.fr>

segPile = function(
		starts,
		ends,
		values,
		groups,
		groupList = unique(groups)
		)
	{
	# All breakpoints
	breaks = c(starts, ends)
	breaks = unique(breaks)
	breaks = breaks[ order(breaks) ]
	
	# Breakpoint index
	index = 1:length(breaks)
	names(index) = breaks

	# To characters
	starts = as.character(starts)
	ends = as.character(ends)
	groups = as.character(groups)
	
	# Mapping matrix
	map = matrix(
		ncol = length(groupList),
		nrow = length(breaks) - 1,
		dimnames = list(
			NULL,
			group = groupList
		)
	)
	
	# Mapping
	for(r in 1:length(starts)) {
		if(starts[r] != ends[r]) {
			map[ index[starts[r]] : (index[ends[r]]-1) , groups[r] ] = values[r]
		}
	}
	
	# Reshapping
	tab = data.frame(
		start = breaks[ 1:(length(breaks)-1) ],
		end = breaks[ 2:length(breaks) ],
		map,
		stringsAsFactors = FALSE,
		check.names = FALSE
	)
	
	return(tab)
}
