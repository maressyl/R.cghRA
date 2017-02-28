# Filter and update a copy of 'track' using a map as produced by "map2design" with the same tracks
# Author : Sylvain Mareschal <mareschal@ovsa.fr>

applyMap <- function(track, map, design) {
	# Checks
	if(!inherits(track, "track.table"))   stop("'track' must inherit 'track.table' class")
	if(!inherits(map, "segmentMap"))      stop("'map' must inherit 'segmentMap' class")
	if(!inherits(design, "cghRA.design")) stop("'design' must inherit 'cghRA.design' class")
	
	# Copy dataset
	newTrack <- track$copy()
	
	# Limit to CNVs retained in map (visible), order as 'map'
	indexes <- as.integer(unlist(strsplit(rownames(map$map), split=",")))
	newTrack$rowOrder(indexes)
	
	# Remap to included probe locations
	newStarts <- design$extract(rep(map$map[,"start"], map$map[,"count"]), "start")
	newEnds <- design$extract(rep(map$map[,"end"], map$map[,"count"]), "end")
	newTrack$fill(, "start", newStarts)
	newTrack$fill(, "end", newEnds)
	
	# Keep original index
	newTrack$addColumn(indexes, "index")
	
	# Genomic order
	newTrack$rowOrder(c("chrom", "start"), na.last=TRUE, decreasing=FALSE)
	
	return(newTrack)
}

