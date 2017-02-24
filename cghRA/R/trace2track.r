# Convert a path trace as returned by cnvScore(trace=TRUE) into an Rgb track.table
# Author : Sylvain Mareschal <mareschal@ovsa.fr>

trace2track <- function(paths, dgv.map, dgv.track) {
	# Checks
	if(!is.data.frame(paths))                                     stop("'paths' must be a data.frame, as returned by cnvScore(trace=TRUE)")
	if(!is.matrix(dgv.map) || storage.mode(dgv.map) != "integer") stop("'dgv.map' must be an integer matrix, as returned by map2design()")
	if(!identical(colnames(dgv.map), c("start", "end", "count"))) stop("'dgv.map' columns must be 'start', 'end' and 'count'")
	if(!inherits(dgv.track, "track.table"))                       stop("'dgv.track' must inherit 'track.table' class")
	
	if(nrow(paths) > 0) {
		# Get dgv.map indexes of mapped CNVs
		mapIndexes <- strsplit(paths$path.cnvList, split=", ")
	
		# Get dgv indexes of mapped CNVs (first index of similar CNVs)
		dgvIndexes <- lapply(mapIndexes, function(x, dgvMap) { as.integer(sub(",.+$", "", rownames(dgv.map)[ as.integer(x) ])) }, dgvMap=dgv.map)
	
		# Get dgv.track indexes of mapped CNVs (first index of similar CNVs)
		trackIndexes <- lapply(dgvIndexes, match, table=dgv.track$extract(,"index"))
	
		# Labels (must be unique)
		pathIndexes <- rep(1:length(trackIndexes), sapply(trackIndexes, length))
		labels <- with(paths[pathIndexes,], sprintf("#%i : %g x %i", pathIndexes, signif(path.score, 3), path.count))
	
		# Build custom track
		track <- dgv.track$copy()
		track$rowOrder(unlist(trackIndexes))
		track$addColumn(labels, "path.score")
		track$addColumn(unlist(lapply(trackIndexes, function(x){1:length(x)})), "path.pos")
		track$addColumn(unlist(lapply(trackIndexes, function(x){rep(length(x), length(x))})), "path.cnvCount")
		
		# Retrieve genomic order
		track$rowOrder(c("chrom", "start", "end"))
	} else {
		# Empty track
		track <- dgv.track$extract(FALSE, asObject=TRUE)
		track$addColumn(character(0), "path.score")
		track$addColumn(integer(0), "path.pos")
		track$addColumn(integer(0), "path.cnvCount")
	}
	
	# Drawing parameters
	track$setParam("groupBy", "path.score")
	track$setParam("grouPosition", "path.pos")
	track$setParam("groupSize", "path.cnvCount")
	track$setParam("border", "#000000")
	track$setParam("label", TRUE)
	track$setParam("labelAdj", "center")
	track$setParam("colorFun",
		function(){
		    output <- rep("#88888899", nrow(slice))
		    output[slice$type %in% c("Duplication", "Gain", "Insertion")] <- "#8888FF99"
		    output[slice$type == "Deletion"] <- "#FF888899"
		    output[slice$type == "Loss"] <- "#88000099"
		    output[slice$type == "Complex"] <- "#BB774499"
		    return(output)
		}
	)
	track$name <- "CNV paths"
	
	return(track)
}

