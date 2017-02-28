# Remap a track to what can be expected from a given CGH-array design
# Author : Sylvain Mareschal <mareschal@ovsa.fr>

map2design <- function(track, design, minProbes=1, quiet=FALSE, warn=TRUE) {
	# Checks
	if(!inherits(track, "track.table"))                                   stop("'track' must inherit 'track.table' class")
	if(!inherits(design, "cghRA.design"))                                 stop("'design' must inherit 'cghRA.design' class")
	if(length(intersect(design$chromosomes(), track$chromosomes())) == 0) stop("'track' and 'design' chromosomes do not overlap")
	
	# Grouping identical segments (bp coordinates)
	if(!isTRUE(quiet)) message("Grouping segments with identical original boundaries ...")
	groupIDs <- sprintf("%s:%i-%i", track$extract(,"chrom"), track$extract(,"start"), track$extract(,"end"))
	groups <- tapply(X=1:track$rowCount, INDEX=groupIDs, FUN=c)
	groupData <- track$extract(sapply(groups, "[", 1L), c("chrom", "start", "end"))
	
	# Faster access to columns
	gChrom <- groupData$chrom
	gStart <- groupData$start
	gEnd <- groupData$end
	gIndexes <- sapply(groups, paste, collapse=",")
	
	# Simplified design
	if(!isTRUE(quiet)) message("Simplifying design ...")
	newDesign <- design$copy()
	newDesign$delColumns(setdiff(design$colNames, c("chrom", "start", "end")))
	newDesign$addColumn(1:design$rowCount, "id")
	
	# To store new coordinates
	mtx <- matrix(as.integer(NA), nrow=length(gIndexes), ncol=2, dimnames=list(gIndexes, c("start", "end")))
	
	# Loop over CNV
	if(!isTRUE(quiet)) message("Remapping segments ...")
	if(!isTRUE(quiet)) message(0, "/", nrow(mtx))
	for(i in 1:nrow(mtx)) {
		if(!isTRUE(quiet) && i %% 5000 == 0) message(i, "/", nrow(mtx))
		probes <- newDesign$slice(chrom=gChrom[i], start=gStart[i], end=gEnd[i])$id
		if(length(probes) > 0) {
			mtx[i,"start"] <- probes[ 1L ]
			mtx[i,"end"] <- probes[ length(probes) ]
		}
	}
	if(!isTRUE(quiet)) message(i, "/", nrow(mtx))
	
	# Filter on visible CNV
	if(!isTRUE(quiet)) message("Filtering segments ...")
	n <- mtx[,"end"] - mtx[,"start"] + 1
	mtx <- mtx[ !is.na(n) & n >= minProbes , , drop=FALSE ]
	
	# Grouping identical segments, keeping original 'track' indexes
	if(!isTRUE(quiet)) message("Grouping segments with identical remapped boundaries ...")
	indexes <- tapply(X=rownames(mtx), INDEX=sprintf("%i-%i", mtx[,"start"], mtx[,"end"]), FUN=paste, collapse=",")
	
	# Reshaping into a matrix, keeping original 'track' index list as row names
	members <- as.character(indexes)
	boundaries <- strsplit(names(indexes), split="-")
	map <- matrix(as.integer(NA), nrow=length(indexes), ncol=3L, dimnames=list(members, c("start", "end", "count")))
	map[,"start"] <- as.integer(sapply(boundaries, "[", 1L))
	map[,"end"] <- as.integer(sapply(boundaries, "[", 2L))
	
	# Segment count
	count <- sapply(gregexpr(",", members, fixed=TRUE), length) + 1L
	count[ grep(",", members, invert=TRUE, fixed=TRUE) ] <- 1L
	map[,"count"] <- as.integer(count)
	
	# Genomic order
	map <- map[ order(map[,"start"], map[,"end"]) , , drop=FALSE ]
	
	# Create an object
	obj <- new("segmentMap")
	obj$map <- map
	obj$trackName <- track$name
	obj$trackSize <- track$rowCount
	obj$designName <- design$name
	obj$designSize <- design$rowCount
	obj$check(warn)
	
	return(obj)
}

