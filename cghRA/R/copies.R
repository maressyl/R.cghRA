# Convert log-ratios into Log-ratio related Copy Numbers
# Author : Sylvain Mareschal <mareschal@ovsa.fr>

LCN = function(x, exact=TRUE) {
	out <- 2*(2^x)
	if(!isTRUE(exact)) out <- round(out)
	return(out)
}


# Compute copy numbers from various log-ratio related input types
# Author : Sylvain Mareschal <mareschal@ovsa.fr>

copies = function(
		x,
		model = NA,
		center = model['center'],
		width = model['width'],
		ploidy = model['ploidy'],
		exact = TRUE,
		from = c("logRatios", "LCN", "copies")
		)
	{
	# Checks
	from <- match.arg(from)
	if(length(ploidy) != 1 || is.na(ploidy) || !is.numeric(ploidy)) stop("'ploidy' must be a single non NA numeric value")
	
	# Convert back copies into logRatios
	if(from == "copies") {
		if(ploidy == 0) {
			ploidy <- 2
			warning("logRatios produced from relative copy numbers assuming an original ploidy of 2")
		}
		return(log(x/ploidy, 2))
	}
	
	# Checks
	if(length(center) != 1 || is.na(center) || !is.numeric(center)) stop("'center' must be a single non NA numeric value")
	if(length(width) != 1 || is.na(width) || !is.numeric(width))    stop("'width' must be a single non NA numeric value")
	
	# Starting from LCN or logRatios
	if(from == "LCN") out <- x
	else              out <- LCN(x, exact=TRUE)
	
	# Apply model
	out <- ((out - center) / width) + ploidy
	
	# Rounding
	if(!isTRUE(exact)) out <- round(out)
	
	return(out)
}
