# Novel implementation of the Wave A-CGH Correction Algorithm (Lepretre et al NAR 2010, PMID:20071741)
# Author : Sylvain Mareschal <mareschal@ovsa.fr>

WACA = function(
		probeNames,
		probeLogRatios,
		bias,
		forceBiasOrdering = TRUE
		)
	{
	# Checks
	probeNames = as.character(probeNames)
	probeLogRatios = as.double(probeLogRatios)
	if(length(probeNames) != length(probeLogRatios)) {
		stop("'probeNames' and 'probeLogRatios' lengths differ")
	}
	if(!is.data.frame(bias) || any(dim(bias) == 0)) {
		stop("'bias' must be a non empty 'data.frame'")
	}
	if(!identical(names(bias), c("wGC150", "wGC500", "wGCprobe", "wGCfrag", "wFragSize"))) {
		warning("Genuine WACA needs 'bias' columns to be 'wGC150', 'wGC500', 'wGCprobe', 'wGCfrag' and 'wFragSize'")
	}
	# if(any(sapply(bias, typeof) != "integer")) {
	#	warning("Integer 'bias' values should be prefered")
	# }
	
	# Ordering 'bias' according to 'probeNames'
	if (isTRUE(forceBiasOrdering) || nrow(bias) != length(probeNames)) {
		if(any(!probeNames %in% row.names(bias))) {
			stop("All 'probeNames' must be present in 'bias' row names")
		}
		bias = bias[ factor(x=probeNames, levels=row.names(bias), ordered=TRUE) , , drop=FALSE ]
	}
	
	for(biasName in names(bias)) {
		# Consider only non NA logRatios and biases
		selection = !is.na(probeLogRatios) & !is.na(bias[, biasName])
		
		# LOWESS fit and simplification 
		model = stats::lowess(x=bias[selection, biasName], y=probeLogRatios[selection], f=2/3)
		model = as.data.frame(model)
		model = model[ !duplicated(model$x) ,]
		
		# Substract bias from each selected probe ('float' bias may be problematic here)
		probeLogRatios[selection] = probeLogRatios[selection] - model$y[ factor(x=bias[selection, biasName], levels=model$x, ordered=TRUE) ]
		
		rm(model)
	}
	
	return(probeLogRatios)
}
