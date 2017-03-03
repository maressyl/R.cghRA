# Implements the "Gene Expression and Dosage Integrator" CGH/transcriptome correlation, as described by Lenz et al. in PNAS 2008 (Supplemental Informations)
# Author : Sylvain Mareschal <mareschal@ovsa.fr>

GEDI <- function(
		cgh,
		cgh.chrom,
		cgh.start,
		cgh.end,
		cgh.genes,
		expr,
		expr.genes,
		permutations = 1000,
		type = c("amplifications", "deletions"),
		quiet = FALSE
	) {
	# Type checks
	if(!is.matrix(cgh) | !is.logical(cgh))                                           stop("'cgh' must be a logical matrix")
	if(!is.character(cgh.chrom) | length(cgh.chrom) == 0 | any(is.na(cgh.chrom)))    stop("'cgh.chrom' must be a non-NA non-empty character vector")
	if(!is.integer(cgh.start) | length(cgh.start) == 0 | any(is.na(cgh.start)))      stop("'cgh.start' must be a non-NA non-empty integer vector")
	if(!is.integer(cgh.end) | length(cgh.end) == 0 | any(is.na(cgh.end)))            stop("'cgh.end' must be a non-NA non-empty integer vector")
	if(!is.character(cgh.genes) | length(cgh.genes) == 0 | any(is.na(cgh.genes)))    stop("'cgh.genes' must be a non-NA non-empty character vector")
	if(!is.matrix(expr) | !is.numeric(expr))                                         stop("'expr' must be a numeric matrix")
	if(!is.character(expr.genes) | length(expr.genes) == 0 | any(is.na(expr.genes))) stop("'expr.genes' must be a non-NA non-empty character vector")
	
	# Alteration type
	type <- match.arg(type)
	if(type == "amplifications") alternative <- "greater"
	else                         alternative <- "less"
	
	# Consistency checks
	if(length(unique(c(length(cgh.chrom), length(cgh.start), length(cgh.end), length(cgh.genes), nrow(cgh)))) > 1) stop("'cgh.chrom', 'cgh.start, 'cgh.end', 'cgh.genes' lengths and 'cgh' row count must unique")
	if(nrow(expr) != length(expr.genes))                                                                           stop("'expr.genes' length and 'expr' row count must be equals")
	
	# Samples
	samples <- intersect(colnames(cgh), colnames(expr))
	if(length(samples) == 0)                                         stop("'cgh' and 'expr' colnames must intersect")
	if(length(samples) != ncol(cgh) | length(samples) != ncol(expr)) warning("Some 'cgh' and / or 'expr' columns are ignored as not present in each dataset")
	if(!quiet) message("Computing on ", length(samples), " samples :")
	
	# Gene lists splitting
	cgh.genes <- strsplit(cgh.genes, split=", ")
	if(all(sapply(cgh.genes, length) <= 1)) warning("Are you sure 'cgh.genes' genes are separated by ', ' ?")
	
	# Output
	gediScore <- double(nrow(cgh))
	gediGenes <- character(nrow(cgh))
	
	# For each region
	for(region in 1:nrow(cgh)) {
		if(!quiet) message("MCR ", region, " / ", nrow(cgh))
		
		# Common genes
		probeSelection <- which(expr.genes %in% cgh.genes[[region]])
		gediGenes[region] <- paste(unique(expr.genes[ probeSelection ]), collapse=", ")
		
		# Sample groups
		altr <- names(which(cgh[region, samples]))
		germ <- names(which(!cgh[region, samples]))
		
		if(length(probeSelection) > 0) {
			if(length(altr) > 1 && length(germ) > 1) {
				# t-test for each involved probeset
				pvalues <- t.tests(x=expr[ probeSelection , altr , drop=FALSE ], y=expr[ probeSelection , germ , drop=FALSE ], alternative=alternative)[,"pval"]
				
				# Best p-value for each gene
				pvalues <- tapply(X=pvalues, INDEX=expr.genes[ probeSelection ], FUN=min, na.rm=TRUE)
				
				# True score computation
				trueScore <- sum(-log(pvalues))
				
				# Permutated scores
				permScores <- double(permutations)
				for(i in 1:permutations) {
					# Group redefinition after permutation
					perm <- cgh[region, samples]
					names(perm) <- sample(names(perm))
					altr <- names(which(perm))
					germ <- names(which(!perm))
					
					# t-test for each involved probeset
					pvalues <- t.tests(x=expr[ probeSelection , altr , drop=FALSE ], y=expr[ probeSelection , germ , drop=FALSE ], alternative=alternative)[,"pval"]
				
					# Best p-value for each gene
					pvalues <- tapply(X=pvalues, INDEX=expr.genes[ probeSelection ], FUN=min, na.rm=TRUE)
					
					# True score computation
					permScores[i] <- sum(-log(pvalues))
				}
				
				# Final score
				gediScore[region] <- sum(trueScore > permScores) / permutations
			} else {
				# No alteration
				gediScore[region] <- NA
			}
		} else {
			# No common gene
			gediScore[region] <- NA
		}
	}
	
	return(list(gediScore=gediScore, gediGenes=gediGenes))
}

# Parallel t-test, for numeric matrixes
t.tests <- function(x, y, MARGIN=1, alternative=c("two.sided", "less", "greater"), mu=0, paired=FALSE, var.equal=FALSE) {
	# Checks
	if(!is.matrix(x) | !is.numeric(x)) stop("'x' must be a numeric matrix")
	if(!is.matrix(y) | !is.numeric(y)) stop("'x' must be a numeric matrix")
	if(mu != 0)                        stop("mu != 0 not implemented")
	if(paired)                         stop("paired = TRUE not implemented")
	if(var.equal)                      stop("var.equal = TRUE not implemented")
	
	# Probe means
	sx <- apply(x, MARGIN, sum, na.rm=TRUE)
	sy <- apply(y, MARGIN, sum, na.rm=TRUE)
	nx <- apply(!is.na(x), MARGIN, sum)
	ny <- apply(!is.na(y), MARGIN, sum)
	mx <- sx / nx
	my <- sy / ny
	
	# T statistics
	vx <- apply(x, MARGIN, stats::var, na.rm=TRUE)
	vy <- apply(y, MARGIN, stats::var, na.rm=TRUE)
	tstat <- (mx - my) / sqrt(vx/nx + vy/ny)
	
	# Degrees of freedom
	df <- ((vx/nx + vy/ny)^2) / (((vx/nx)^2)/(nx - 1) + ((vy/ny)^2)/(ny - 1))
	
	# P-values
	if(alternative == "less") {
		pval <- stats::pt(tstat, df)
	} else if(alternative == "greater") {
		pval <- stats::pt(tstat, df, lower.tail = FALSE)
	} else {
		pval <- 2 * stats::pt(-abs(tstat), df)
	}
	
	return(cbind(tstat, df, pval))
}
