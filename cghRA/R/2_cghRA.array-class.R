# Reference class for cghRA array objects (merged probe values and annotation)
# Author : Sylvain Mareschal <mareschal@ovsa.fr>

# R5 class definition
setRefClass(
	Class = "cghRA.array",
	fields = list(
		design = "cghRA.design",
		probes = "cghRA.probes",
		organism = function(newValue) {
			if(missing(newValue)) { return(design$organism)
			} else                { stop("Read-only field, mutate design$organism to achieve the same purpose")
			}
		},
		assembly = function(newValue) {
			if(missing(newValue)) { return(design$assembly)
			} else                { stop("Read-only field, mutate design$assembly to achieve the same purpose")
			}
		}
	),
	contains = "crossable",
	methods = list(

as.CNA = function() {
"Returns a CNA object (DNAcopy) with the object content."
	
	return(
		suppressWarnings(
			CNA(
				genomdat = probes$extract(,"logRatio"),
				chrom = design$extract(,"chrom"),
				maploc = design$extract(,"start"),
				data.type = "logratio"
			)
		)
	)
},

as.profileCGH = function(chrom=c("merged", "levels"), quiet=FALSE) {
"Returns a profileCGH object (GLAD) with the object content.
- chrom   : single character value defining how to deal with chromosome names :
            'merged' forces chromosome arms to be merged (as chromosome arms are not handled)
            'levels' converts chromosome to integers (can be deceiving for factors)
- quiet   : single logical value, whether to warn for factor to integer conversion or not."
	
	suppressMessages(library(GLAD, quietly=TRUE))
	
	# Chromosomes
	chrom <- match.arg(chrom)
	if(chrom == "merged") {
		# Always merged
		chrom <- design$eraseArms(temp=TRUE)$extract(,"chrom")
	} else if(chrom == "levels") {
		# To integer
		chrom <- design$extract(,"chrom")
		if(!isTRUE(quiet) && is.factor(chrom)) warning("Converting factor to integer, correspondance may be broken")
		chrom <- as.integer(chrom)
	}
	
	return(
		GLAD::as.profileCGH(
			data.frame(
				LogRatio = probes$extract(,"logRatio"),
				PosOrder = 1:probes$rowCount,
				PosBase = design$extract(,"start"),
				Chromosome = chrom,
				Clone = design$extract(,"name"),
				stringsAsFactors = FALSE
			)
		)
	)

},

check = function(warn=TRUE) {
"Raises an error if the object is not valid, else returns TRUE"
	
	# Elements
	design$check(warn=warn)
	probes$check(warn=warn)
	
	# Consistency
	if(design$getRowCount() != probes$getRowCount())             stop("'design' and 'probes' row counts does not match")
	
	# Column ambiguity
	columns <- intersect(probes$colNames, design$colNames)
	for(k in columns) {
		if(!identical(design$extract(,k), probes$extract(,k)))   stop(sprintf("Column \"%s\" is ambiguous (must have identical content in 'design' and 'probes' or different names)", k))
	}
	
	return(TRUE)
},

chromosomes = function() {
"Returns the chromosome list as a vector. NULL is valid if non relevant, but should be avoided when possible."
	
	return(design$chromosomes())
},

defaultParams = function(...) {
"Returns class-specific defaults for graphical parameters. Inheriting class should overload it to define their own defaults.
- ...   : may be used by inheriting methods, especially for inter-dependant parameters."
	
	params <- callSuper(...)
	
	params$ylab <- .self$name
	params$ysub <- .self$assembly
	params$drawFun <- "draw.points"
	params$column <- "logRatio"
	params$ylim <- c(-3,3)
	params$yaxt <- "s"
	
	return(params)
},

DLRS = function(method=c("agilent", "original"), na.rm=TRUE) {
"Computes the Derivative Log Ratio Spread from the probes.
- method   : 'agilent' or 'original', implying distinct formulas."
	
	# Arg check
	method = match.arg(method)
	
	# Ordered logRatio
	logRatio = .self$probes$extract(, "logRatio")
	
	# Computation
	if (method == "agilent") {
		return(
			IQR(diff(logRatio), na.rm=na.rm) / (1.349 * sqrt(2))
		)
	} else if (method == "original") {
		return(
			sd(diff(logRatio), na.rm=na.rm) / sqrt(2)
		)
	}
},

DNAcopy = function(smooth = TRUE, ...) {
"Apply the Circular Binary Segmentation, as implemented in DNAcopy, and return a cghRA.regions object.
- smooth      : a list of arguments to be passed to smooth.CNA(), TRUE to use the default parameters or FALSE to skip smoothing.
- ...         : arguments to be passed to segment()."
	
	# CNA conversion
	cna = .self$as.CNA()
	
	# Smoothing
	if(isTRUE(smooth)) smooth = list()
	smooth$x = cna
	if(is.list(smooth)) {
		cna <- do.call(
			what = smooth.CNA,
			args = smooth
		)
	}
	
	# Segmentation
	segmentList = segment(
		x = cna,
		...
	)
	
	# DNAcopy bug fix
	segmentList$out <- subset(segmentList$out, !is.na(seg.mean) & num.mark > 0)
	
	# New object
	output <- cghRA.regions(
		chrom = factor(segmentList$out$chrom, levels=design$chromosomes()),
		strand = factor(rep(NA, nrow(segmentList$out)), levels=c("-","+")),
		start = segmentList$out$loc.start,
		end = segmentList$out$loc.end,
		probes = as.integer(segmentList$out$num.mark),
		logRatio = segmentList$out$seg.mean,
		.name = .self$name,
		.organism = .self$organism,
		.assembly = .self$assembly,
		.makeNames = TRUE
	)
	
	return(output)
},

extract = function(i=NULL, j=NULL) {
"Extracts values from 'probes' and 'design' into a data.frame.
- i   : row selection, see the R5Table method for further details.
- j   : column selection, see the R5Table method for further details."
	
	# Men at work
	if(is.expression(j)) stop("Expression indexing of columns is not handled by the cghRA.array class")
	
	# Temporary environment union for expression evaluation
	if(is.expression(i)) {
		parentEnv <- parent.env(design$values)
		parent.env(design$values) <<- probes$values
		on.exit(parent.env(design$values) <<- parentEnv, add=TRUE)
	}
	
	# Row indexing
	i = .self$design$indexes(i, "row")
	output <- cbind(
		.self$design$extract(i=i, drop=FALSE),
		.self$probes$extract(i=i, drop=FALSE)
	)
	
	# Column indexing
	if(!is.null(j)) output <- output[,j]
	
	return(output)
},

getChromEnd = function(chrom) {
"Returns as a single integer value the ending position of the object description of the given chromosome. NA (integer) is valid if non relevant, but should be avoided when possible.
- chrom   : single integer, numeric or character value, the chromosomal location. NA is not required to be handled."
	
	return(design$getChromEnd(chrom))
},

GLAD = function(chrom=c("merged", "levels"), quiet=FALSE, output=c("regions", "raw", "both"), ...) {
"Apply the Gain and Loss Analysis of Dna, as implemented in GLAD, and return a cghRA.regions object.
- chrom, quiet   : to be passed to the as.profileCGH method.
- output         : single character value defining the returned value :
                   'regions' returns a cghRA.regions object with the segmented genome
                   'raw' returns the glad() output
                   'both' adds a 'cghRA.regions' element to the glad() output list to return both
- ...            : arguments to be passed to glad()."
	
	# Checks
	output <- match.arg(output)
	
	# profileCGH conversion
	prf <- .self$as.profileCGH(chrom=chrom, quiet=quiet)
	
	# Segmentation
	results <- GLAD::glad(
		profileCGH = prf,
		...
	)
	
	# Raw output
	if(output == "raw") return(results)
	
	# cghRA.regions
	filtered <- subset(results$profileValues, OutliersTot==0)
	regions <- cghRA.regions(
		chrom = factor(tapply(X=filtered$Chromosome, INDEX=filtered$Region, FUN=unique)),
		strand = factor(rep(NA, length(unique(filtered$Region))), levels=c("-","+")),
		start = as.integer(tapply(X=filtered$PosBase, INDEX=filtered$Region, FUN=min)),
		end = as.integer(tapply(X=filtered$PosBase, INDEX=filtered$Region, FUN=max)),
		probes = as.integer(tapply(X=filtered$PosOrder, INDEX=filtered$Region, FUN=length)),
		logRatio = as.numeric(tapply(X=filtered$Smoothing, INDEX=filtered$Region, FUN=unique)),
		copies = as.integer(tapply(X=filtered$ZoneGNL, INDEX=filtered$Region, FUN=unique)),
		.name = .self$name,
		.organism = .self$organism,
		.assembly = .self$assembly,
		.makeNames = TRUE
	)
	
	# Regions output
	if(output == "regions") return(regions)
	
	# Both output
	results$cghRA.regions <- regions
	return(results)
},

initialize = function(design=new("cghRA.design"), probes=new("cghRA.probes"), ...) { # "organism" and "assembly" are wrappers, do not initialize !
	callSuper(...)
	initFields(design=design, probes=probes)
},

MAplot = function(pch="+", cex=0.6, xlab="A (log2 mean)", ylab="M (log2 ratio)", ...) {
"MA plot of all the probes.
- ...   : arguments to be passed to plot()."
	
	# Data
	R <- log(.self$probes$extract(,"rFin"), 2)
	G <- log(.self$probes$extract(,"gFin"), 2)
	M <- R - G
	A <- (R + G) / 2
	
	# Plot
	plot(
		x = A,
		y = M,
		xlab = xlab,
		ylab = ylab,
		cex = cex,
		pch = pch,
		...
	)
	abline(h=0, col="#666666")
},

maskByFlag = function(flags="^flag_", pattern=TRUE, multiple=c("any", "all"), na=FALSE) {
"Replaces logRatios of flagged probes by NA.
- flags      : character vector, the columns to coerce as boolean and use as flags.
- pattern    : single logical value, whether to consider 'flags' as regular expressions or fixed values.
- multiple   : mask a probe when 'all' its flag columns are TRUE or when 'any' is."
	
	# Arg check
	multiple = match.arg(multiple)
	
	# Flag collection
	flagColumns = character(0)
	if(isTRUE(pattern)) {
		# Multiple patterns
		for(flag in flags) {
			flagColumns = c(
				flagColumns,
				grep(
					pattern = flag,
					x = .self$probes$colNames,
					value = TRUE
				)
			)
		}
	} else {
		# Fixed column names
		flagColumns = flags
	}
	
	# Flag check
	if(any(!flagColumns %in% .self$probes$colNames)) stop("All 'flags' must correspond to existing 'probes' columns.")
	
	# Boolean vector
	bool = as.matrix(.self$probes$extract(, flagColumns))
	storage.mode(bool) = "logical"
	bool[ is.na(bool) ] = na
	bool = apply(bool, 1, multiple)
	
	# Modification
	.self$probes$fill(bool, "logRatio", as.numeric(NA))
},

replicates = function(fun=median, na.rm=TRUE, ...) {
"Apply 'fun' to replicated probes (same name), masking all members but one.
- fun   : single character value, the function to apply.
- ...   : to be passed to 'fun'."
	
	# Looking for duplicated probes
	prbName <- .self$design$extract(,"name")
	dup <- duplicated(prbName, fromLast=FALSE) | duplicated(prbName, fromLast=TRUE)
	
	if(any(dup)) {
		# All probes
		prbLogRatio <- .self$probes$extract(,"logRatio")
		prbIndex <- 1:.self$probes$getRowCount()
		
		# Focus on replicated probes
		prbName <- prbName[dup]
		prbLogRatio <- prbLogRatio[dup]
		prbIndex <- prbIndex[dup]
		
		# Masking all duplicates
		.self$probes$fill(prbIndex, "logRatio", as.numeric(NA))
		
		# Filling first indexes only
		indexes <- sapply(tapply(X=prbIndex, INDEX=prbName, FUN=c), "[", 1)
		values <- tapply(X=prbLogRatio, INDEX=prbName, FUN=fun, na.rm=na.rm, ...)
		.self$probes$fill(indexes, "logRatio", as.numeric(values))
	}
},

show = function(include=FALSE, fieldWidth=10) {
"Interactive printing
- include   : single logical value, if TRUE class name will not be printed."
	
	# Class name
	if(!isTRUE(include)) {
		cat("\n  \"cghRA.array\" reference class object\n")
	} else {
		cat("\n  Extends \"cghRA.array\"\n")
	}
	
	# Fields
	cat(sprintf("  %-*s : cghRA.design <%d x %d>\n", fieldWidth, "design", design$rowCount, design$colCount))
	cat(sprintf("  %-*s : cghRA.probes <%d x %d>\n", fieldWidth, "probes", probes$rowCount, probes$colCount))
},

slice = function(chrom, start, end) {
"Extracts elements in the specified window as a data.frame.
- chrom   : single integer or character value, the chromosomal location.
- start   : single integer value, inferior boundary of the window.
- end     : single integer value, superior boundary of the window."
	
	if(is.numeric(chrom)) chrom <- as.integer(chrom)
	if(is.numeric(start)) start <- as.integer(start)
	if(is.numeric(end))   end <- as.integer(end)
	
	# Collect columns
	columns <- paste(
		c(
			sprintf("%s=design$values[[\"%s\"]]", design$colNames, design$colReferences),
			sprintf("%s=probes$values[[\"%s\"]]", probes$colNames, probes$colReferences)
		),
		collapse = ", "
	)
	
	# Create call
	if(is.character(chrom)) {
		k <- parse(text=sprintf(".External(\"track\", PACKAGE=\"Rgb\", mode='sub', \"%s\", %d, %d, design$index, %s)", chrom, start, end, columns))[[1]]
	} else if(is.integer(chrom)) {
		k <- parse(text=sprintf(".External(\"track\", PACKAGE=\"Rgb\", mode='sub', %d, %d, %d, design$index, %s)", chrom, start, end, columns))[[1]]
	} else {
		stop("'chrom' must be integer or character")
	}
	
	# Evaluation
	return(eval(k))
},

spatial = function(filename=sprintf("%s.png", .self$name), palSize=254, palEnds=c("#0000FF", "#FFFF00"), ...) {
"Produces a spatial representation of the logRatios, to identify spatial biases.
- filename   : single character value, the path to the PNG output.
- palSize    : single integer value, the amount of color levels for logRatios. Should be lesser or equal to 254 to produce small PNG files.
- palEnds    : character vector to be passed to colorRampPalette() for palette generation."
	
	# Check
	if(any(! c("row","col") %in% .self$design$colNames)) stop("Unable to draw spatial plot without coordinates")
	
	# Color palette
	palSize = as.integer(palSize)
	pal = colorRampPalette(palEnds)(palSize)
	
	# Data extraction
	logRatios = .self$probes$extract(,"logRatio")
	cols = .self$design$extract(,"col")
	rows = .self$design$extract(,"row")
	
	# LogRatio matrix
	intensities = matrix(
		data = as.double(NA),
		ncol = max(cols),
		nrow = max(rows)
	)
	intensities[ cbind(rows, cols) ] = logRatios
	
	# Density limits
	dns = density(logRatios, na.rm=TRUE, n=1024, kernel="gaussian")
	xmin = dns$x[ head(which(dns$y > max(dns$y) / 200), 1) ]
	xmax = dns$x[ tail(which(dns$y > max(dns$y) / 200), 1) ]
	
	# Intensities to [0:255]
	intensities[ intensities < xmin ] = xmin
	intensities[ intensities > xmax ] = xmax
	intensities = (intensities - xmin) / (xmax - xmin)
	intensities = round(intensities*(palSize-1))
	
	# To RGB colors
	raster = matrix(data=pal[ as.integer(intensities)+1L ], nrow=nrow(intensities), ncol=ncol(intensities))
	raster[ is.na(raster) ] = "#000000"
	
	# Graphical parameters
	res = 72
	inchesWidth = 0.2 + ncol(intensities)/res + 0.2
	inchesHeight = (5*0.393700787) + nrow(intensities)/res + 0.2
	png(filename=filename, width=inchesWidth, height=inchesHeight, units="in", res=res)
	layout(matrix(1:2, ncol=1), heights=c(lcm(5), 1))
	
	# Density plot
	par(mai=c(0.2, 0.2, 0.5, 0.2))
	breaks = seq(from=xmin, to=xmax, length.out=(palSize+1))
	breaks[1] = -100
	breaks[palSize+1] = 100
	plot(x=NA, y=NA, xlim=range(dns$x), ylim=range(dns$y), xlab="", ylab="", yaxt="n", cex.axis=0.75, mgp=c(0,0,0), tcl=-0.2)
	rect(
		xleft = breaks[1:palSize],
		xright = breaks[2:(palSize+1)],
		ybottom = -10,
		ytop = 10,
		border = pal,
		col = pal
	)
	par(new=TRUE)
	plot(dns, xlab="", ylab="", xaxt="n", yaxt="n", main=sprintf("%s logRatio", .self$name), zero.line = FALSE, cex.main=1, font.main=1)
	
	# Intensities plot
	par(mai=c(0.2, 0.2, 0, 0.2))
	plot(x=NA, y=NA, xlim=c(0, ncol(intensities)), ylim=c(0, nrow(intensities)), xlab="", ylab="", xaxs="i", yaxs="i", xaxt="n", yaxt="n", bty="n")
	rasterImage(
		image = as.raster(raster),
		xleft = 0,
		ybottom = 0,
		xright = ncol(raster),
		ytop = nrow(raster),
		angle = 0,
		interpolate = FALSE
	)
	
	dev.off()
},

WACA = function() {
"Apply the Waves aCGH Correction Algorithm (Lepretre et al. 2009) to the array logRatios."
	
	.self$probes$fill(
		,
		"logRatio",
		cghRA::WACA(
			probeNames = .self$design$extract(,"name"),
			probeLogRatios = .self$probes$extract(,"logRatio"),
			bias = .self$design$extract(, c("wGC150", "wGC500", "wGCprobe", "wGCfrag", "wFragSize")),
			forceBiasOrdering = FALSE
		)
	)
}

	)
)

# Constructor
cghRA.array <- function(.design, .probes, .name, .parameters, warn=TRUE) {
	# New object
	object <- new("cghRA.array")
	if(!missing(.design))                   object$design <- .design
	if(!missing(.probes))                   object$probes <- .probes
	if(!missing(.parameters))               object$parameters <- .parameters
	
	# Name field ('drawable' field actually)
	if(!missing(.name))          { object$name <- .name
	} else if(!missing(.probes)) { object$name <- .probes$name
	}
	
	# Probe order consistency
	if(!identical(object$design$extract(,"id"), object$probes$extract(,"id"))) {
		if(max(object$design$extract(,"id")) != max(object$probes$extract(,"id"))) warning("Max 'id' does not match, probes and design may not be compatible")
		f <- factor(x=object$design$extract(,"id"), levels=object$probes$extract(,"id"))
		object$probes$rowOrder(as.integer(f))
		object$probes$fill(,"id", object$design$extract(,"id"))
	}
	
	# Check
	object$check(warn=warn)
	
	return(object)
}

