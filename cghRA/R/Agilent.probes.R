# Produces a cghRA.probes object from an Agilent Feature Extraction file
# Author : Sylvain Mareschal <mareschal@ovsa.fr>

Agilent.probes = function(
	file,
	columns = c(
		rFin = "rProcessedSignal",
		gFin = "gProcessedSignal",
		flag_rIsSaturated = "rIsSaturated",
		flag_gIsSaturated = "gIsSaturated",
		flag_rIsFeatNonUnifOL = "rIsFeatNonUnifOL",
		flag_gIsFeatNonUnifOL = "gIsFeatNonUnifOL",
		flag_rIsBGNonUnifOL = "rIsBGNonUnifOL",
		flag_gIsBGNonUnifOL = "gIsBGNonUnifOL",
		flag_rIsFeatPopnOL = "rIsFeatPopnOL",
		flag_gIsFeatPopnOL = "gIsFeatPopnOL",
		flag_rIsBGPopnOL = "rIsBGPopnOL",
		flag_gIsBGPopnOL = "gIsBGPopnOL"			
	),
	...
	) {
	# File connection
	connection = file(
		description = file,
		open = "r"
	)
	
	# Subtables
	nlines = c(1, 1, -1)
	rawData = list()
	for(i in 1:length(nlines)) {
		# Type extraction
		types = scan(
			file = connection,
			what = "",
			sep = "\t",
			nlines = 1,
			quiet = TRUE,
			comment.char = "",
			quote = NULL
		)
		
		# Type conversion
		types[1] = "character"
		types[ types == "boolean" ] = "integer"
		types[ types == "float" ] = "double"
		types[ types == "text" ] = "character"
		
		# Column headers
		headers = scan(
			file = connection,
			what = "",
			sep = "\t",
			nlines = 1,
			quiet = TRUE,
			comment.char = "",
			quote = NULL
		)
		
		# Table read
		rawData[[i]] = utils::read.table(
			file = connection,
			sep = "\t",
			dec = ".",
			header = FALSE,
			nrows = nlines[i],
			col.names = headers,
			colClasses = types,
			comment.char = "",
			quote = ""
		)
		
		# Separator
		scan(
			file = connection,
			what = "",
			sep = "\t",
			nlines = 1,
			quiet = TRUE,
			comment.char = "",
			quote = NULL
		)
	}
	
	# Closing connection
	close(connection)
	
	# New object
	probes = new("cghRA.probes", name=sub("^.*/([^/]+)\\.[^\\.]+$", "\\1", file))
	
	# Probe ID
	probes$addColumn(
		rawData[[3]][,"FeatureNum"],
		"id"
	)
	
	# Columns
	for(colName in names(columns)) {
		probes$addColumn(
			rawData[[3]][, columns[colName]],
			colName
		)
	}
	
	# Log ratios
	if(all(c("gFin", "rFin") %in% probes$getColNames())) {
		probes$addColumn(
			log(probes$extract(,"rFin") / probes$extract(,"gFin"), 2),
			"logRatio"
		)
	}
	
	# Final check
	probes$check()
	
	return(probes)
}
