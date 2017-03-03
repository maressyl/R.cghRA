# Interactive Tcl-Tk track annotation and CNV-scoring
# Author : Sylvain Mareschal <mareschal@ovsa.fr>

tk.annotate = function(
		globalTopLevel,
		localTopLevel
		)
	{
	
	## FUNCTIONS ##

	trackFilesBrowse <- function() {
		trackFiles <<- tk.file(
			title = "Choose Rgb-compliant track files to process",
			typeNames = c("cghRA modelized region files", "cghRA raw region files", "Rgb track file"),
			typeExt = c(".copies.rdt", ".regions.rdt", ".rdt"),
			multiple = TRUE,
			mandatory = FALSE,
			type = "open",
			parent = localTopLevel
		)
		tcltk::tclvalue(trackFilesValue) <- sprintf("%i file(s) selected", length(trackFiles))
	}
	
	annoTrackBrowse <- function() {
		tcltk::tclvalue(annoTrackValue) <- tk.file(
			title = "Choose a track to use for annotation",
			typeNames = "Rgb track file",
			typeExt = ".rdt",
			multiple = FALSE,
			mandatory = FALSE,
			type = "open",
			parent = localTopLevel
		)
	}
	
	designBrowse <- function() {
		tcltk::tclvalue(designValue) <- tk.file(
			title = "Choose a processed design file",
			typeNames = "R data file",
			typeExt = ".rdt",
			multiple = FALSE,
			mandatory = FALSE,
			type = "open",
			parent = localTopLevel
		)
	}
	
	dgvTrackBrowse <- function() {
		tcltk::tclvalue(dgvTrackValue) <- tk.file(
			title = "Choose a track listing polymorphisms",
			typeNames = "Rgb track file",
			typeExt = ".rdt",
			multiple = FALSE,
			mandatory = FALSE,
			type = "open",
			parent = localTopLevel
		)
	}
	
	inputError <- function(message, icon="error", panel) {
		# Status
		if(panel == "cnv") tcltk::tclvalue(cnvStatusValue) <- "Error"
		
		# Dialog box
		tcltk::tkmessageBox(
			parent = localTopLevel,
			icon = icon,
			type = "ok",
			title = "Input error",
			message = message
		)
		
	}
	
	annoProcess <- function() {
		# Data collection
		annoTrackFile <- tcltk::tclvalue(annoTrackValue)
		type <- tcltk::tclvalue(annoTypeValue)
		fuzziness <- as.integer(tcltk::tclvalue(fuzzinessValue))
		maxElements <- as.integer(tcltk::tclvalue(maxElementsValue))
		
		# Data check
		if(!file.exists(annoTrackFile)) return(inputError(message=sprintf("Annotation track file can not be found"), panel="anno"))
		if(length(trackFiles) > 0) {
			for(f in trackFiles) {
				if(!file.exists(f)) return(inputError(message=sprintf("A region file can not be found :\n%s", f), panel="anno"))
			}
		} else {
			return(inputError(message="No region file provided", panel="anno"))
		}
		
		# Processing
		tcltk::tkconfigure(localTopLevel, cursor="watch")
		handle(
			expr = {
				# Parse annotation
				annotation <- Rgb::readRDT(annoTrackFile)
				colname <- paste(annotation$name, type, sep=".")
				
				# Set progression bar
				tcltk::tkconfigure(annoBar, maximum=length(trackFiles))
				tcltk::tclvalue(annoProgression) <- 0
				
				# Process region files
				for(trackFile in trackFiles) {
					# Read file
					track <- Rgb::readRDT(trackFile)
					
					# Process
					values <- track$cross(
						annotation = annotation,
						type = type,
						fuzziness = fuzziness,
						maxElements = maxElements
					)
					
					# Replace or create column
					if(colname %in% track$getColNames()) { track$fill(NULL, colname, values)
					} else                               { track$addColumn(values, colname)
					}
					
					# Save back
					Rgb::saveRDT(track, trackFile)
					
					# Increment bar
					tcltk::tclvalue(annoProgression) <- as.numeric(tcltk::tclvalue(annoProgression)) + 1
					tcltk::tcl("update", "idletasks")
				}
			},
			messageHandler = NULL,
			warningHandler = function(w) {
				tcltk::tkmessageBox(parent=localTopLevel, icon="warning", type="ok", title="R warning", message=conditionMessage(w))
			},
			errorHandler = function(e) {
				tcltk::tkconfigure(localTopLevel, cursor="arrow")
				return(tcltk::tkmessageBox(parent=localTopLevel, icon="error", type="ok", title="R error", message=conditionMessage(e)))
			}
		)
		
		# Switch back the cursor
		tcltk::tkconfigure(localTopLevel, cursor="arrow")
		
		# Done
		tcltk::tkmessageBox(
			parent = localTopLevel,
			icon = "info",
			type = "ok",
			title = "Annotation complete",
			message = sprintf("Column \"%s\" was added or updated in %i file(s).", colname, length(trackFiles))
		)
	}
	
	cnvProcess <- function() {
		# Data collection
		designFile <- tcltk::tclvalue(designValue)
		dgvTrackFile <- tcltk::tclvalue(dgvTrackValue)
		filter <- as.logical(as.integer(tcltk::tclvalue(filtApplyValue)))
		fillGaps <- as.logical(as.integer(tcltk::tclvalue(filtFillValue)))
		threshold <- as.numeric(tcltk::tclvalue(filtThresholdValue))
		
		# Data check
		if(!file.exists(designFile))   return(inputError(message=sprintf("Design file can not be found"), panel="cnv"))
		if(!file.exists(dgvTrackFile)) return(inputError(message=sprintf("Polymorphism track file can not be found"), panel="cnv"))
		if(length(trackFiles) > 0) {
			for(f in trackFiles) {
				if(!file.exists(f)) return(inputError(message=sprintf("A region file can not be found :\n%s", f), panel="cnv"))
			}
		} else {
			return(inputError(message="No region file provided", panel="cnv"))
		}
		
		# Status
		tcltk::tclvalue(cnvStatusValue) <- "Launched"
		tcltk::tkconfigure(localTopLevel, cursor="watch")
		tcltk::tclvalue(cnvProgression) <- 0
		tcltk::tcl("update", "idletasks")
		
		# Processing
		handle(
			expr = {
				# Parse design
				design <- Rgb::readRDT(designFile)
				
				# Map file
				mapFile <- sprintf("%s_map.rds", sub("\\.rdt$", "", dgvTrackFile))
				if(file.exists(mapFile)) {
					# Existing one, check compatibility
					dgvMap <- readRDS(mapFile)
					if(!inherits(dgvMap, "segmentMap"))      stop(sprintf("'%s' file exists but do not contain a 'segmentMap' object", basename(mapFile)))
					if(dgvMap$designName != design$name)     stop(sprintf("The map stored in '%s' was built with another design file (name mistmatch)", basename(mapFile)))
					if(dgvMap$designSize != design$rowCount) stop(sprintf("The map stored in '%s' was built with another design file (size mistmatch)", basename(mapFile)))
				} else {
					# Status
					tcltk::tclvalue(cnvStatusValue) <- "Bank parsing"
					tcltk::tcl("update", "idletasks")
					
					# Parse the polymorphism track
					dgvTrack <- Rgb::readRDT(dgvTrackFile)
					
					# Status
					tcltk::tclvalue(cnvStatusValue) <- "Bank mapping"
					tcltk::tcl("update", "idletasks")
					
					# Set progression bar
					tcltk::tkconfigure(cnvBar, maximum=1L)
					
					# Map the polymorphism track
					dgvMap <- map2design(dgvTrack, design)
					
					# Save the map
					saveRDS(dgvMap, file=mapFile)
				}
				
				# Status
				tcltk::tclvalue(cnvStatusValue) <- "Sample mapping"
				tcltk::tcl("update", "idletasks")
				
				# Set progression bar
				tcltk::tkconfigure(cnvBar, maximum=length(trackFiles))
				tcltk::tclvalue(cnvProgression) <- 0
				
				# Process region files
				for(trackFile in trackFiles) {
					# Read file
					track <- Rgb::readRDT(trackFile)
					
					# Remap segments to current design
					trackMap <- map2design(track, design)
	
					# Compute CNV score for each segment
					score <- cnvScore(trackMap, dgvMap, expand=TRUE)
					if("cnvScore" %in% track$getColNames()) { track$fill(NULL, "cnvScore", score)
					} else                                  { track$addColumn(score, "cnvScore")
					}
					
					# Filter out polymorphisms
					if(isTRUE(filter)) {
						track$rowOrder(which(track$extract(,"cnvScore") < threshold))
						if(isTRUE(fillGaps)) track$fillGaps()
					}
					
					# Save back
					Rgb::saveRDT(track, trackFile)
				}
				
				# Status
				tcltk::tclvalue(cnvStatusValue) <- "Done"
				tcltk::tkconfigure(localTopLevel, cursor="arrow")
				tcltk::tcl("update", "idletasks")
		
				# Done
				tcltk::tkmessageBox(
					parent = localTopLevel,
					icon = "info",
					type = "ok",
					title = "CNV score computation complete",
					message = sprintf("Column \"cnvScore\" was added or updated in %i file(s).", length(trackFiles))
				)
			},
			messageHandler = function(m) {
				text <- conditionMessage(m)
				if(grepl("^([0-9]+)/([0-9]+)", text)) {
					step <- floor(as.numeric(tcltk::tclvalue(cnvProgression)))
					current <- as.integer(sub("^([0-9]+)/([0-9]+)", "\\1", text))
					total <- as.integer(sub("^([0-9]+)/([0-9]+)", "\\2", text))
					tcltk::tclvalue(cnvProgression) <- min(step + current / total, total)
					tcltk::tcl("update", "idletasks")
				}
			},
			warningHandler = function(w) {
				tcltk::tkmessageBox(parent=localTopLevel, icon="warning", type="ok", title="R warning", message=conditionMessage(w))
			},
			errorHandler = function(e) {
				tcltk::tkconfigure(localTopLevel, cursor="arrow")
				tcltk::tclvalue(cnvStatusValue) <- "Error"
				return(tcltk::tkmessageBox(parent=localTopLevel, icon="error", type="ok", title="R error", message=conditionMessage(e)))
			}
		)
	}
	
	exportProcess <- function() {
		# Data collection
		sep <- tcltk::tclvalue(exportSepValue)
		dec <- tcltk::tclvalue(exportDecValue)
		quote <- as.logical(as.integer(tcltk::tclvalue(exportQuoteValue)))
		header <- as.logical(as.integer(tcltk::tclvalue(exportHeaderValue)))
		
		# Data check
		if(length(trackFiles) > 0) {
			for(f in trackFiles) {
				if(!file.exists(f)) return(inputError(message=sprintf("A region file can not be found :\n%s", f), panel="cnv"))
			}
		} else {
			return(inputError(message="No region file provided", panel="cnv"))
		}
		
		# Status
		tcltk::tkconfigure(localTopLevel, cursor="watch")
		tcltk::tkconfigure(exportBar, maximum=length(trackFiles))
		tcltk::tclvalue(exportProgression) <- 0
		tcltk::tcl("update", "idletasks")
		
		# Processing
		handle(
			expr = {
				# Process region files
				for(trackFile in trackFiles) {
					# Read file
					track <- Rgb::readRDT(trackFile)
					
					# Extract content
					content <- track$extract()
					
					# Exported file name
					outFile <- sub("\\.rdt$", ".csv", trackFile, ignore.case=TRUE)
					
					# Header
					if(isTRUE(header)) {
						# Collect fields
						metaData <- character(0)
						fields <- setdiff(names(track$getRefClass()$fields()), names(getRefClass("refTable")$fields()))
						for(fieldName in fields) {
							fieldContent <- track$field(fieldName)
							if(is.atomic(fieldContent) && length(fieldContent) == 1) {
								metaData[ fieldName ] <- as.character(fieldContent)
							}
						}
					
						# Write header
						cat(sprintf("# %s = %s\n", names(metaData), metaData), file=outFile, append=FALSE, sep="")
					}
					
					# Content
					suppressWarnings(utils::write.table(content, file=outFile, sep=sep, dec=dec, quote=quote, row.names=FALSE, col.names=TRUE, append=TRUE))
					
					# Increment bar
					tcltk::tclvalue(exportProgression) <- as.numeric(tcltk::tclvalue(exportProgression)) + 1
					tcltk::tcl("update", "idletasks")
				}
				
				# Switch back the cursor
				tcltk::tkconfigure(localTopLevel, cursor="arrow")
				
				# Done
				tcltk::tkmessageBox(
					parent = localTopLevel,
					icon = "info",
					type = "ok",
					title = "Export complete",
					message = sprintf("%i RDT files exported to CSV", length(trackFiles))
				)
			},
			messageHandler = NULL,
			warningHandler = function(w) {
				tcltk::tkmessageBox(parent=localTopLevel, icon="warning", type="ok", title="R warning", message=conditionMessage(w))
			},
			errorHandler = function(e) {
				tcltk::tkconfigure(localTopLevel, cursor="arrow")
				tcltk::tclvalue(cnvStatusValue) <- "Error"
				return(tcltk::tkmessageBox(parent=localTopLevel, icon="error", type="ok", title="R error", message=conditionMessage(e)))
			}
		)
	}
	
	
	
	## INTERFACE ##

	embedded <- !missing(globalTopLevel) && !missing(localTopLevel)
	if(!embedded) {
		# Linux default style
		if(.Platform$OS.type == "unix") try(tcltk::tcl("ttk::style", "theme", "use", "clam"), silent=TRUE)
		
		# Top level
		globalTopLevel <- localTopLevel <- tcltk::tktoplevel(class="cghRA")
		tcltk::tktitle(localTopLevel) <- "cghRA - Annotate regions"
		icon16 <- tcltk::tcl("image", "create", "photo", file=system.file("cghRA_16x16.gif", package="cghRA"))
		icon32 <- tcltk::tcl("image", "create", "photo", file=system.file("cghRA_32x32.gif", package="cghRA"))
		tcltk::tcl("wm", "iconphoto", localTopLevel, "-default", icon16, icon32)
	}
	
	# Evolutive width
	tcltk::tkgrid.columnconfigure(localTopLevel, 1, weight=1)
	
		# Input frame
		inputFrame <- tcltk::ttklabelframe(parent=localTopLevel, relief="groove", borderwidth=2, text="Common parameters")
		tcltk::tkgrid(inputFrame, column=1, row=1, padx=5, pady=5, sticky="we")
		tcltk::tkgrid.columnconfigure(inputFrame, 2, weight=1)
			
			# Region files
			trackFiles <- character(0)
			trackFilesValue <- tcltk::tclVar("0 file(s) selected")
			trackFilesButton <- tcltk::tkbutton(parent=inputFrame, text="Track(s) to process", command=trackFilesBrowse, width=20)
			trackFilesLabel <- tcltk::tklabel(parent=inputFrame, textvariable=trackFilesValue)
			tcltk::tkgrid(trackFilesButton, column=1, row=2, padx=5, pady=5)
			tcltk::tkgrid(trackFilesLabel, column=2, row=2, padx=5, pady=5, sticky="w")
			
		# Annotate frame
		annotateFrame <- tcltk::ttklabelframe(parent=localTopLevel, relief="groove", borderwidth=2, text="Add an annotation column")
		tcltk::tkgrid(annotateFrame, column=1, row=2, padx=5, pady=5, sticky="we")
		tcltk::tkgrid.columnconfigure(annotateFrame, 2, weight=1)
		
			# Annotation track
			annoTrackValue <- tcltk::tclVar("")
			annoTrackButton <- tcltk::tkbutton(parent=annotateFrame, command=annoTrackBrowse, text="Annotation track", width=20)
			annoTrackEntry <- tcltk::tkentry(parent=annotateFrame, textvariable=annoTrackValue, width=80)
			tcltk::tkgrid(annoTrackButton, column=1, row=1, padx=5, pady=5)
			tcltk::tkgrid(annoTrackEntry, column=2, row=1, padx=5, pady=5, sticky="we")
			
			# Annotation parameters
			annoParamLabel <- tcltk::tklabel(parent=annotateFrame, text="Parameters", width=20)
			annoParamFrame <- tcltk::tkframe(parent=annotateFrame)
			tcltk::tkgrid(annoParamLabel, column=1, row=2, padx=5, pady=5)
			tcltk::tkgrid(annoParamFrame, column=2, row=2, sticky="we")
			tcltk::tkgrid.columnconfigure(annoParamFrame, 7, weight=1)
			
				# Type combobox
				annoTypeValue <- tcltk::tclVar("name")
				annoTypeLabel <- tcltk::tklabel(parent=annoParamFrame, text="Crossing type :")
				annoTypeCombo <- tcltk::ttkcombobox(parent=annoParamFrame, values=c("name", "cover", "count", "cytoband"), textvariable=annoTypeValue, width=10, justify="center")
				tcltk::tkgrid(annoTypeLabel, column=1, row=1, padx=5, pady=5)
				tcltk::tkgrid(annoTypeCombo, column=2, row=1, padx=c(5,20), pady=5, sticky="w")
				
				# Fuzziness input
				fuzzinessValue <- tcltk::tclVar("50000")
				fuzzinessLabel <- tcltk::tklabel(parent=annoParamFrame, text="Fuzziness (bp) :")
				fuzzinessEntry <- tcltk::tkentry(parent=annoParamFrame, textvariable=fuzzinessValue, width=10, justify="center")
				tcltk::tkgrid(fuzzinessLabel, column=3, row=1, padx=5, pady=5)
				tcltk::tkgrid(fuzzinessEntry, column=4, row=1, padx=c(5,20), pady=5, sticky="w")
				
				# Max elements input
				maxElementsValue <- tcltk::tclVar("30")
				maxElementsLabel <- tcltk::tklabel(parent=annoParamFrame, text="Max elements :")
				maxElementsEntry <- tcltk::tkentry(parent=annoParamFrame, textvariable=maxElementsValue, width=5, justify="center")
				tcltk::tkgrid(maxElementsLabel, column=5, row=1, padx=5, pady=5)
				tcltk::tkgrid(maxElementsEntry, column=6, row=1, padx=c(5,20), pady=5, sticky="w")
			
			# Annotate button
			annoButton <- tcltk::tkbutton(parent=annotateFrame, text="Compute annotation", background="#A9BEE1", font="Verdana 8 bold", command=annoProcess, width=20)
			tcltk::tkgrid(annoButton, column=1, row=3, padx=5, pady=5)
			
			# Annotate progression bar
			annoProgression <- tcltk::tclVar("0")
			annoBar <- tcltk::ttkprogressbar(parent=annotateFrame, variable=annoProgression, maximum=1)
			tcltk::tkgrid(annoBar, column=2, row=3, padx=5, pady=5, sticky="ew")
			
		# cnvScore frame
		cnvFrame <- tcltk::ttklabelframe(parent=localTopLevel, relief="groove", borderwidth=2, text="Polymorphism likelihood scoring")
		tcltk::tkgrid(cnvFrame, column=1, row=3, padx=5, pady=5, sticky="we")
		tcltk::tkgrid.columnconfigure(cnvFrame, 3, weight=1)
			
			# Design file
			designValue <- tcltk::tclVar("")
			designButton <- tcltk::tkbutton(parent=cnvFrame, text="Design file", command=designBrowse, width=20)
			designEntry <- tcltk::tkentry(parent=cnvFrame, textvariable=designValue)
			tcltk::tkgrid(designButton, column=1, row=1, padx=5, pady=5)
			tcltk::tkgrid(designEntry, column=2, columnspan=2, row=1, padx=5, pady=5, sticky="we")
			
			# Polymorphism track
			dgvTrackValue <- tcltk::tclVar("")
			dgvTrackButton <- tcltk::tkbutton(parent=cnvFrame, command=dgvTrackBrowse, text="Polymorphism track", width=20)
			dgvTrackEntry <- tcltk::tkentry(parent=cnvFrame, textvariable=dgvTrackValue, width=80)
			tcltk::tkgrid(dgvTrackButton, column=1, row=2, padx=5, pady=5)
			tcltk::tkgrid(dgvTrackEntry, column=2, columnspan=2, row=2, padx=5, pady=5, sticky="we")
			
			# Filtering parameters
			filtParamLabel <- tcltk::tklabel(parent=cnvFrame, text="Filtering", width=20)
			filtParamFrame <- tcltk::tkframe(parent=cnvFrame)
			tcltk::tkgrid(filtParamLabel, column=1, row=3, padx=5, pady=5)
			tcltk::tkgrid(filtParamFrame, column=2, columnspan=2, row=3, sticky="we")
			tcltk::tkgrid.columnconfigure(filtParamFrame, 7, weight=1)
				
				# Filtering checkbox
				filtApplyValue <- tcltk::tclVar("1")
				filtApplyCheck <- tcltk::tkcheckbutton(parent=filtParamFrame, variable=filtApplyValue)
				filtApplyLabel <- tcltk::tklabel(parent=filtParamFrame, text="Filter :")
				tcltk::tkgrid(filtApplyLabel, column=1, row=1, padx=5, pady=5)
				tcltk::tkgrid(filtApplyCheck, column=2, row=1, padx=c(5,20), pady=5, sticky="w")
				
				# Threshold input
				filtThresholdValue <- tcltk::tclVar("1")
				filtThresholdLabel <- tcltk::tklabel(parent=filtParamFrame, text="Max score :")
				filtThresholdEntry <- tcltk::tkentry(parent=filtParamFrame, textvariable=filtThresholdValue, width=5, justify="center")
				tcltk::tkgrid(filtThresholdLabel, column=3, row=1, padx=5, pady=5)
				tcltk::tkgrid(filtThresholdEntry, column=4, row=1, padx=c(5,20), pady=5, sticky="w")
				
				# Fill gap checkbox
				filtFillValue <- tcltk::tclVar("1")
				filtFillCheck <- tcltk::tkcheckbutton(parent=filtParamFrame, variable=filtFillValue)
				filtFillLabel <- tcltk::tklabel(parent=filtParamFrame, text="Fill gaps :")
				tcltk::tkgrid(filtFillLabel, column=5, row=1, padx=5, pady=5)
				tcltk::tkgrid(filtFillCheck, column=6, row=1, padx=c(5,20), pady=5, sticky="w")
			
			# cnvScore button
			cnvButton <- tcltk::tkbutton(parent=cnvFrame, text="Compute CNV score", background="#A9BEE1", font="Verdana 8 bold", command=cnvProcess, width=20)
			tcltk::tkgrid(cnvButton, column=1, row=4, padx=5, pady=5)
			
			# cnvScore progression bar
			cnvStatusValue <- tcltk::tclVar("Ready")
			cnvProgression <- tcltk::tclVar("0")
			cnvStatusLabel <- tcltk::tklabel(parent=cnvFrame, textvariable=cnvStatusValue, width=15)
			cnvBar <- tcltk::ttkprogressbar(parent=cnvFrame, variable=cnvProgression, maximum=1)
			tcltk::tkgrid(cnvStatusLabel, column=2, row=4, padx=5, pady=5, sticky="ew")
			tcltk::tkgrid(cnvBar, column=3, row=4, padx=5, pady=5, sticky="ew")
		
		# Export frame
		exportFrame <- tcltk::ttklabelframe(parent=localTopLevel, relief="groove", borderwidth=2, text="Export RDT files to a spread-sheet editor")
		tcltk::tkgrid(exportFrame, column=1, row=4, padx=5, pady=5, sticky="we")
		tcltk::tkgrid.columnconfigure(exportFrame, 2, weight=1)
		
			# Export parameters
			exportParamLabel <- tcltk::tklabel(parent=exportFrame, text="File format", width=20)
			exportParamFrame <- tcltk::tkframe(parent=exportFrame)
			tcltk::tkgrid(exportParamLabel, column=1, row=1, padx=5, pady=5)
			tcltk::tkgrid(exportParamFrame, column=2, row=1, sticky="we")
			tcltk::tkgrid.columnconfigure(exportParamFrame, 8, weight=1)
			
				# Column separator combobox
				exportSepValue <- tcltk::tclVar(",")
				exportSepLabel <- tcltk::tklabel(parent=exportParamFrame, text="Column separator :")
				exportSepCombo <- tcltk::ttkcombobox(parent=exportParamFrame, values=c(",", ";", "\\t"), textvariable=exportSepValue, width=2, justify="center")
				tcltk::tkgrid(exportSepLabel, column=1, row=1, padx=5, pady=5)
				tcltk::tkgrid(exportSepCombo, column=2, row=1, padx=c(5,20), pady=5, sticky="w")
			
				# Decimal separator combobox
				exportDecValue <- tcltk::tclVar(".")
				exportDecLabel <- tcltk::tklabel(parent=exportParamFrame, text="Decimal separator :")
				exportDecCombo <- tcltk::ttkcombobox(parent=exportParamFrame, values=c(".", ","), textvariable=exportDecValue, width=2, justify="center")
				tcltk::tkgrid(exportDecLabel, column=3, row=1, padx=5, pady=5)
				tcltk::tkgrid(exportDecCombo, column=4, row=1, padx=c(5,20), pady=5, sticky="w")
				
				# Filtering checkbox
				exportQuoteValue <- tcltk::tclVar("1")
				exportQuoteCheck <- tcltk::tkcheckbutton(parent=exportParamFrame, variable=exportQuoteValue)
				exportQuoteLabel <- tcltk::tklabel(parent=exportParamFrame, text="Cell quoting :")
				tcltk::tkgrid(exportQuoteLabel, column=5, row=1, padx=5, pady=5)
				tcltk::tkgrid(exportQuoteCheck, column=6, row=1, padx=c(5,20), pady=5, sticky="w")
			
				# Header checkbox
				exportHeaderValue <- tcltk::tclVar("1")
				exportHeaderCheck <- tcltk::tkcheckbutton(parent=exportParamFrame, variable=exportHeaderValue)
				exportHeaderLabel <- tcltk::tklabel(parent=exportParamFrame, text="Header :")
				tcltk::tkgrid(exportHeaderLabel, column=7, row=1, padx=5, pady=5)
				tcltk::tkgrid(exportHeaderCheck, column=8, row=1, padx=c(5,20), pady=5, sticky="w")
			
			# Export button
			exportButton <- tcltk::tkbutton(parent=exportFrame, text="Export to CSV", background="#A9BEE1", font="Verdana 8 bold", command=exportProcess, width=20)
			tcltk::tkgrid(exportButton, column=1, row=2, padx=5, pady=5)
			
			# Annotate progression bar
			exportProgression <- tcltk::tclVar("0")
			exportBar <- tcltk::ttkprogressbar(parent=exportFrame, variable=exportProgression, maximum=1)
			tcltk::tkgrid(exportBar, column=2, row=2, padx=5, pady=5, sticky="ew")
		
}

