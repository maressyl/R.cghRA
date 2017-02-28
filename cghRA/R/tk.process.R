# Interactive Tcl-Tk array processing
# Author : Sylvain Mareschal <mareschal@ovsa.fr>

tk.process = function(
		globalTopLevel,
		localTopLevel		
		)
	{
	
	## FUNCTIONS ##
	
	probeFilesBrowse <- function() {
		probeFiles <<- tk.file(
			title = "Choose Agilent raw data files",
			typeNames = c("Feature Extraction output", "Processed probe files", "Any file type"),
			typeExt = c(".txt .txt.gz .gz", ".probes.rdt", "*"),
			multiple = TRUE,
			mandatory = FALSE,
			type = "open",
			parent = localTopLevel
		)
		tcltk::tclvalue(probeFilesValue) <- sprintf("%i file(s) selected", length(probeFiles))
	}
	
	outputBrowse <- function() {
		tcltk::tclvalue(outputValue) <- tk.folder(
			title = "Where to save processed files",
			mustexist = TRUE,
			mandatory = FALSE
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
	
	maskAction <- function() {
		if(tcltk::tclvalue(maskCheckValue) == "1") {
			# Enable
			tcltk::tkconfigure(maskEntry, state="normal")
			tcltk::tclvalue(maskEntryValue) <- maskEntryDefault
		} else {
			# Disable
			tcltk::tkconfigure(maskEntry, state="disabled")
			tcltk::tclvalue(maskEntryValue) <- ""
		}
	}
	
	replicateAction <- function() {
		if(tcltk::tclvalue(replicateCheckValue) == "1") {
			# Enable
			tcltk::tkconfigure(replicateEntry, state="normal")
			tcltk::tclvalue(replicateFunValue) <- "median"
		} else {
			# Disable
			tcltk::tkconfigure(replicateEntry, state="disabled")
			tcltk::tclvalue(replicateFunValue) <- ""
		}
	}
	
	segmentAction <- function() {
		if(tcltk::tclvalue(segValue) == "1") {
			# Enable widgets
			tcltk::tkconfigure(cbsCombo, state="readonly")
			tcltk::tkconfigure(cbsText, state="normal")
			tcltk::tkconfigure(gapsCheck, state="normal")
			tcltk::tkconfigure(modelizeCheck, state="normal")
			
			# Reset values to default
			tcltk::tclvalue(cbsComboValue) <- "fast"
			tcltk::tclvalue(gapsValue) <- "1"
			tcltk::tclvalue(modelizeValue) <- "1"
		} else {
			# Disable widgets
			tcltk::tkconfigure(cbsCombo, state="disabled")
			tcltk::tkconfigure(cbsText, state="disabled")
			tcltk::tkconfigure(gapsCheck, state="disabled")
			tcltk::tkconfigure(modelizeCheck, state="disabled")
			
			# Empty values
			tcltk::tclvalue(cbsComboValue) <- ""
			tcltk::tclvalue(gapsValue) <- "0"
			tcltk::tclvalue(modelizeValue) <- "0"
		}
		
		# Propagation
		cbsComboAction()
		modelizeAction()
	}
	
	modelizeAction <- function() {
		if(tcltk::tclvalue(modelizeValue) == "1") {
			# Enable
			tcltk::tkconfigure(modelEntry, state="normal")
			tcltk::tkconfigure(fittestSegmCheck, state="normal")
			tcltk::tclvalue(modelValue) <- process.default("modelizeArgs")
			tcltk::tclvalue(fittestSegmValue) <- "1"
		} else {
			# Disable
			tcltk::tkconfigure(modelEntry, state="disabled")
			tcltk::tkconfigure(fittestSegmCheck, state="disabled")
			tcltk::tclvalue(modelValue) <- ""
			tcltk::tclvalue(fittestSegmValue) <- "0"
		}
	}
	
	cbsComboAction <- function() {
		profile <- tcltk::tclvalue(cbsComboValue)
		if(profile == "custom") {
			tcltk::tkconfigure(cbsText, state="normal")
		} else {
			tcltk::tkconfigure(cbsText, state="normal")
			tcltk::tkdelete(cbsText, "1.0", "end")
			if(profile != "") tcltk::tkinsert(cbsText, "end", paste(process.default("segmentArgs", profile), collapse="\n"))
			tcltk::tkconfigure(cbsText, state="disabled")
		}
	}
	
	inputError <- function(message, icon="error") {
		tcltk::tkmessageBox(
			parent = localTopLevel,
			icon = icon,
			type = "ok",
			title = "Input error",
			message = message
		)
	}
	
	processAction = function() {
		
		# Extract from interface
		designFile <- tcltk::tclvalue(designValue)
		parserName <- tcltk::tclvalue(parserValue)
		outDirectory <- tcltk::tclvalue(outputValue)
		replicateFun <- get(tcltk::tclvalue(replicateFunValue))
		modelizeArgs <- tcltk::tclvalue(modelValue)
		
		# Data check
		if(!exists(parserName, mode="function")) return(inputError("The probe file parser must be the name of an existing function"))
		
		# Extract segmentation args
		segmentArgs <- as.character(tcltk::tcl(cbsText, "dump", text="1.0", "end"))[ c(FALSE, TRUE, FALSE) ]
		segmentArgs <- gsub("\n", "", segmentArgs)
		segmentArgs <- segmentArgs[ segmentArgs != "" ]
		
		# Extract cluster specifications
		cluster <- tcltk::tclvalue(clusterValue)
		if(cluster == "1") { cluster <- FALSE
		} else             { cluster <- list(spec=as.integer(cluster))
		}
		
		logFile <- sprintf("%s/process.log", outDirectory)
		
		
		# Probe processing
		steps <- "parse"
		if(tcltk::tclvalue(maskCheckValue) == "1")      steps <- c(steps, "mask")
		if(tcltk::tclvalue(replicateCheckValue) == "1") steps <- c(steps, "replicates")
		if(tcltk::tclvalue(wacaValue) == "1")           steps <- c(steps, "waca")
		steps <- c(steps, "export")
		
		# Segmentation
		if(tcltk::tclvalue(spatialValue) == "1")        steps <- c(steps, "spatial")
		if(tcltk::tclvalue(segValue) == "1")            steps <- c(steps, "segment")
		if(tcltk::tclvalue(gapsValue) == "1")           steps <- c(steps, "fill")
		if(tcltk::tclvalue(modelizeValue) == "1")       steps <- c(steps, "modelize")
		if(tcltk::tclvalue(segValue) == "1")            steps <- c(steps, "export")
		
		# Modelization
		if(tcltk::tclvalue(fittestSegmValue) == "1")    steps <- c(steps, "fittest", "export")
		if(tcltk::tclvalue(modelizeValue) == "1")       steps <- c(steps, "applyModel", "export")
		
		
		
		# Agilent.probes arguments
		columns <- tcltk::tclvalue(maskEntryValue)
		columns <- strsplit(columns, split=" *, *")[[1]]
		columns <- columns[ columns != "" ]
		names(columns) <- sprintf("flag_%s", columns)
		columns <- c(rFin="rProcessedSignal", gFin="gProcessedSignal", columns)
		probeArgs <- list(columns=columns)
		
		# Cursor
		tcltk::tkconfigure(localTopLevel, cursor="watch")
		tcltk::tcl("update")
		
		# Processing (catching what occurs ahead of process.log)
		handle(
			expr = {
				process(
					inputs = probeFiles,
					cluster = cluster,
					steps = steps,
					design = designFile,
					outDirectory = outDirectory,
					probeParser = parserName,
					probeArgs = probeArgs,
					replicateFun = replicateFun,
					segmentArgs = segmentArgs,
					modelizeArgs = modelizeArgs
				)
			},
			messageHandler = NULL,
			warningHandler = function(w) {
				tcltk::tkmessageBox(parent=localTopLevel, icon="warning", type="ok", title="R warning in process()", message=conditionMessage(w))
			},
			errorHandler = function(e) {
				tcltk::tkmessageBox(parent=localTopLevel, icon="error", type="ok", title="R error in process()", message=conditionMessage(e))
			}
		)
		
		# Cursor
		tcltk::tkconfigure(localTopLevel, cursor="arrow")
		
		# Open log file
		sts <- system2("open", logFile, stdout=FALSE, stderr=FALSE)                     ### Windows & Mac OS
		if(sts != 0L) sts <- system2("xdg-open", logFile, stdout=FALSE, stderr=FALSE)   ### Linux
	}
	
	
	## INTERFACE ##
	
	embedded <- !missing(globalTopLevel) && !missing(localTopLevel)
	if(!embedded) {
		# Linux default style
		if(.Platform$OS.type == "unix") try(tcltk::tcl("ttk::style", "theme", "use", "clam"), silent=TRUE)
		
		# Top level
		globalTopLevel <- localTopLevel <- tcltk::tktoplevel(class="cghRA")
		tcltk::tktitle(localTopLevel) <- "cghRA - Array processing"
		icon16 <- tcltk::tcl("image", "create", "photo", file=system.file("cghRA_16x16.gif", package="cghRA"))
		icon32 <- tcltk::tcl("image", "create", "photo", file=system.file("cghRA_32x32.gif", package="cghRA"))
		tcltk::tcl("wm", "iconphoto", localTopLevel, "-default", icon16, icon32)
		
		# Restrict to horizontal resizing
		tcltk::tkwm.resizable(localTopLevel, 1, 0)
	}
	
	# Evolutive width
	tcltk::tkgrid.columnconfigure(localTopLevel, 1, weight=1)
	
		# Preprocessing frame
		preproFrame <- tcltk::ttklabelframe(parent=localTopLevel, relief="groove", borderwidth=2, text="Preprocessing")
		tcltk::tkgrid(preproFrame, column=1, columnspan=2, row=1, padx=5, pady=5, sticky="we")
		tcltk::tkgrid.columnconfigure(preproFrame, 3, weight=1)
		
			# Probe files
			probeFiles <- character(0)
			probeFilesValue <- tcltk::tclVar("0 file(s) selected")
			probeFilesButton <- tcltk::tkbutton(parent=preproFrame, text="Probe files", command=probeFilesBrowse, width=20)
			probeFilesLabel <- tcltk::tklabel(parent=preproFrame, textvariable=probeFilesValue)
			tcltk::tkgrid(probeFilesButton, column=1, row=1, padx=5)
			tcltk::tkgrid(probeFilesLabel, column=2, columnspan=2, row=1, padx=5, sticky="w")
			
			# Processing button
			processButton <- tcltk::tkbutton(parent=preproFrame, text="Process arrays", command=processAction, background="#A9BEE1", font="Verdana 8 bold")
			tcltk::tkgrid(processButton, column=3, row=1, padx=c(5,10), sticky="e")
			
			# Parser combobox
			parserValue <- tcltk::tclVar("Agilent.probes")
			parserLabel <- tcltk::tklabel(parent=preproFrame, text="Probe file parser", width=20)
			parserCombo <- tcltk::ttkcombobox(parent=preproFrame, values=c("Agilent.probes", "custom.probes"), textvariable=parserValue, width=20, justify="center")
			tcltk::tkgrid(parserLabel, column=1, row=2, padx=5)
			tcltk::tkgrid(parserCombo, column=2, columnspan=2, row=2, padx=5, sticky="w")
			
			# Output directory
			outputValue <- tcltk::tclVar(".")
			outputButton <- tcltk::tkbutton(parent=preproFrame, text="Output directory", command=outputBrowse, width=20)
			outputEntry <- tcltk::tkentry(parent=preproFrame, textvariable=outputValue)
			tcltk::tkgrid(outputButton, column=1, row=3, padx=5)
			tcltk::tkgrid(outputEntry, column=2, columnspan=2, row=3, padx=5, sticky="we")
			
			# Design file
			designValue <- tcltk::tclVar("")
			designButton <- tcltk::tkbutton(parent=preproFrame, text="Design file", command=designBrowse, width=20)
			designEntry <- tcltk::tkentry(parent=preproFrame, textvariable=designValue)
			tcltk::tkgrid(designButton, column=1, row=4, padx=5)
			tcltk::tkgrid(designEntry, column=2, columnspan=2, row=4, padx=5, sticky="we")
			
			# Mask by flag
			maskCheckValue <- tcltk::tclVar("1")
			maskEntryDefault <- "rIsSaturated, gIsSaturated, rIsFeatNonUnifOL, gIsFeatNonUnifOL, rIsBGNonUnifOL, gIsBGNonUnifOL, rIsFeatPopnOL, gIsFeatPopnOL, rIsBGPopnOL, gIsBGPopnOL"
			maskEntryValue <- tcltk::tclVar(maskEntryDefault)
			maskLabel <- tcltk::tklabel(parent=preproFrame, text="Mask flagged probes", width=20)
			maskCheck <- tcltk::tkcheckbutton(parent=preproFrame, variable=maskCheckValue, command=maskAction)
			maskEntry <- tcltk::tkentry(parent=preproFrame, textvariable=maskEntryValue)
			tcltk::tkgrid(maskLabel, column=1, row=5, padx=5)
			tcltk::tkgrid(maskCheck, column=2, row=5, padx=0, sticky="w")
			tcltk::tkgrid(maskEntry, column=3, row=5, padx=5, sticky="we")
			
			# Replicates
			replicateCheckValue <- tcltk::tclVar("1")
			replicateFunValue <- tcltk::tclVar("median")
			replicateLabel <- tcltk::tklabel(parent=preproFrame, text="Replicated probes", width=20)
			replicateCheck <- tcltk::tkcheckbutton(parent=preproFrame, variable=replicateCheckValue, command=replicateAction)
			replicateEntry <- tcltk::tkentry(parent=preproFrame, textvariable=replicateFunValue, width=15, justify="center")
			tcltk::tkgrid(replicateLabel, column=1, row=6, padx=5)
			tcltk::tkgrid(replicateCheck, column=2, row=6, padx=0, sticky="w")
			tcltk::tkgrid(replicateEntry, column=3, row=6, padx=5, sticky="w")
			
			# WACA check
			wacaValue <- tcltk::tclVar("1")
			wacaLabel <- tcltk::tklabel(parent=preproFrame, text="WACA correction", width=20)
			wacaCheck <- tcltk::tkcheckbutton(parent=preproFrame, variable=wacaValue)
			tcltk::tkgrid(wacaLabel, column=1, row=7, padx=5)
			tcltk::tkgrid(wacaCheck, column=2, row=7, padx=0, sticky="w")
			
			# Spatial check
			spatialValue <- tcltk::tclVar("1")
			spatialLabel <- tcltk::tklabel(parent=preproFrame, text="Spatial distributions", width=20)
			spatialCheck <- tcltk::tkcheckbutton(parent=preproFrame, variable=spatialValue)
			tcltk::tkgrid(spatialLabel, column=1, row=8, padx=5)
			tcltk::tkgrid(spatialCheck, column=2, row=8, padx=0, sticky="w")
			
			# Cluster parameters
			clusterValue <- tcltk::tclVar( if(length(find.package("parallel", quiet=TRUE)) == 1) { parallel::detectCores() } else { NA } )
			clusterLabel <- tcltk::tklabel(parent=preproFrame, text="CPU cores", width=20)
			clusterEntry <- tcltk::tkentry(parent=preproFrame, textvariable=clusterValue, width=5, justify="center")
			tcltk::tkgrid(clusterLabel, column=1, row=9, padx=5, pady=c(0,5))
			tcltk::tkgrid(clusterEntry, column=2, columnspan=2, row=9, padx=5, pady=c(0,5), sticky="w")
		
		# Segmentation frame
		segFrame <- tcltk::ttklabelframe(parent=localTopLevel, relief="groove", borderwidth=2, text="Segmentation")
		tcltk::tkgrid(segFrame, column=1, columnspan=2, row=2, padx=5, pady=5, sticky="we")
		tcltk::tkgrid.columnconfigure(segFrame, 2, weight=1)
		
			# Segment check
			segValue <- tcltk::tclVar("1")
			segLabel <- tcltk::tklabel(parent=segFrame, text="Segment", width=20)
			segCheck <- tcltk::tkcheckbutton(parent=segFrame, variable=segValue, command=segmentAction)
			tcltk::tkgrid(segLabel, column=1, row=1, padx=5)
			tcltk::tkgrid(segCheck, column=2, row=1, padx=0, sticky="w")
			
			# fillGaps check
			gapsValue <- tcltk::tclVar("1")
			gapsLabel <- tcltk::tklabel(parent=segFrame, text="Fill inter-segment gaps", width=20)
			gapsCheck <- tcltk::tkcheckbutton(parent=segFrame, variable=gapsValue)
			tcltk::tkgrid(gapsLabel, column=1, row=2, padx=5)
			tcltk::tkgrid(gapsCheck, column=2, row=2, padx=0, sticky="w")
			
			# Fittest segmentation
			cbsLabel <- tcltk::tklabel(parent=segFrame, text="DNAcopy profiles", width=20)
			cbsFrame <- tcltk::tkframe(parent=segFrame)
			tcltk::tkgrid(cbsLabel, column=1, row=3, padx=5, sticky="n")
			tcltk::tkgrid(cbsFrame, column=2, row=3, sticky="we")
			tcltk::tkgrid.columnconfigure(cbsFrame, 5, weight=1)
			
				# Profile combobox
				cbsComboValue <- tcltk::tclVar("fast")
				cbsCombo <- tcltk::ttkcombobox(parent=cbsFrame, values=c(process.default()$segmentArgs, "custom"), textvariable=cbsComboValue, width=20, justify="center", state="readonly")
				tcltk::tkgrid(cbsCombo, column=1, columnspan=2, row=1, padx=5, sticky="w")
				
				# Details
				cbsScroll <- tcltk::tkscrollbar(parent=cbsFrame, repeatinterval=5, command=function(...) { tcltk::tkyview(cbsText,...) })
				tcltk::tkgrid(cbsScroll, column=2, row=2, padx=c(1, 5), pady=c(0,5), sticky="nsew")
				cbsText <- tcltk::tktext(parent=cbsFrame, height=3, yscrollcommand=function(...) { tcltk::tkset(cbsScroll, ...) })
				tcltk::tkinsert(cbsText, "end", paste(process.default("segmentArgs"), collapse="\n"))
				tcltk::tkgrid(cbsText, column=1, row=2, padx=c(5, 1), pady=c(0,5), sticky="nsew")
			
		# Auto model frame
		modelFrame <- tcltk::ttklabelframe(parent=localTopLevel, relief="groove", borderwidth=2, text="Copy number")
		tcltk::tkgrid(modelFrame, column=1, row=3, padx=5, pady=5, sticky="we")
		tcltk::tkgrid.columnconfigure(modelFrame, 2, weight=1)
		
			# Auto model check
			modelizeValue <- tcltk::tclVar("1")
			modelizeLabel <- tcltk::tklabel(parent=modelFrame, text="Modelize copy numbers", width=20)
			modelizeCheck <- tcltk::tkcheckbutton(parent=modelFrame, variable=modelizeValue, command=modelizeAction)
			tcltk::tkgrid(modelizeLabel, column=1, row=1, padx=5)
			tcltk::tkgrid(modelizeCheck, column=2, row=1, padx=0, sticky="w")
			
			# Fittest segmentation
			fittestLabel <- tcltk::tklabel(parent=modelFrame, text="Fittest segmentation", width=20)
			fittestSegmValue <- tcltk::tclVar("1")
			fittestSegmCheck <- tcltk::tkcheckbutton(parent=modelFrame, variable=fittestSegmValue)
			tcltk::tkgrid(fittestLabel, column=1, row=2, padx=5)
			tcltk::tkgrid(fittestSegmCheck, column=2, row=2, padx=0, sticky="w")
			
			# Auto model arguments
			modelValue <- tcltk::tclVar(process.default("modelizeArgs"))
			modelLabel <- tcltk::tklabel(parent=modelFrame, text="Custom parameters", width=20)
			modelEntry <- tcltk::tkentry(parent=modelFrame, textvariable=modelValue)
			tcltk::tkgrid(modelLabel, column=1, row=3, padx=5, pady=c(0,5))
			tcltk::tkgrid(modelEntry, column=2, row=3, padx=5, pady=c(0,5), sticky="we")
	
	# Events
	tcltk::tkbind(cbsCombo, "<<ComboboxSelected>>", cbsComboAction)
}

