# Interactive Tcl-Tk series processing
# Author : Sylvain Mareschal <mareschal@ovsa.fr>

tk.series = function(
		globalTopLevel,
		localTopLevel
		)
	{
	
	## FUNCTIONS ##

	regionFilesBrowse <- function() {
		regionFiles <<- tk.file(
			title = "Choose cghRA processed region files",
			typeNames = c("cghRA modelized region files", "cghRA raw region files"),
			typeExt = c(".copies.rdt", ".regions.rdt"),
			multiple = TRUE,
			mandatory = FALSE,
			type = "open",
			parent = localTopLevel
		)
		tcltk::tclvalue(regionFilesValue) <- sprintf("%i file(s) selected", length(regionFiles))
	}

	outputDirBrowse <- function() {
		tcltk::tclvalue(outputDirValue) <- tk.folder(
			title = "Where to save the produced files",
			mustexist = TRUE,
			mandatory = FALSE
		)
	}
	
	mcrFileBrowse <- function() {
		tcltk::tclvalue(mcrFileValue) <- tk.file(
			title = "Choose a penetrance file",
			typeNames = "Penetrance file",
			typeExt = ".rdt",
			multiple = FALSE,
			mandatory = FALSE,
			type = "open",
			parent = localTopLevel
		)
	}
	
	bandTrackBrowse <- function() {
		tcltk::tclvalue(bandTrackValue) <- tk.file(
			title = "Choose a Cytoband annotation track",
			typeNames = "Cytoband track",
			typeExt = ".rdt",
			multiple = FALSE,
			mandatory = FALSE,
			type = "open",
			parent = localTopLevel
		)
	}
	
	commonData <- function() {
		# Data collection
		name <- tcltk::tclvalue(nameValue)
		output <- tcltk::tclvalue(outputDirValue)
		rawStates <- tcltk::tclvalue(statesValue)
		value <- tcltk::tclvalue(statesTypeValue)
		
		# Data check
		if(name == "")                                           return(tcltk::tkmessageBox(parent=localTopLevel, icon="error", type="ok", title="Input error", message=sprintf("Series name is empty")))
		if(!file.exists(output) || !file.info(output)['isdir'])  return(tcltk::tkmessageBox(parent=localTopLevel, icon="error", type="ok", title="Input error", message=sprintf("Output directory can not be found")))
		if(length(regionFiles) > 0) {
			for(f in regionFiles) {
				if(!file.exists(f))                              return(tcltk::tkmessageBox(parent=localTopLevel, icon="error", type="ok", title="Processing error", message=sprintf("A region file can not be found :\n%s", f)))
			}
		} else {
			return(tcltk::tkmessageBox(parent=localTopLevel, icon="error", type="ok", title="Processing error", message=sprintf("No region file provided")))
		}
		
		# States
		rawStates <- strsplit(strsplit(rawStates, split=" ")[[1]], split="[\\(\\),]")
		if(any(sapply(rawStates, length) != 3))                  return(tcltk::tkmessageBox(parent=localTopLevel, icon="error", type="ok", title="Input error", message="Invalid copy level description"))
		states <- lapply(lapply(rawStates, "[", 2:3), as.numeric)
		names(states) <- sapply(rawStates, "[", 1)
		
		# Processing
		tcltk::tkconfigure(localTopLevel, cursor="watch")
		handle(
			expr = {
				# Series aggregation
				series <- cghRA.series(.name=name)
				for(regionFile in regionFiles) {
					series$add(Rgb::readRDT(regionFile))
				}
				series$check()
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
		tcltk::tkconfigure(localTopLevel, cursor="arrow")
		
		return(
			list(
				name = name,
				series = series,
				output = output,
				states = states,
				value = value
			)
		)
	}
	
	trackProcess <- function() {
		# Common data check
		commons <- commonData()
		if(is.list(commons)) {
			# Specific data
			para <- as.logical(as.integer(tcltk::tclvalue(trackParaValue)))
			pene <- as.logical(as.integer(tcltk::tclvalue(trackPeneValue)))
			pool <- as.logical(as.integer(tcltk::tclvalue(trackPoolValue)))
			plot <- as.logical(as.integer(tcltk::tclvalue(trackPlotPlotValue)))
			width <- as.integer(tcltk::tclvalue(trackPlotWidthValue))
			height <- as.integer(tcltk::tclvalue(trackPlotHeightValue))
			res <- as.integer(tcltk::tclvalue(trackPlotResValue))
			
			# Processing
			tcltk::tkconfigure(localTopLevel, cursor="watch")
			handle(
				expr = {
					# Parallelization
					if(para) {
						parallelization <- commons$series$parallelize(value=commons$value)
						saveRDT(parallelization, sprintf("%s/%s.rdt", commons$output, parallelization$name))
					}
					
					# Penetrance
					if(pene) {
						# Produce tracks
						penetrance <- commons$series$penetrance(value=commons$value, states=commons$states)
						for(p in penetrance) saveRDT(p, sprintf("%s/%s.rdt", commons$output, p$name))
						
						# Summary plot
						if(plot) {
							# Build drawable list
							dl <- drawable.list()
							for(p in penetrance) dl$add(file=NA, track=p)
							
							# Plot
							grDevices::png(sprintf("%s/%s penetrance.png", commons$output, commons$series$name), width=width, height=height, res=res)
							singlePlot(dl)
							grDevices::dev.off()							
						}
					}
					
					# Pool
					if(pool) {
						# Produce track
						poolTrack <- commons$series$pool(value=commons$value, states=commons$states)
						saveRDT(poolTrack, sprintf("%s/%s.rdt", commons$output, poolTrack$name))
						
						# Summary plot
						if(plot) {
							# Build drawable list
							dl <- drawable.list()
							dl$add(file=NA, track=poolTrack)
							
							# Plot
							grDevices::png(sprintf("%s/%s pool.png", commons$output, commons$series$name), width=width, height=height, res=res)
							singlePlot(dl)
							grDevices::dev.off()							
						}
					}
					
					# Finished
					tcltk::tkmessageBox(parent=localTopLevel, icon="info", type="ok", title="Done", message="Track production complete")
				},
				messageHandler = NULL,
				warningHandler = function(w) {
					tcltk::tkmessageBox(parent=localTopLevel, icon="warning", type="ok", title="R warning", message=conditionMessage(w))
				},
				errorHandler = function(e) {
					tcltk::tkmessageBox(parent=localTopLevel, icon="error", type="ok", title="R error", message=conditionMessage(e))
				}
			)
			tcltk::tkconfigure(localTopLevel, cursor="arrow")
		}
	}
	
	stepsProcess <- function() {
		# Common data check
		commons <- commonData()
		if(is.list(commons)) {
			# Specific data
			dpen <- as.numeric(tcltk::tclvalue(stepsDpenValue))
			vpen <- as.numeric(tcltk::tclvalue(stepsVpenValue))
			gpen <- as.numeric(tcltk::tclvalue(stepsGpenValue))
			nested <- tcltk::tclvalue(stepsNestedValue)
			
			# Processing
			tcltk::tkconfigure(localTopLevel, cursor="watch")
			handle(
				expr = {
					# MCR computation
					STEPS <- commons$series$STEPS(
						value = commons$value,
						states = commons$states,
						dpen = dpen,
						vpen = vpen,
						gpen = gpen,
						threshold = NA,
						nested = nested
					)
					
					# RDT files
					for(state in names(STEPS)) {
						saveRDT(STEPS[[ state ]], sprintf("%s/%s.rdt", commons$output, STEPS[[ state ]]$name))
					}
					
					# Finished
					tcltk::tkmessageBox(parent=localTopLevel, icon="info", type="ok", title="Done", message="STEPS computation complete")
				},
				messageHandler = NULL,
				warningHandler = function(w) {
					tcltk::tkmessageBox(parent=localTopLevel, icon="warning", type="ok", title="R warning", message=conditionMessage(w))
				},
				errorHandler = function(e) {
					tcltk::tkmessageBox(parent=localTopLevel, icon="error", type="ok", title="R error", message=conditionMessage(e))
				}
			)
			tcltk::tkconfigure(localTopLevel, cursor="arrow")
		}
	}
	
	sraProcess <- function() {
		# Common data check
		commons <- commonData()
		if(is.list(commons)) {
			# Specific data
			sra <- as.logical(as.integer(tcltk::tclvalue(lenzSraValue)))
			lra <- as.logical(as.integer(tcltk::tclvalue(lenzLraValue)))
			if(sra || lra) {
				# Processing
				tcltk::tkconfigure(localTopLevel, cursor="watch")
				handle(
					expr = {
						if(sra) {
							# SRA computation
							SRA <- commons$series$SRA(
								value = commons$value,
								states = commons$states,
								quiet = TRUE
							)
							
							# RDT files
							for(state in names(SRA)) {
								saveRDT(SRA[[ state ]], sprintf("%s/%s.rdt", commons$output, SRA[[ state ]]$name))
							}
						}
						
						if(lra) {
							# LRA computation
							LRA <- commons$series$LRA(
								value = commons$value,
								states = commons$states,
								quiet = TRUE
							)
							
							# RDT files
							for(state in names(LRA)) {
								saveRDT(LRA[[ state ]], sprintf("%s/%s.rdt", commons$output, LRA[[ state ]]$name))
							}
						}
						
						# Finished
						tcltk::tkmessageBox(parent=localTopLevel, icon="info", type="ok", title="Done", message="Reccurent Abnormalities computation complete")
					},
					messageHandler = NULL,
					warningHandler = function(w) {
						tcltk::tkmessageBox(parent=localTopLevel, icon="warning", type="ok", title="R warning", message=conditionMessage(w))
					},
					errorHandler = function(e) {
						tcltk::tkmessageBox(parent=localTopLevel, icon="error", type="ok", title="R error", message=conditionMessage(e))
					}
				)
				tcltk::tkconfigure(localTopLevel, cursor="arrow")
			}
		}
	}
	


	## INTERFACE ##

	embedded <- !missing(globalTopLevel) && !missing(localTopLevel)
	if(!embedded) {
		# Linux default style
		if(.Platform$OS.type == "unix") try(tcltk::tcl("ttk::style", "theme", "use", "clam"), silent=TRUE)
		
		# Top level
		globalTopLevel <- localTopLevel <- tcltk::tktoplevel(class="cghRA")
		tcltk::tktitle(localTopLevel) <- "cghRA - Series processing"
		icon16 <- tcltk::tcl("image", "create", "photo", file=system.file("cghRA_16x16.gif", package="cghRA"))
		icon32 <- tcltk::tcl("image", "create", "photo", file=system.file("cghRA_32x32.gif", package="cghRA"))
		tcltk::tcl("wm", "iconphoto", localTopLevel, "-default", icon16, icon32)
	}
	
	# Evolutive width
	tcltk::tkgrid.columnconfigure(localTopLevel, 1, weight=1)
	
		# Input frame
		inputFrame <- tcltk::ttklabelframe(parent=localTopLevel, relief="groove", borderwidth=2, text="Series to analyze")
		tcltk::tkgrid(inputFrame, column=1, row=1, padx=5, pady=5, sticky="we")
		tcltk::tkgrid.columnconfigure(inputFrame, 2, weight=1)
			
			# Name input
			nameValue <- tcltk::tclVar("")
			nameLabel <- tcltk::tklabel(parent=inputFrame, text="Series name", width=20)
			nameEntry <- tcltk::tkentry(parent=inputFrame, textvariable=nameValue)
			tcltk::tkgrid(nameLabel, column=1, row=1, padx=5, pady=5)
			tcltk::tkgrid(nameEntry, column=2, row=1, padx=5, pady=5, sticky="w")
			
			# Region files
			regionFiles <- character(0)
			regionFilesValue <- tcltk::tclVar("0 file(s) selected")
			regionFilesButton <- tcltk::tkbutton(parent=inputFrame, text="Region files", command=regionFilesBrowse, width=17)
			regionFilesLabel <- tcltk::tklabel(parent=inputFrame, textvariable=regionFilesValue)
			tcltk::tkgrid(regionFilesButton, column=1, row=2, padx=5, pady=5)
			tcltk::tkgrid(regionFilesLabel, column=2, row=2, padx=5, pady=5, sticky="w")
			
			# Output directory
			outputDirValue <- tcltk::tclVar(".")
			outputDirButton <- tcltk::tkbutton(parent=inputFrame, text="Output directory", command=outputDirBrowse, width=17)
			outputDirEntry <- tcltk::tkentry(parent=inputFrame, textvariable=outputDirValue)
			tcltk::tkgrid(outputDirButton, column=1, row=3, padx=5, pady=5)
			tcltk::tkgrid(outputDirEntry, column=2, columnspan=2, row=3, padx=5, pady=5, sticky="we")
			
			# States label
			statesValue <- tcltk::tclVar("deletion(-Inf,-0.5) gain(0.5,Inf)")
			statesLabel <- tcltk::tklabel(parent=inputFrame, text="Alteration intervals", width=20)
			statesEntry <- tcltk::tkentry(parent=inputFrame, textvariable=statesValue)
			tcltk::tkgrid(statesLabel, column=1, row=4, padx=5, pady=5)
			tcltk::tkgrid(statesEntry, column=2, row=4, padx=5, pady=5, sticky="we")
			
			# State type combobox
			statesTypeValue <- tcltk::tclVar("copies")
			statesCombo <- tcltk::ttkcombobox(parent=inputFrame, values=c("copies", "logRatio"), textvariable=statesTypeValue, justify="center", width=8, state="readonly")
			tcltk::tkgrid(statesCombo, column=3, row=4, padx=5, pady=5)
		
		# Track frame
		trackFrame <- tcltk::ttklabelframe(parent=localTopLevel, relief="groove", borderwidth=2, text="Produce Rgb track files")
		tcltk::tkgrid(trackFrame, column=1, row=2, padx=5, pady=5, sticky="we")
		tcltk::tkgrid.columnconfigure(trackFrame, 2, weight=1)
			
			# Parameters label
			trackParamLabel <- tcltk::tklabel(parent=trackFrame, text="Selection", width=20)
			trackParamFrame <- tcltk::tkframe(parent=trackFrame)
			tcltk::tkgrid(trackParamLabel, column=1, row=2, padx=5, pady=5)
			tcltk::tkgrid(trackParamFrame, column=2, row=2, sticky="we")
			tcltk::tkgrid.columnconfigure(trackParamFrame, 7, weight=1)
				
				# Parallelize checkbox
				trackParaValue <- tcltk::tclVar("0")
				trackParaCheck <- tcltk::tkcheckbutton(parent=trackParamFrame, variable=trackParaValue)
				trackParaLabel <- tcltk::tklabel(parent=trackParamFrame, text="Copy number matrix", anchor="w", justify="left")
				tcltk::tkgrid(trackParaCheck, column=1, row=1, padx=0, pady=5)
				tcltk::tkgrid(trackParaLabel, column=2, row=1, padx=c(0,20), pady=5)
				
				# Penetrance checkbox
				trackPeneValue <- tcltk::tclVar("1")
				trackPeneCheck <- tcltk::tkcheckbutton(parent=trackParamFrame, variable=trackPeneValue)
				trackPeneLabel <- tcltk::tklabel(parent=trackParamFrame, text="Penetrance", anchor="w", justify="left")
				tcltk::tkgrid(trackPeneCheck, column=3, row=1, padx=0, pady=5)
				tcltk::tkgrid(trackPeneLabel, column=4, row=1, padx=c(0,20), pady=5)
				
				# Pool checkbox
				trackPoolValue <- tcltk::tclVar("1")
				trackPoolCheck <- tcltk::tkcheckbutton(parent=trackParamFrame, variable=trackPoolValue)
				trackPoolLabel <- tcltk::tklabel(parent=trackParamFrame, text="Altered segment pool", anchor="w", justify="left")
				tcltk::tkgrid(trackPoolCheck, column=5, row=1, padx=0, pady=5)
				tcltk::tkgrid(trackPoolLabel, column=6, row=1, padx=c(0,20), pady=5)
				
				# Processing button
				trackButton <- tcltk::tkbutton(parent=trackParamFrame, text="Compute", background="#A9BEE1", font="Verdana 8 bold", command=trackProcess)
				tcltk::tkgrid(trackButton, column=8, row=1, padx=5, pady=5, sticky="e")
			
			# Parameters label
			trackPlotLabel <- tcltk::tklabel(parent=trackFrame, text="Summary plots", width=20)
			trackPlotFrame <- tcltk::tkframe(parent=trackFrame)
			tcltk::tkgrid(trackPlotLabel, column=1, row=3, padx=5, pady=5)
			tcltk::tkgrid(trackPlotFrame, column=2, row=3, sticky="we")
			tcltk::tkgrid.columnconfigure(trackParamFrame, 7, weight=1)
				
				# Produce plot checkbox
				trackPlotPlotValue <- tcltk::tclVar("1")
				trackPlotPlotCheck <- tcltk::tkcheckbutton(parent=trackPlotFrame, variable=trackPlotPlotValue)
				trackPlotPlotLabel <- tcltk::tklabel(parent=trackPlotFrame, text="Plot", anchor="w", justify="left")
				tcltk::tkgrid(trackPlotPlotCheck, column=1, row=1, padx=0, pady=5)
				tcltk::tkgrid(trackPlotPlotLabel, column=2, row=1, padx=c(0,20), pady=5)
				
				# Width and height inputs
				trackPlotWidthValue <- tcltk::tclVar("1600")
				trackPlotHeightValue <- tcltk::tclVar("1200")
				trackPlotWidthLabel <- tcltk::tklabel(parent=trackPlotFrame, text="size (px) :")
				trackPlotWidthEntry <- tcltk::tkentry(parent=trackPlotFrame, textvariable=trackPlotWidthValue, width=5, justify="center")
				trackPlotHeightLabel <- tcltk::tklabel(parent=trackPlotFrame, text="x")
				trackPlotHeightEntry <- tcltk::tkentry(parent=trackPlotFrame, textvariable=trackPlotHeightValue, width=5, justify="center")
				tcltk::tkgrid(trackPlotWidthLabel, column=3, row=1, padx=c(20,5), pady=5)
				tcltk::tkgrid(trackPlotWidthEntry, column=4, row=1, padx=5, pady=5, sticky="w")
				tcltk::tkgrid(trackPlotHeightLabel, column=5, row=1, padx=0, pady=5)
				tcltk::tkgrid(trackPlotHeightEntry, column=6, row=1, padx=5, pady=5, sticky="w")
				
				# Res input
				trackPlotResValue <- tcltk::tclVar("100")
				trackPlotResLabel <- tcltk::tklabel(parent=trackPlotFrame, text="res (ppi) :")
				trackPlotResEntry <- tcltk::tkentry(parent=trackPlotFrame, textvariable=trackPlotResValue, width=5, justify="center")
				tcltk::tkgrid(trackPlotResLabel, column=7, row=1, padx=c(20,5), pady=5)
				tcltk::tkgrid(trackPlotResEntry, column=8, row=1, padx=5, pady=5, sticky="w")
		
		# MCR frame
		mcrFrame <- tcltk::ttklabelframe(parent=localTopLevel, relief="groove", borderwidth=2, text="Compute regions of interest")
		tcltk::tkgrid(mcrFrame, column=1, row=3, padx=5, pady=5, sticky="we")
		tcltk::tkgrid.columnconfigure(mcrFrame, 2, weight=1)
			
			# STEPS
			stepsLabel <- tcltk::tklabel(parent=mcrFrame, text="STEPS", width=20)
			stepsFrame <- tcltk::tkframe(parent=mcrFrame)
			tcltk::tkgrid(stepsLabel, column=1, row=1, padx=5, pady=5)
			tcltk::tkgrid(stepsFrame, column=2, row=1, sticky="we")
			tcltk::tkgrid.columnconfigure(stepsFrame, 9, weight=1)
				
				# dpen input
				stepsDpenValue <- tcltk::tclVar("2")
				stepsDpenLabel <- tcltk::tklabel(parent=stepsFrame, text="dpen :")
				stepsDpenEntry <- tcltk::tkentry(parent=stepsFrame, textvariable=stepsDpenValue, width=5, justify="center")
				tcltk::tkgrid(stepsDpenLabel, column=1, row=1, padx=5, pady=5)
				tcltk::tkgrid(stepsDpenEntry, column=2, row=1, padx=5, pady=5, sticky="w")
				
				# dpen input
				stepsVpenValue <- tcltk::tclVar("0.8")
				stepsVpenLabel <- tcltk::tklabel(parent=stepsFrame, text="vpen :")
				stepsVpenEntry <- tcltk::tkentry(parent=stepsFrame, textvariable=stepsVpenValue, width=5, justify="center")
				tcltk::tkgrid(stepsVpenLabel, column=3, row=1, padx=c(20,5), pady=5)
				tcltk::tkgrid(stepsVpenEntry, column=4, row=1, padx=5, pady=5, sticky="w")
				
				# gpen input
				stepsGpenValue <- tcltk::tclVar("0.3")
				stepsGpenLabel <- tcltk::tklabel(parent=stepsFrame, text="gpen :")
				stepsGpenEntry <- tcltk::tkentry(parent=stepsFrame, textvariable=stepsGpenValue, width=5, justify="center")
				tcltk::tkgrid(stepsGpenLabel, column=5, row=1, padx=c(20,5), pady=5)
				tcltk::tkgrid(stepsGpenEntry, column=6, row=1, padx=5, pady=5, sticky="w")
				
				# nested combobox
				stepsNestedValue <- tcltk::tclVar("merge")
				stepsNestedLabel <- tcltk::tklabel(parent=stepsFrame, text="nested :")
				stepsNestedCombo <- tcltk::ttkcombobox(parent=stepsFrame, values=c("merge", "flag", "none"), textvariable=stepsNestedValue, justify="center", width=6, state="readonly")
				tcltk::tkgrid(stepsNestedLabel, column=7, row=1, padx=c(20,5), pady=5)
				tcltk::tkgrid(stepsNestedCombo, column=8, row=1, padx=5, pady=5, sticky="w")
				
				# Processing button
				stepsButton <- tcltk::tkbutton(parent=stepsFrame, text="Compute", background="#A9BEE1", font="Verdana 8 bold", command=stepsProcess)
				tcltk::tkgrid(stepsButton, column=10, row=1, padx=5, pady=5, sticky="e")
			
			# Lenz
			lenzLabel <- tcltk::tklabel(parent=mcrFrame, text="Lenz et al, PNAS 2008", width=20)
			lenzFrame <- tcltk::tkframe(parent=mcrFrame)
			tcltk::tkgrid(lenzLabel, column=1, row=2, padx=5, pady=5)
			tcltk::tkgrid(lenzFrame, column=2, row=2, sticky="we")
			tcltk::tkgrid.columnconfigure(lenzFrame, 5, weight=1)
				
				# SRA checkbox
				lenzSraValue <- tcltk::tclVar("1")
				lenzSraCheck <- tcltk::tkcheckbutton(parent=lenzFrame, variable=lenzSraValue)
				lenzSraLabel <- tcltk::tklabel(parent=lenzFrame, text="SRA", justify="left")
				tcltk::tkgrid(lenzSraCheck, column=1, row=1, padx=0, pady=5)
				tcltk::tkgrid(lenzSraLabel, column=2, row=1, padx=c(0,20), pady=5)
				
				# LRA checkbox
				lenzLraValue <- tcltk::tclVar("1")
				lenzLraCheck <- tcltk::tkcheckbutton(parent=lenzFrame, variable=lenzSraValue)
				lenzLraLabel <- tcltk::tklabel(parent=lenzFrame, text="LRA", justify="left")
				tcltk::tkgrid(lenzLraCheck, column=3, row=1, padx=0, pady=5)
				tcltk::tkgrid(lenzLraLabel, column=4, row=1, padx=c(0,20), pady=5)
			
				# Processing button
				lenzButton <- tcltk::tkbutton(parent=lenzFrame, text="Compute", background="#A9BEE1", font="Verdana 8 bold", command=sraProcess)
				tcltk::tkgrid(lenzButton, column=6, row=1, padx=5, pady=5, sticky="e")
			
}

