# Interactive Tcl-Tk design processing
# Author : Sylvain Mareschal <mareschal@ovsa.fr>

tk.design = function(
		organism = "Human",
		assembly = "GRCh37",
		chromosomes = "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y",
		chromFiles = "",
		restrictionSites = "AluI=AG|CT, RsaI=GT|AC",
		globalTopLevel,
		localTopLevel
		)
	{
	
	## FUNCTIONS ##

	inputFileBrowse <- function() {
		tcltk::tclvalue(inputFileValue) <- tk.file(
			title = "Choose a design file",
			typeNames = c("Agilent TDT file", "Any custom design file"),
			typeExt = c(".txt", ".*"),
			multiple = FALSE,
			mandatory = FALSE,
			type = "open",
			parent = localTopLevel
		)
	}

	outputDirBrowse <- function() {
		tcltk::tclvalue(outputDirValue) <- tk.folder(
			title = "Where to save the design file",
			mustexist = TRUE,
			mandatory = FALSE
		)
	}

	probeFileBrowse <- function() {
		tcltk::tclvalue(probeFileValue) <- tk.file(
			title = "Choose a probe sequence file",
			typeNames = "Fasta file",
			typeExt = ".fa .fst .fasta .txt",
			multiple = FALSE,
			mandatory = FALSE,
			type = "open",
			parent = localTopLevel
		)
	}
	
	bandFileBrowse <- function() {
		tcltk::tclvalue(bandFileValue) <- tk.file(
			title = "Choose a Cytoband annotation track",
			typeNames = "Cytoband track",
			typeExt = ".rdt",
			multiple = FALSE,
			mandatory = FALSE,
			type = "open",
			parent = localTopLevel
		)
	}
	
	chromFilesBrowse <- function() {
		chromFiles <<- tk.file(
			title = "Choose chromosome sequence files",
			typeNames = "Fasta file",
			typeExt = ".fa .fst .fasta .txt",
			multiple = TRUE,
			mandatory = FALSE,
			type = "open",
			parent = localTopLevel
		)
		tcltk::tclvalue(chromFilesValue) <- sprintf("%i file(s) selected", length(chromFiles))
	}

	remapCheck <- function(x) {
		if(tcltk::tclvalue(remapValue) == "1") {
			# Enable
			tcltk::tkconfigure(probeFileButton, state="normal")
			tcltk::tkconfigure(probeFileEntry, state="normal")
			tcltk::tkconfigure(chromFilesButton, state="normal")
		} else {
			# Disable if bias disabled
			if(tcltk::tclvalue(biasValue) != "1") {
				tcltk::tkconfigure(probeFileButton, state="disabled")
				tcltk::tkconfigure(probeFileEntry, state="disabled")
				tcltk::tkconfigure(chromFilesButton, state="disabled")
			}
		}
	}

	biasCheck <- function() {
		if(tcltk::tclvalue(biasValue) == "1") {
			# Enable
			tcltk::tkconfigure(sitesEntry, state="normal")
			tcltk::tkconfigure(probeFileButton, state="normal")
			tcltk::tkconfigure(probeFileEntry, state="normal")
			tcltk::tkconfigure(chromFilesButton, state="normal")
		} else {
			# Disable
			tcltk::tkconfigure(sitesEntry, state="disabled")
			if(tcltk::tclvalue(remapValue) != "1") {
				tcltk::tkconfigure(probeFileButton, state="disabled")
				tcltk::tkconfigure(probeFileEntry, state="disabled")
				tcltk::tkconfigure(chromFilesButton, state="disabled")
			}
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
	
	process = function() {
		
		# Data collection
		doRemap <- as.logical(as.integer(tcltk::tclvalue(remapValue)))
		doBias <- as.logical(as.integer(tcltk::tclvalue(biasValue)))
		parserName <- tcltk::tclvalue(parserValue)
		designFile <- tcltk::tclvalue(inputFileValue)
		bandFile <- tcltk::tclvalue(bandFileValue)
		outDir <- tcltk::tclvalue(outputDirValue)
		organism <- tcltk::tclvalue(organismValue)
		assembly <- tcltk::tclvalue(assemblyValue)
		chromosomes <- tcltk::tclvalue(chromValue)
		
		# Data check
		if(!file.exists(designFile))                return(inputError("Input file can not be found"))
		if(organism == "")                          return(inputError("Organism is empty"))
		if(assembly == "")                          return(inputError("Assembly is empty"))
		if(!grepl("^[^,]+(,[^,]+)*$", chromosomes)) return(inputError("Chromosome names must be separated by commas"))
		if(!exists(parserName, mode="function"))    return(inputError("The input file parser must be the name of an existing function"))
		
		# Data reshaping
		chromosomes <- strsplit(chromosomes, split=",")[[1]]
		
		if(doRemap) {
			# Data collection
			probeFile <- tcltk::tclvalue(probeFileValue)
			
			# Check probe file
			if(!file.exists(probeFile)) return(inputError("Probe file can not be found"))
		}
		
		if(doRemap || doBias) {
			# Chromosome file existence
			if(length(chromFiles) > 0) {
				for(f in chromFiles) {
					if(!file.exists(f)) return(inputError(sprintf("A chromosome file can not be found :\n%s", f)))
				}
			} else {
				return(inputError("No chromosome file provided"))
			}
			
			# Warn about name inconsistency
			chromNames <- sub("^(.+)\\.[^\\.]+$", "\\1", basename(chromFiles))
			if(!setequal(chromNames, chromosomes)) {
				# Build warning message
				message <- ""
				
				# Chromosome without file
				missingFile <- setdiff(chromosomes, chromNames)
				if(length(missingFile) > 0) message <- sprintf("%s%i chromosome(s) without FASTA file: %s\n\n", message, length(missingFile), paste(missingFile, collapse=", "))
				
				# File without chromosome
				missingChrom <- basename(chromFiles)[ ! chromNames %in% chromosomes ]
				if(length(missingChrom) > 0) message <- sprintf("%s%i FASTA file(s) of unknown chromosome: %s\n\n", message, length(missingChrom), paste(missingChrom, collapse=", "))
				
				# Reminder
				message <- sprintf("%sPlease remind that file names (without extension) must perfectly match provided chromosomes names.", message)
				
				if(length(intersect(chromNames, chromosomes)) == 0) {
					# No chromosome, stop here
					return(inputError(message, icon="error"))
				} else {
					# Some overlap, continue
					inputError(message, icon="warning")
				}
			}
		}
		
		if(doBias) {
			# Data collection
			fragSites <- tcltk::tclvalue(sitesValue)
			
			# Data check
			if(!grepl("^([^=]+=[ACGT]+\\|[ACGT]+)(, [^=]+=[ACGT]+\\|[ACGT]+)*$", fragSites)) return(inputError("Restriction sites must be provided as in the following example :\nAluI=G|CT, RsaI=GT|AC"))
		
			# Data reshaping
			sites <- strsplit(strsplit(fragSites, split=", ")[[1]], split="=")
			fragSites <- sapply(sites, "[", 2)
			names(fragSites) <- sapply(sites, "[", 1)
			
			# Default values
			digits <- 6
		}
		
		# Cursor
		tcltk::tkconfigure(localTopLevel, cursor="watch")
		tcltk::tcl("update")
		
		# Immediate warnings
		ops <- options(warn = 1)
		on.exit(options(ops))
		
		# Log
		logFile <- sprintf("%s/design.log", outDir)
		unlink(logFile)
		cat(sprintf("%s\n\n", Sys.time()), file=logFile, append=FALSE)
		cat(sprintf("R %s.%s (%s)\n", R.version$major, R.version$minor, R.version$platform), file=logFile, append=TRUE)
		for(package in c("Rgb", "cghRA")) {
			if(length(find.package(package, quiet=TRUE)) > 0) cat(sprintf("%-10s version %s\n", package, utils::packageVersion(package)), file=logFile, append=TRUE)
		}
		cat("\n\n", file=logFile, append=TRUE)
		
		# Processing
		handle(
			expr = {
				# Parse design file
				message("<Extraction>\n")
				args <- list(
					designFile,
					organism = organism,
					assembly = assembly,
					chromosomes = chromosomes
				)
				design <- do.call(what=parserName, args=args)
				if(!is(design, "cghRA.design")) stop("The input file parser did not returned a cghRA.design object")
				
				if(doRemap) {
					# Probe remapping
					message("\n<Remapping>\n")
					design$remap(
						probeFile = probeFile,
						chromFiles = chromFiles,
						chromPattern = "^(.+)\\.[^\\.]+$",
						blatArgs = character(0),
						rawOutput = FALSE,
						noMulti = TRUE,
						noOverlap = TRUE,
						noPartial = TRUE,
						verbose = 2
					)
				}
				
				if(doBias) {
					# WACA bias computation
					message("\n<WACA bias computation>\n")
					design$bias(
						chromFiles = chromFiles,
						chromPattern = "^(.+)\\.[^\\.]+$",
						fragSites = fragSites,
						digits = digits,
						verbose = 1
					)
				}
				
				# Centromere splitting
				if(file.exists(bandFile)) {
					bands <- readRDT(bandFile)
					bands$eraseArms()
					bands <- bands$extract(expression(stain == "acen"))
					centromeres <- tapply(bands$start, bands$chrom, max)
					design$addArms(centromeres)
				}
				
				# Write to .rdt file
				message("\n<Writing>\n")
				Rgb::saveRDT(design, sprintf("%s/%s.design.rdt", outDir, design$name))
				
				message("\n<Done>", appendLF=FALSE)
			},
			messageHandler = function(m) {
				# Log
				content <- conditionMessage(m)
				cat(content, file=logFile, append=TRUE)
			},
			warningHandler = function(w) {
				# Log
				content <- sprintf("\nWARNING : %s\n", conditionMessage(w))
				cat(content, file=logFile, append=TRUE)
			},
			errorHandler = function(e) {
				# Log
				content <- sprintf("ERROR : %s", conditionMessage(e))
				delim <- paste(rep("=", nchar(content)), collapse="")
				content <- sprintf("\n%s\n%s\n%s\n", delim, content, delim)
				cat(content, file=logFile, append=TRUE)
			}
		)
		
		# Cursor back to normal
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
		tcltk::tktitle(localTopLevel) <- "cghRA - Design processing"
		icon16 <- tcltk::tcl("image", "create", "photo", file=system.file("cghRA_16x16.gif", package="cghRA"))
		icon32 <- tcltk::tcl("image", "create", "photo", file=system.file("cghRA_32x32.gif", package="cghRA"))
		tcltk::tcl("wm", "iconphoto", localTopLevel, "-default", icon16, icon32)
		
		# Restrict to horizontal resizing
		tcltk::tkwm.resizable(localTopLevel, 1, 0)
	}
	
	# Evolutive width
	tcltk::tkgrid.columnconfigure(localTopLevel, 1, weight=1)
		
		# Core frame
		coreFrame <- tcltk::ttklabelframe(parent=localTopLevel, relief="groove", borderwidth=2, text="Extraction")
		tcltk::tkgrid(coreFrame, column=1, columnspan=2, row=1, padx=5, pady=5, sticky="we")
		tcltk::tkgrid.columnconfigure(coreFrame, 3, weight=1)
		
			# Organism
			organismValue <- tcltk::tclVar(organism)
			organismLabel <- tcltk::tklabel(parent=coreFrame, text="Organism", width=20)
			organismEntry <- tcltk::tkentry(parent=coreFrame, textvariable=organismValue, width=20)
			tcltk::tkgrid(organismLabel, column=1, row=1, padx=5, pady=5)
			tcltk::tkgrid(organismEntry, column=2, columnspan=2, row=1, padx=5, pady=5, sticky="w")
			
			# Processing button
			processButton <- tcltk::tkbutton(parent=coreFrame, text="Process design", command=process, background="#A9BEE1", font="Verdana 8 bold")
			tcltk::tkgrid(processButton, column=3, row=1, padx=c(5,10), pady=5, sticky="e")
			
			# Assembly
			assemblyValue <- tcltk::tclVar(assembly)
			assemblyLabel <- tcltk::tklabel(parent=coreFrame, text="Assembly", width=20)
			assemblyEntry <- tcltk::tkentry(parent=coreFrame, textvariable=assemblyValue, width=20)
			tcltk::tkgrid(assemblyLabel, column=1, row=2, padx=5, pady=5)
			tcltk::tkgrid(assemblyEntry, column=2, columnspan=2, row=2, padx=5, pady=5, sticky="w")
			
			# Chromosomes
			chromValue <- tcltk::tclVar(chromosomes)
			chromLabel <- tcltk::tklabel(parent=coreFrame, text="Chromosomes", width=20)
			chromEntry <- tcltk::tkentry(parent=coreFrame, textvariable=chromValue)
			tcltk::tkgrid(chromLabel, column=1, row=3, padx=5, pady=5)
			tcltk::tkgrid(chromEntry, column=2, columnspan=2, row=3, padx=5, pady=5, sticky="we")
			
			# Input file
			inputFileValue <- tcltk::tclVar("")
			inputFileButton <- tcltk::tkbutton(parent=coreFrame, text="Input file", command=inputFileBrowse, width=20)
			inputFileEntry <- tcltk::tkentry(parent=coreFrame, textvariable=inputFileValue)
			tcltk::tkgrid(inputFileButton, column=1, row=4, padx=5, pady=5)
			tcltk::tkgrid(inputFileEntry, column=2, columnspan=2, row=4, padx=5, pady=5, sticky="we")
			
			# Parser combobox
			parserValue <- tcltk::tclVar("Agilent.design")
			parserLabel <- tcltk::tklabel(parent=coreFrame, text="Input file parser", width=20)
			parserCombo <- tcltk::ttkcombobox(parent=coreFrame, values=c("Agilent.design", "custom.design"), textvariable=parserValue, width=20, justify="center")
			tcltk::tkgrid(parserLabel, column=1, row=5, padx=5, pady=5)
			tcltk::tkgrid(parserCombo, column=2, columnspan=2, row=5, padx=5, pady=5, sticky="w")
			
			# Output file
			outputDirValue <- tcltk::tclVar(".")
			outputDirButton <- tcltk::tkbutton(parent=coreFrame, text="Output directory", command=outputDirBrowse, width=20)
			outputDirEntry <- tcltk::tkentry(parent=coreFrame, textvariable=outputDirValue)
			tcltk::tkgrid(outputDirButton, column=1, row=6, padx=5, pady=5)
			tcltk::tkgrid(outputDirEntry, column=2, columnspan=2, row=6, padx=5, pady=5, sticky="we")
			
		# Update frame
		updateFrame <- tcltk::ttklabelframe(parent=localTopLevel, relief="groove", borderwidth=2, text="Enhancement")
		tcltk::tkgrid(updateFrame, column=1, columnspan=2, row=2, padx=5, pady=5, sticky="we")
		tcltk::tkgrid.columnconfigure(updateFrame, 2, weight=1)
		
			# Remap check
			remapValue <- tcltk::tclVar("1")
			remapLabel <- tcltk::tklabel(parent=updateFrame, text="Remap probes", width=20)
			remapCheck <- tcltk::tkcheckbutton(parent=updateFrame, variable=remapValue, command=remapCheck)
			tcltk::tkgrid(remapLabel, column=1, row=1, padx=5, pady=5)
			tcltk::tkgrid(remapCheck, column=2, row=1, padx=5, pady=5, sticky="w")
			
			# Bias check
			biasValue <- tcltk::tclVar("1")
			biasLabel <- tcltk::tklabel(parent=updateFrame, text="Compute WACA biases", width=20)
			biasCheck <- tcltk::tkcheckbutton(parent=updateFrame, variable=biasValue, command=biasCheck)
			tcltk::tkgrid(biasLabel, column=1, row=2, padx=5, pady=5)
			tcltk::tkgrid(biasCheck, column=2, row=2, padx=5, pady=5, sticky="w")
			
			# Bands file
			bandFileValue <- tcltk::tclVar("")
			bandFileButton <- tcltk::tkbutton(parent=updateFrame, text="Cytoband file (arm split)", command=bandFileBrowse, width=20)
			bandFileEntry <- tcltk::tkentry(parent=updateFrame, textvariable=bandFileValue)
			tcltk::tkgrid(bandFileButton, column=1, row=3, padx=5, pady=5)
			tcltk::tkgrid(bandFileEntry, column=2, row=3, padx=5, pady=5, sticky="we")
			
			# Probe file
			probeFileValue <- tcltk::tclVar("")
			probeFileButton <- tcltk::tkbutton(parent=updateFrame, text="Probe file", command=probeFileBrowse, width=20)
			probeFileEntry <- tcltk::tkentry(parent=updateFrame, textvariable=probeFileValue)
			tcltk::tkgrid(probeFileButton, column=1, row=4, padx=5, pady=5)
			tcltk::tkgrid(probeFileEntry, column=2, row=4, padx=5, pady=5, sticky="we")
			
			# Chrom files
			chromFiles <- character(0)
			chromFilesValue <- tcltk::tclVar("0 file(s) selected")
			chromFilesButton <- tcltk::tkbutton(parent=updateFrame, text="Chromosome files", command=chromFilesBrowse, width=20)
			chromFilesLabel <- tcltk::tklabel(parent=updateFrame, textvariable=chromFilesValue)
			tcltk::tkgrid(chromFilesButton, column=1, row=5, padx=5, pady=5)
			tcltk::tkgrid(chromFilesLabel, column=2, row=5, padx=5, pady=5, sticky="w")
			
			# Fragmentation sites
			sitesValue <- tcltk::tclVar(restrictionSites)
			sitesLabel <- tcltk::tklabel(parent=updateFrame, text="Restriction sites", width=20)
			sitesEntry <- tcltk::tkentry(parent=updateFrame, textvariable=sitesValue)
			tcltk::tkgrid(sitesLabel, column=1, row=6, padx=5, pady=5)
			tcltk::tkgrid(sitesEntry, column=2, row=6, padx=5, pady=5, sticky="we")
	
	# Default width
	if(!embedded) tcltk::tkwm.geometry(localTopLevel, sprintf("600x%i", as.integer(tcltk::tkwinfo("height", localTopLevel))))
		
}

