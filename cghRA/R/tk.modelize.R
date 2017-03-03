# Interactive Tcl-Tk copy number modelization
# Author : Sylvain Mareschal <mareschal@ovsa.fr>

tk.modelize = function(
		compress = "gzip",
		compression_level = 9,
		exclude = c("X", "Y", "Xp", "Xq", "Yp", "Yq"),
		globalTopLevel,
		localTopLevel,
		render = c("auto", "png", "tkrplot"),
		tkrplot.scale = 1,
		png.res = 100,
		png.file = tempfile(fileext=".png")
		)
	{
	# Check renderer
	render <- match.arg(render)
	tcltkVer <- as.double(sub("^([0-9]+\\.[0-9]+).+$", "\\1", tcltk::tclVersion()))
	if(tcltkVer >= 8.6) {
		# PNG is not compatible
		if(render == "auto") render <- "png"
	} else {
		# PNG is compatible
		if(render == "auto")       { render <- "tkrplot"
		} else if(render == "png") { stop("PNG rendering requires tcltk 8.6 or above")
		}
	}
	
	# Check if tkrplot is installed
	if(length(find.package("tkrplot", quiet=TRUE)) == 0L) {
		tcltk::tkmessageBox(
			icon = "info",
			type = "ok",
			title = "tkrplot package missing",
			message = "The optional 'tkrplot' package is required to use tk.browse(), please install it and try again."
		)
	}
	
	# Current file
	arrayFiles <- character(0)
	extension <- as.character(NA)
	index <- 0L
	
	# Regions storage
	segStarts <- NULL
	segEnds <- NULL
	segChroms <- NULL
	segLogRatios <- NULL
	segProbes <- NULL
	
	# Model storage
	model <- NA
	
	# cghRA.regions storage
	regions <- NULL
	
	
	## FUNCTIONS ##
	
	savePar <- list()
	xConvert = function(x) {
		x <- as.integer(x)
		
		# Plot area coordinates in pixels
		xMin <- savePar$plt[1] * width
		xMax <- savePar$plt[2] * width
		
		# From pixel area to x range
		return((x - xMin) / (xMax - xMin) * (savePar$usr[2] - savePar$usr[1]) + savePar$usr[1])
	}
	yConvert = function(y) {
		y <- as.integer(y)
		
		# Plot area coordinates in pixels
		yMin <- (1-savePar$plt[3]) * height
		yMax <- (1-savePar$plt[4]) * height
		
		# From pixel area to x range
		return((y - yMin) / (yMax - yMin) * (savePar$usr[4] - savePar$usr[3]) + savePar$usr[3])
	}
	
	closest <- function(x, y, all.x, all.y) {
		distances <- abs((all.x - x) / max(all.x)) + abs((all.y - y) / max(all.y))
		return(which.min(distances))
	}
	
	click.x <- NULL
	click.y <- NULL
	click.id <- 0
	mousePress = function(x, y) {
		if(index != 0L) {
			# Nearest point on the plot
			target <- closest(xConvert(x), yConvert(y), LCN(segLogRatios, exact=TRUE), segProbes)
			if(target == click.id) {
				# Deselect
				click.x <<- NULL
				click.y <<- NULL
				click.id <<- 0
				
				# Update plot
				replot()
			} else {		
				# Copy number, if available
				if(!is.na(model['center']) && !is.na(model['width'])) {
					copies <- copies(segLogRatios[target], model=model)
				} else {
					copies <- NA
				}
				
				# Retain click coordinates
				click.x <<- LCN(segLogRatios, exact=TRUE)[ target ]
				click.y <<- segProbes[ target ]
				click.id <<- target
				
				# Update plot
				replot()
				
				# Info box
				tcltk::tkmessageBox(
					parent = localTopLevel,
					icon = "info",
					type = "ok",
					title = "Nearest segment",
					message = sprintf(
						"logRatio :\t%.5f\nCopies :\t%.2f\n\nProbes :\t%d\nChrom :\t%s\nStart :\t%.6f Mb\nEnd :\t%.6f Mb",
						segLogRatios[target],
						copies,
						segProbes[target],
						segChroms[target],
						segStarts[target] / 1e6,
						segEnds[target] / 1e6		
					)
				)
			}
		}
	}
	
	arrayFilesBrowse <- function() {
		# Update file list
		tmpArrayFiles <- tk.file(
			title = "Choose a region list",
			typeNames = c("cghRA region file", "Tab-separated region list"),
			typeExt = c(".regions.rdt", ".txt"),
			multiple = TRUE,
			mandatory = FALSE,
			type = "open",
			parent = localTopLevel
		)
		
		if(length(tmpArrayFiles) > 0) {
			# Data update
			arrayFiles <<- tmpArrayFiles
			
			# Interface update
			tcltk::tclvalue(arrayFilesValue) <- sprintf("1/%i selected files", length(arrayFiles))
			
			# Start with the new first
			index <<- 1L
			importFile()
		}
	}
	
	updateSliders = function() {
		# bandwidth
		if(is.na(model['bw'])) {
			tcltk::tclvalue(bwSliderValue) <- 0
			tcltk::tclvalue(bwSliderText) <- "NA"
		} else {
			tcltk::tclvalue(bwSliderValue) <- model['bw']
			tcltk::tclvalue(bwSliderText) <- as.character(round(model['bw'], 3))
		}
		
		# center
		if(is.na(model['center'])) {
			tcltk::tclvalue(centerSliderValue) <- 0
			tcltk::tclvalue(centerSliderText) <- "NA"
		} else {
			tcltk::tclvalue(centerSliderValue) <- model['center']
			tcltk::tclvalue(centerSliderText) <- as.character(round(model['center'], 3))
		}
		
		# width
		if(is.na(model['width'])) {
			tcltk::tclvalue(widthSliderValue) <- 0
			tcltk::tclvalue(widthSliderText) <- "NA"
		} else {
			tcltk::tclvalue(widthSliderValue) <- model['width']
			tcltk::tclvalue(widthSliderText) <- as.character(round(model['width'], 3))
		}
		
		# peakFrom
		if(is.na(model['peakFrom'])) {
			tcltk::tclvalue(fromSliderValue) <- -10
			tcltk::tclvalue(fromSliderText) <- "NA"
		} else {
			tcltk::tclvalue(fromSliderValue) <- model['peakFrom']
			tcltk::tclvalue(fromSliderText) <- as.character(round(model['peakFrom'], 2))
		}
		
		# peakTo
		if(is.na(model['peakTo'])) {
			tcltk::tclvalue(toSliderValue) <- 1.5
			tcltk::tclvalue(toSliderText) <- "NA"
		} else {
			tcltk::tclvalue(toSliderValue) <- model['peakTo']
			tcltk::tclvalue(toSliderText) <- as.character(round(model['peakTo'], 2))
		}
	}
	
	updateBw <- function(...) {
		if(index > 0) {
			newValue <- as.double(list(...)[[1]])
			if(newValue == 0) {
				model['bw'] <<- NA
				tcltk::tclvalue(bwSliderText) <- "NA"
			} else {
				model['bw'] <<- newValue
				tcltk::tclvalue(bwSliderText) <- as.character(newValue)
			}
			replot()
		}
	}
	
	updateCenter <- function(...) {
		if(index > 0) {
			newValue <- as.double(list(...)[[1]])
			if(newValue == 0) {
				model['center'] <<- NA
				tcltk::tclvalue(centerSliderText) <- "NA"
			} else {
				model['center'] <<- newValue
				tcltk::tclvalue(centerSliderText) <- as.character(newValue)
			}
			replot()
		}
	}
	
	updateWidth <- function(...) {
		if(index > 0) {
			newValue <- as.double(list(...)[[1]])
			if(newValue == 0) {
				model['width'] <<- NA
				tcltk::tclvalue(widthSliderText) <- "NA"
			} else {
				model['width'] <<- newValue
				tcltk::tclvalue(widthSliderText) <- as.character(newValue)
			}
			replot()
		}
	}
	
	updateFrom <- function(...) {
		if(index > 0) {
			newValue <- as.double(list(...)[[1]])
			if(newValue == -10) {
				model['peakFrom'] <<- NA
				tcltk::tclvalue(fromSliderText) <- "NA"
			} else {
				model['peakFrom'] <<- newValue
				tcltk::tclvalue(fromSliderText) <- as.character(newValue)
			}
			replot()
		}
	}
	
	updateTo <- function(...) {
		if(index > 0) {
			newValue <- as.double(list(...)[[1]])
			if(newValue == 1.5) {
				model['peakTo'] <<- NA
				tcltk::tclvalue(toSliderText) <- "NA"
			} else {
				model['peakTo'] <<- newValue
				tcltk::tclvalue(toSliderText) <- as.character(newValue)
			}
			replot()
		}
	}
	
	# Compute current plot area height, in pixels
	autoHeight <- function() {
		out <- as.integer(tcltk::tclvalue(tcltk::tkwinfo("height", plotFrame))) - 30L
		return(out)
	}
	
	# Compute current plot area width, in pixels
	autoWidth <- function() {
		out <- as.integer(tcltk::tclvalue(tcltk::tkwinfo("width", plotFrame))) - 10L
		return(out)
	}
	
	# model.test() call to produce the plot
	plot.core <- function() {
		graphics::par(bg="#FFFFFF", cex=1)
		savePar <<- model.test(
			segLogRatios = segLogRatios,
			segChroms = segChroms,
			segLengths = segProbes,
			minDensity = as.double(tcltk::tclvalue(minDensityValue)),
			model = model,
			returnPar = TRUE,
			exclude = exclude,
			title = basename(arrayFiles[index])
		)
		if(!is.null(click.x)) graphics::points(x=click.x, y=click.y, pch=1, col="#FF0000")
	}
	
	# Welcome screen
	plot.empty <- function() {
		graphics::par(bg="#FFFFFF", mar=c(0,0,0,0))
		graphics::plot(x=NA, y=NA, xlim=0:1, ylim=0:1, xlab="", ylab="", xaxt="n", yaxt="n", bty="n")
		graphics::text(x=0.5, y=0.5, labels="Welcome to cghRA !\n\nClick \"Select files\" and select *.regions.rdt files to begin.")
	}
	
	# Replot using 'png' rendered
	plot.png <- function(empty) {
		# Produce image file
		grDevices::png(png.file, width=width, height=height, res=png.res)
		if(isTRUE(empty)) { plot.empty()
		} else            { plot.core()
		}
		grDevices::dev.off()
		
		# Refresh image
		tcltk::tkconfigure(plotImage, file=png.file, width=width, height=height)
		tcltk::tkconfigure(plotWidget, width=width, height=height)
	}
	
	# Replot using 'tkrplot' rendered
	plot.tkrplot <- function(empty) {
		tkrplot::tkrreplot(
			lab = plotWidget,
			fun = if(isTRUE(empty)) { plot.empty } else { plot.core },
			hscale = hscale * tkrplot.scale,
			vscale = vscale * tkrplot.scale
		)
	}
	
	# Correct tkrplot scale factor
	changeScale <- function() {
		tkrplot.scale <<- tk.value(
			parent = NULL,
			type = "double",
			title = "Expansion factor (1 = 100%)",
			default = tkrplot.scale,
			allowEmpty = FALSE
		)
		replot()
	}
	
	# Pixel / tkrplot "scale" unit conversion factor
	scaleFactor <- NA
	
	# Replot common workflow
	replot <- function(empty=FALSE) {
		# Check coordinates
		if(!isTRUE(empty)) empty <- index == 0L
		
		# Guess scale factor from 2 x 1 empty plot (1 scale unit = x un-resized pixels)
		if(render == "tkrplot" && is.na(scaleFactor)) {
			scaleFactor <<- as.integer(tcltk::tclvalue(tcltk::tkwinfo("width", plotWidget))) / 2
			if(scaleFactor < 100) stop("Scale factor detection seems to have failed (", scaleFactor, ")")
		}
		
		# Adjust size
		height <<- autoHeight()
		width <<- autoWidth()
		if(render == "tkrplot") {
			vscale <<- height / scaleFactor
			hscale <<- width / scaleFactor
		}
		
		# Grab focus to avoid keyboard shortcuts quirks
		tcltk::tkfocus(plotWidget)
		
		# Replot
		handle(
			expr = {
				if(render == "png")            { plot.png(empty=empty)
				} else if(render == "tkrplot") { plot.tkrplot(empty=empty)
				}
			},
			# Silently ignore message()
			messageHandler = NULL,
			# Pass warning() but continue execution
			warningHandler = function(w) {
				tcltk::tkmessageBox(
					parent = localTopLevel,
					icon = "warning",
					type = "ok",
					title = "Warning in model.test()",
					message = conditionMessage(w)
				)
			},
			# Pass stop() and stop execution
			errorHandler = function(e) {
				tcltk::tkmessageBox(
					parent = localTopLevel, 
					icon = "error",
					type = "ok",
					title = "Error in model.test()",
					message = conditionMessage(e)
				)
			}					
		)
		
		# Adjust for Windows's magnifying factor
		if(render == "tkrplot" && tkrplot.scale > 1) {
			tcltk::tcl("update", "idletasks")
			tcltk::tkconfigure(plotWidget, width=width - 10L)
			tcltk::tkconfigure(plotWidget, height=height)
		}
	
		# Model was updated
		enableUpdate()
	}
	
	autoAction <- function() {
		model <<- model.auto(
			segLogRatios = segLogRatios,
			segChroms = segChroms,
			segLengths = segProbes,
			from = as.double(tcltk::tclvalue(fromValue)),
			to = as.double(tcltk::tclvalue(toValue)),
			by = as.double(tcltk::tclvalue(byValue)),
			precision = as.integer(tcltk::tclvalue(precisionValue)),
			maxPeaks = as.integer(tcltk::tclvalue(maxPeaksValue)),
			minWidth = as.double(tcltk::tclvalue(minWidthValue)),
			maxWidth = as.double(tcltk::tclvalue(maxWidthValue)),
			minDensity = as.double(tcltk::tclvalue(minDensityValue)),
			peakFrom = model['peakFrom'],
			peakTo = model['peakTo'],
			ploidy = as.double(tcltk::tclvalue(ploidyValue)),
			method = tcltk::tclvalue(methodValue),
			discreet = TRUE,
			exclude = exclude
		)
		updateSliders()
		replot()
	}
	
	importRdt = function() {
		if(length(find.package("cghRA", quiet=TRUE)) > 0) {
			# File loading
			regions <<- Rgb::readRDT(arrayFiles[index])
			if(is(regions, "cghRA.regions")) {
				if(regions$rowCount > 0) {
					# Data extraction
					segStarts <<- regions$extract(, "start")
					segEnds <<- regions$extract(, "end")
					segChroms <<- regions$extract(, "chrom")
					segLogRatios <<- regions$extract(, "logRatio")
					segProbes <<- regions$extract(, "probes")
					model <<- regions$model
					
					# Interface update
					updateSliders()
					replot()
					
					return(TRUE)
				} else {
					tcltk::tkmessageBox(parent=localTopLevel, icon="error", type="ok", title="Importation error", message="The cghRA.regions object contains no region")
					return(FALSE)
				}
			} else {
				tcltk::tkmessageBox(parent=localTopLevel, icon="error", type="ok", title="Importation error", message="The .rdt file does not contain a valid cghRA.regions object")
				return(FALSE)
			}
		} else {
			tcltk::tkmessageBox(parent=localTopLevel, icon="error", type="ok", title="Importation error", message=".rdt files need the 'cghRA' package")
			return(FALSE)
		}
	}
	
	exportRdt <- function() {
		if(!is.null(regions)) {
			# Model update
			regions$model <<- model
			
			# Save regions to file
			saveRDT(regions, arrayFiles[index])
			
			# Save copies to file
			if(regions$modelized()) {
				copies <- regions$model.apply()
				saveRDT(copies, sub("\\.regions\\.rdt$", ".copies.rdt", arrayFiles[index]))
			} else {
				unlink(sub("\\.regions\\.rdt$", ".copies.rdt", arrayFiles[index]))
			}
			
			return(TRUE)
		} else {
			tcltk::tkmessageBox(parent=localTopLevel, icon="error", type="ok", title="Exportation error", message="No valid cghRA.regions object in memory")
			return(FALSE)
		}
	}
	
	importTxt <- function() {
		tab <- try(utils::read.table(file=arrayFiles[index], sep="\t", dec=".", header=TRUE, stringsAsFactors=FALSE, comment.char="#"), silent=TRUE)
		if(!is(tab, "try-error")) {
			if(identical(names(tab), c("chrom", "start", "end", "probes", "logRatio"))) {
				# Region extraction
				segChroms <<- tab$chrom
				segStarts <<- tab$start
				segEnds <<- tab$end
				segLogRatios <<- tab$logRatio
				segProbes <<- tab$probes
				
				# Model extraction
				modelLine <- scan(file=arrayFiles[index], what="", sep="\n", n=1, quiet=TRUE)
				if(grepl("^#[A-Za-z]+=(([0-9\\.e-]+)|(NA))(, [A-Za-z]+=(([0-9\\.e-]+)|(NA)))*$", modelLine)) {
					modelLine <- substr(modelLine, 2, nchar(modelLine))
					modelLine <- strsplit(strsplit(modelLine, split=", ")[[1]], split="=")
					newModel <- suppressWarnings(as.numeric(sapply(modelLine, "[", 2)))
					names(newModel) <- sapply(modelLine, "[", 1)
					model <<- newModel
				} else {
					model <<- NA
				}
				
				# Interface update
				updateSliders()
				replot()
				
				return(TRUE)
			} else {
				tcltk::tkmessageBox(parent=localTopLevel, icon="error", type="ok", title="Importation error", message="Columns in tab-separated file must be :\n\nchrom\nstart\nend\nlogRatio\nprobes")
				return(FALSE)
			}
		} else {
			tcltk::tkmessageBox(parent=localTopLevel, icon="error", type="ok", title="Importation error", message=sprintf("Unable to read tab-separated file :\n\n%s", conditionMessage(attr(tab, "condition"))))
			return(FALSE)
		}
	}
	
	exportTxt <- function() {
		# Model
		modelLine <- sprintf("#%s", paste(sprintf("%s=%g", names(model), model), collapse=", "))
		
		# Save regions to file
		write(modelLine, file=arrayFiles[index], append=FALSE)
		suppressWarnings(
			utils::write.table(
				data.frame(
					chrom = segChroms,
					start = segStarts,
					end = segEnds,
					logRatio = segLogRatios,
					probes = segProbes,
					stringsAsFactors = FALSE
				),
				file = arrayFiles[index],
				sep = "\t",
				dec = ".",
				col.names = TRUE,
				row.names = FALSE,
				append = TRUE
			)
		)
		
		# Save copies to file
		if(!is.na(model['center']) && !is.na(model['width'])) {
			fileName <- sub("\\.txt$", ".copies.txt", arrayFiles[index])
			write(modelLine, file=fileName, append=FALSE)
			suppressWarnings(
				utils::write.table(
					model.apply(
						segStarts = segStarts,
						segEnds = segEnds,
						segChroms = segChroms,
						segLogRatios = segLogRatios,
						segLengths = segProbes,
						model = model,
						exact = FALSE,
						merge = TRUE
					),
					file = fileName,
					sep = "\t",
					dec = ".",
					col.names = TRUE,
					row.names = FALSE,
					append = TRUE
				)
			)
		}
		
		return(TRUE)
	}	
	
	importFile = function() {
		# File type
		extension <<- sub("^.*(\\.[^\\.]+)$", "\\1", arrayFiles[index])
		if(extension == ".rdt")        { status <- importRdt()
		} else if(extension == ".txt") { status <- importTxt()
		} else                         { tcltk::tkmessageBox(parent=localTopLevel, icon="error", type="ok", title="Importation error", message=sprintf("Unknown file extension \"%s\" (.rdt or .txt expected)", extension))
		}
		
		# Status
		if(status) { disableUpdate()
		} else     { tcltk::tkmessageBox(parent=localTopLevel, icon="error", type="ok", title="Importation error", message="Unable to import file, may be corrupted or with a wrong file extension")
		}
	}
	
	exportFile <- function() {
		# Cursor
		tcltk::tkconfigure(localTopLevel, cursor="watch")
		
		# Export
		if(extension == ".rdt")        { exportRdt()
		} else if(extension == ".txt") { exportTxt()
		} else                         { stop("Unexpected error in exportFile()")
		}
		
		# No need to update again
		disableUpdate()
		
		# Cursor
		tcltk::tkconfigure(localTopLevel, cursor="arrow")
	}
	
	nextFile <- function() {
		# Index update
		if(index >= length(arrayFiles)) {
			tcltk::tkmessageBox(parent=localTopLevel, icon="info", type="ok", title="Series complete", message="There is no next file to load")
		} else {
			# Deselect
			click.x <<- NULL
			click.y <<- NULL
			click.id <<- 0
			
			# Update
			index <<- index + 1L
			tcltk::tclvalue(arrayFilesValue) <- sprintf("%i/%i selected files", index, length(arrayFiles))
			importFile()
		}
	}
	
	previousFile <- function() {
		# Index update
		if(index <= 1) {
			tcltk::tkmessageBox(parent=localTopLevel, icon="info", type="ok", title="Series complete", message="There is no previous file to load")
		} else {
			# Deselect
			click.x <<- NULL
			click.y <<- NULL
			click.id <<- 0
			
			# Update
			index <<- index - 1L
			tcltk::tclvalue(arrayFilesValue) <- sprintf("%i/%i selected files", index, length(arrayFiles))
			importFile()
		}
	}
	
	enableUpdate <- function() {
		tcltk::tkconfigure(updateFileButton, background="#CC6666")
	}
	
	disableUpdate <- function() {
		tcltk::tkconfigure(updateFileButton, background=defaultBackground)
	}
	
	
	## INTERFACE ##
	
	embedded <- !missing(globalTopLevel) && !missing(localTopLevel)
	if(!embedded) {
		# Linux default style
		if(.Platform$OS.type == "unix") try(tcltk::tcl("ttk::style", "theme", "use", "clam"), silent=TRUE)
		
		# Top level
		globalTopLevel <- localTopLevel <- tcltk::tktoplevel(class="cghRA")
		tcltk::tktitle(localTopLevel) <- "cghRA - Modelization"
		icon16 <- tcltk::tcl("image", "create", "photo", file=system.file("cghRA_16x16.gif", package="cghRA"))
		icon32 <- tcltk::tcl("image", "create", "photo", file=system.file("cghRA_32x32.gif", package="cghRA"))
		tcltk::tcl("wm", "iconphoto", localTopLevel, "-default", icon16, icon32)
	}
	
	# Evolutive width
	tcltk::tkgrid.columnconfigure(localTopLevel, 1, weight=1)
	tcltk::tkgrid.rowconfigure(localTopLevel, 1, weight=1)
		
		# Plot frame
		plotFrame <- tcltk::ttklabelframe(parent=localTopLevel, relief="groove", borderwidth=2, text="Graphical representation")
		tcltk::tkgrid(plotFrame, column=1, columnspan=2, row=1, padx=5, pady=5, sticky="nsew")
		
			# R-Plot widget (wait for maximization to apply and propagate)
			tcltk::.Tcl("update idletasks")
			
			# Default size
			height <- 300L
			width <- autoWidth()
			
			if(render == "png") {
				# Display (empty) PNG image
				plotImage <- tcltk::tkimage.create("photo", width=width, height=height)
				plotWidget <- tcltk::tkcanvas(plotFrame, width=width, height=height)
				tcltk::tkcreate(plotWidget, "image", 0, 0, anchor="nw", image=plotImage)
			} else if(render == "tkrplot") {
				# tkrplot widget (fixed size image to guess scaleFactor)
				hscale <- 2
				vscale <- 0.6
				plotWidget <- tkrplot::tkrplot(parent=plotFrame, fun=plot.empty, hscale=hscale, vscale=vscale)
			}
			
			tcltk::tkgrid(plotWidget, column=1, row=1, padx=5, pady=5)
			
		# Slider frame
		sliderFrame <- tcltk::ttklabelframe(parent=localTopLevel, relief="groove", borderwidth=2, text="Model parameters")
			
			# Bandwidth slider
			bwSliderValue <- tcltk::tclVar("0")
			bwSliderText <- tcltk::tclVar("NA")
			bwSliderLabel1 <- tcltk::tklabel(sliderFrame, text="bandwidth", width=10)
			bwSliderLabel2 <- tcltk::tklabel(sliderFrame, textvariable=bwSliderText, width=10)
			bwSliderScale <- tcltk::tkscale(sliderFrame, from=0, to=0.5, showvalue=FALSE, variable=bwSliderValue, resolution=0.001, orient="horizontal", sliderlength=30, length=400, command=updateBw)
			tcltk::tkgrid(bwSliderLabel1, column=1, row=1, padx=2, pady=c(2,1))
			tcltk::tkgrid(bwSliderLabel2, column=2, row=1, padx=2, pady=1)
			tcltk::tkgrid(bwSliderScale, column=3, row=1, padx=2, pady=c(1,2), sticky="we")

			# Center slider
			centerSliderValue <- tcltk::tclVar("0")
			centerSliderText <- tcltk::tclVar("NA")
			centerSliderLabel1 <- tcltk::tklabel(sliderFrame, text="center", width=10)
			centerSliderLabel2 <- tcltk::tklabel(sliderFrame, textvariable=centerSliderText, width=10)
			centerSliderScale <- tcltk::tkscale(sliderFrame, from=0, to=4, showvalue=FALSE, variable=centerSliderValue, resolution=0.001, orient="horizontal", sliderlength=30, length=400, command=updateCenter)
			tcltk::tkgrid(centerSliderLabel1, column=1, row=2, padx=2, pady=c(2,1))
			tcltk::tkgrid(centerSliderLabel2, column=2, row=2, padx=2, pady=1)
			tcltk::tkgrid(centerSliderScale, column=3, row=2, padx=2, pady=c(1,2), sticky="we")
		
			# Width slider
			widthSliderValue <- tcltk::tclVar("0")
			widthSliderText <- tcltk::tclVar("NA")
			widthSliderLabel1 <- tcltk::tklabel(sliderFrame, text="width", width=10)
			widthSliderLabel2 <- tcltk::tklabel(sliderFrame, textvariable=widthSliderText, width=10)
			widthSliderScale <- tcltk::tkscale(sliderFrame, from=0, to=2, showvalue=FALSE, variable=widthSliderValue, resolution=0.001, orient="horizontal", sliderlength=30, length=400, command=updateWidth)
			tcltk::tkgrid(widthSliderLabel1, column=1, row=3, padx=2, pady=c(2,1))
			tcltk::tkgrid(widthSliderLabel2, column=2, row=3, padx=2, pady=1)
			tcltk::tkgrid(widthSliderScale, column=3, row=3, padx=2, pady=c(1,2), sticky="we")
		
		tcltk::tkgrid(sliderFrame, column=1, row=2, padx=5, pady=5, sticky="nsew")
		tcltk::tkgrid.columnconfigure(sliderFrame, 3, weight=1)
		
		# Peak range frame
		peakFrame <- tcltk::ttklabelframe(parent=localTopLevel, relief="groove", borderwidth=2, text="Peak range")
			
			# PeakFrom slider
			default <- paste(as.character(formals(model.auto)$peakFrom), collapse="")
			fromSliderValue <- tcltk::tclVar(default)
			fromSliderText <- tcltk::tclVar(default)
			fromSliderLabel1 <- tcltk::tklabel(peakFrame, text="peakFrom", width=10)
			fromSliderLabel2 <- tcltk::tklabel(peakFrame, textvariable=fromSliderText, width=10)
			fromSliderScale <- tcltk::tkscale(peakFrame, from=-10, to=1.5, showvalue=FALSE, variable=fromSliderValue, resolution=0.01, orient="horizontal", sliderlength=30, length=400, command=updateFrom)
			tcltk::tkgrid(fromSliderLabel1, column=1, row=1, padx=2, pady=c(2,1))
			tcltk::tkgrid(fromSliderLabel2, column=2, row=1, padx=2, pady=1)
			tcltk::tkgrid(fromSliderScale, column=3, row=1, padx=2, pady=c(1,2), sticky="we")
			
			# PeakTo slider
			default <- paste(as.character(formals(model.auto)$peakTo), collapse="")
			toSliderValue <- tcltk::tclVar(default)
			toSliderText <- tcltk::tclVar(default)
			toSliderLabel1 <- tcltk::tklabel(peakFrame, text="peakTo", width=10)
			toSliderLabel2 <- tcltk::tklabel(peakFrame, textvariable=toSliderText, width=10)
			toSliderScale <- tcltk::tkscale(peakFrame, from=-10, to=1.5, showvalue=FALSE, variable=toSliderValue, resolution=0.01, orient="horizontal", sliderlength=30, length=400, command=updateTo)
			tcltk::tkgrid(toSliderLabel1, column=1, row=2, padx=2, pady=c(2,1))
			tcltk::tkgrid(toSliderLabel2, column=2, row=2, padx=2, pady=1)
			tcltk::tkgrid(toSliderScale, column=3, row=2, padx=2, pady=c(1,2), sticky="we")
		
		tcltk::tkgrid(peakFrame, column=1, row=3, padx=5, pady=5, sticky="nsew")
		tcltk::tkgrid.columnconfigure(peakFrame, 3, weight=1)
		
		# Auto frame
		autoFrame <- tcltk::ttklabelframe(parent=localTopLevel, relief="groove", borderwidth=2, text="Auto model")
			
			# Auto model
			autoButton <- tcltk::tkbutton(parent=autoFrame, text="Modelize", command=autoAction, width=10)
			tcltk::tkgrid(autoButton, column=1, columnspan=3, row=1, padx=5, pady=c(10,15))
			
			# From
			fromValue <- tcltk::tclVar(formals(model.auto)$from)
			fromLabel <- tcltk::tklabel(parent=autoFrame, padx=5, text="bandwidth (from, to)")
			fromEntry <- tcltk::tkentry(parent=autoFrame, width=6, textvariable=fromValue)
			tcltk::tkgrid(fromLabel, column=1, row=2)
			tcltk::tkgrid(fromEntry, column=2, row=2)
			
			# To
			toValue <- tcltk::tclVar(formals(model.auto)$to)
			toEntry <- tcltk::tkentry(parent=autoFrame, width=6, textvariable=toValue)
			tcltk::tkgrid(toEntry, column=3, row=2, padx=c(0,5))
			
			# By
			byValue <- tcltk::tclVar(formals(model.auto)$by)
			byLabel <- tcltk::tklabel(parent=autoFrame, padx=5, text="bandwidth (by)")
			byEntry <- tcltk::tkentry(parent=autoFrame, width=6, textvariable=byValue)
			tcltk::tkgrid(byLabel, column=1, row=3)
			tcltk::tkgrid(byEntry, column=2, row=3)
			
			# Precision
			precisionValue <- tcltk::tclVar(formals(model.auto)$precision)
			precisionLabel <- tcltk::tklabel(parent=autoFrame, padx=5, text="precision")
			precisionEntry <- tcltk::tkentry(parent=autoFrame, width=6, textvariable=precisionValue)
			tcltk::tkgrid(precisionLabel, column=1, row=4)
			tcltk::tkgrid(precisionEntry, column=2, row=4)
			
			# Min density
			minDensityValue <- tcltk::tclVar(formals(model.auto)$minDensity)
			minDensityLabel <- tcltk::tklabel(parent=autoFrame, padx=5, text="minDensity*")
			minDensityEntry <- tcltk::tkentry(parent=autoFrame, width=6, textvariable=minDensityValue)
			tcltk::tkgrid(minDensityLabel, column=1, row=5)
			tcltk::tkgrid(minDensityEntry, column=2, row=5)
			
			# maxPeaks
			maxPeaksValue <- tcltk::tclVar(formals(model.auto)$maxPeaks)
			maxPeaksLabel <- tcltk::tklabel(parent=autoFrame, padx=5, text="maxPeaks")
			maxPeaksEntry <- tcltk::tkentry(parent=autoFrame, width=6, textvariable=maxPeaksValue)
			tcltk::tkgrid(maxPeaksLabel, column=1, row=6)
			tcltk::tkgrid(maxPeaksEntry, column=2, row=6)
			
			# minWidth
			minWidthValue <- tcltk::tclVar(formals(model.auto)$minWidth)
			minWidthLabel <- tcltk::tklabel(parent=autoFrame, padx=5, text="width (min, max)")
			minWidthEntry <- tcltk::tkentry(parent=autoFrame, width=6, textvariable=minWidthValue)
			tcltk::tkgrid(minWidthLabel, column=1, row=7)
			tcltk::tkgrid(minWidthEntry, column=2, row=7)
			
			# maxWidth
			maxWidthValue <- tcltk::tclVar(formals(model.auto)$maxWidth)
			maxWidthEntry <- tcltk::tkentry(parent=autoFrame, width=6, textvariable=maxWidthValue)
			tcltk::tkgrid(maxWidthEntry, column=3, row=7, padx=c(0,5))
			
			# Method
			methodValue <- tcltk::tclVar("stm")
			methodLabel <- tcltk::tklabel(parent=autoFrame, padx=5, text="method")
			methodEntry <- tcltk::tkentry(parent=autoFrame, width=6, textvariable=methodValue)
			tcltk::tkgrid(methodLabel, column=1, row=8)
			tcltk::tkgrid(methodEntry, column=2, row=8)
			
			# Ploidy
			ploidyValue <- tcltk::tclVar(formals(model.auto)$ploidy)
			ploidyLabel <- tcltk::tklabel(parent=autoFrame, padx=5, text="ploidy")
			ploidyEntry <- tcltk::tkentry(parent=autoFrame, width=6, textvariable=ploidyValue)
			tcltk::tkgrid(ploidyLabel, column=1, row=9)
			tcltk::tkgrid(ploidyEntry, column=2, row=9)
			
		tcltk::tkgrid(autoFrame, column=2, row=2, rowspan=3, padx=5, pady=5, sticky="nsew")
		
		# File frame
		fileFrame <- tcltk::ttklabelframe(parent=localTopLevel, relief="groove", borderwidth=2, text="Actions")
		tcltk::tkgrid.columnconfigure(fileFrame, 2, weight=1)
		tcltk::tkgrid.columnconfigure(fileFrame, 6, weight=1)
			
			# Resize button
			resizeButton <- tcltk::tkbutton(parent=fileFrame, text="Adjust plot size", command=changeScale)
			tcltk::tkgrid(resizeButton, column=1, row=1, padx=5, pady=5)
			
			# Previous file
			previousFileButton <- tcltk::tkbutton(parent=fileFrame, text="Previous", command=previousFile)
			tcltk::tkgrid(previousFileButton, column=3, row=1, padx=5, pady=5)
			
			# Array files
			arrayFilesValue <- tcltk::tclVar("Select files")
			arrayFilesButton <- tcltk::tkbutton(parent=fileFrame, textvariable=arrayFilesValue, command=arrayFilesBrowse, width=15)
			tcltk::tkgrid(arrayFilesButton, column=4, row=1, padx=5, pady=5)
			
			# Next file
			nextFileButton <- tcltk::tkbutton(parent=fileFrame, text="Next", command=nextFile)
			tcltk::tkgrid(nextFileButton, column=5, row=1, padx=5, pady=5)
			
			# Update file
			updateFileButton <- tcltk::tkbutton(parent=fileFrame, text="Update files", command=exportFile)
			tcltk::tkgrid(updateFileButton, column=7, row=1, padx=5, pady=5)
			defaultBackground <- as.character(tcltk::tkcget(updateFileButton, "-background"))
		
		tcltk::tkgrid(fileFrame, column=1, row=4, padx=5, pady=5, sticky="nsew")
	
	# Events
	tcltk::tkbind(plotWidget, "<ButtonPress-1>", mousePress)
}

