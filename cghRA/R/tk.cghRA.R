# Launcher for cghRA tools
# Author : Sylvain Mareschal <mareschal@ovsa.fr>

tk.cghRA = function(
		blocking = FALSE,
		tkrplot.scale = 1
		)
	{
	
	## FUNCTIONS ##
	
	tabChanged <- function() {
		index <- as.integer(tcltk::tclvalue(tcltk::tcl(notebook, "index", "current"))) + 1L
		tcltk::tkfocus(force=notebookTab[[index]])
	}
	
	
	## INTERFACE ##
	
	# Linux default style
	if(.Platform$OS.type == "unix") try(tcltk::tcl("ttk::style", "theme", "use", "clam"), silent=TRUE)
	
	# Top level
	topLevel <- tcltk::tktoplevel(class="cghRA")
	tcltk::tktitle(topLevel) <- "cghRA"
	icon16 <- tcltk::tcl("image", "create", "photo", file=system.file("cghRA_16x16.gif", package="cghRA"))
	icon32 <- tcltk::tcl("image", "create", "photo", file=system.file("cghRA_32x32.gif", package="cghRA"))
	tcltk::tcl("wm", "iconphoto", topLevel, "-default", icon16, icon32)
		
		# Notebook
		notebook <- tcltk::ttknotebook(topLevel)
		notebookTab <- list()
		tcltk::tkgrid(notebook, column=1, row=1, sticky="nsew")
			
			# Design
			notebookTab[[1]] <- tcltk::ttkframe(topLevel)
			tk.design(globalTopLevel=topLevel, localTopLevel=notebookTab[[1]])
			tcltk::tkadd(notebook, notebookTab[[1]], text="Design processing", sticky="nsew")
			
			# Process
			notebookTab[[2]] <- tcltk::ttkframe(topLevel)
			tk.process(globalTopLevel=topLevel, localTopLevel=notebookTab[[2]])
			tcltk::tkadd(notebook, notebookTab[[2]], text="Array processing", sticky="nsew")
			
			# Modelize
			notebookTab[[3]] <- tcltk::ttkframe(topLevel)
			modelizeTopLevel <- tk.modelize(globalTopLevel=topLevel, localTopLevel=notebookTab[[3]], tkrplot.scale=tkrplot.scale)
			tcltk::tkadd(notebook, notebookTab[[3]], text="Copy number modelization", sticky="nsew")
			
			# Annotate
			notebookTab[[4]] <- tcltk::ttkframe(topLevel)
			tk.annotate(globalTopLevel=topLevel, localTopLevel=notebookTab[[4]])
			tcltk::tkadd(notebook, notebookTab[[4]], text="Annotate regions", sticky="nsew")
			
			# Series
			notebookTab[[5]] <- tcltk::ttkframe(topLevel)
			tk.series(globalTopLevel=topLevel, localTopLevel=notebookTab[[5]])
			tcltk::tkadd(notebook, notebookTab[[5]], text="Series processing", sticky="nsew")
			
	# Horizontal resizing
	tcltk::tkgrid.columnconfigure(topLevel, 1, weight=1)
	tcltk::tkgrid.rowconfigure(topLevel, 1, weight=1)
	
	
	## LAUNCH ##
	
	# Events
	tcltk::tkbind(notebook, "<<NotebookTabChanged>>", tabChanged)
	
	# Wait for closing
	if(isTRUE(blocking)) tcltk::tkwait.window(topLevel)
}

