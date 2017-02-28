#!/bin/Rscript

### Author : Sylvain Mareschal <mareschal@ovsa.fr>
### Tcl-tk window to type a value

tk.value <- function(parent=NULL, type=c("character", "integer", "double"), title="Enter a value", default="", allowEmpty=FALSE) {
	
	# Check args
	type <- match.arg(type)
	
	# Command
	out <- character(0)
	doneCommand <- function() {
		# Get value
		out <<- tcltk::tclvalue(value)
		if(!isTRUE(allowEmpty) && out == "") {
			# Empty string
			tcltk::tkmessageBox(
				parent = topLevel, 
				icon = "error",
				type = "ok",
				title = "Incorrect value",
				message = "Please enter a value."
			)
			return()
		}
		
		if(type == "integer") {
			# Integer output
			out <<- tcltk::tclvalue(value)
			if(grepl("[^0-9]", out)) {
				# Not an integer
				tcltk::tkmessageBox(
					parent = topLevel, 
					icon = "error",
					type = "ok",
					title = "Incorrect value",
					message = "The provided value is not an integer."
				)
				return()
			}
			out <<- as.integer(out)			
		} else if(type == "double") {
			# Double output
			out <<- tcltk::tclvalue(value)
			out <<- sub(",", ".", out, fixed=TRUE)
			if(grepl("[^0-9\\.]", out)) {
				# Not an integer
				tcltk::tkmessageBox(
					parent = topLevel, 
					icon = "error",
					type = "ok",
					title = "Incorrect value",
					message = "The provided value is not a decimal number."
				)
				return()
			}
			out <<- as.double(out)			
		}
		
		# Return success
		tcltk::tkdestroy(topLevel)
	}
	
	# Top level
	topLevel <- tcltk::tktoplevel(class="Rgb")
	tcltk::tktitle(topLevel) <- "Prompt"
	
	# Make slave
	if(!is.null(parent)) {
		tcltk::tcl("wm", "transient", topLevel, parent)
		tcltk::tcl("wm", "withdraw", topLevel)
		tcltk::tcl("wm", "deiconify", topLevel)
	}
	
	# Resizable
	tcltk::tkgrid.columnconfigure(topLevel, 1, weight=1)
	tcltk::tkwm.resizable(topLevel, 1, 0)
	
	# Value
	value <- tcltk::tclVar(default)
	label <- tcltk::tklabel(parent=topLevel, text=title)
	entry <- tcltk::tkentry(parent=topLevel, width=30, textvariable=value, justify="center")
	tcltk::tkgrid(label, column=1, row=1, padx=10, pady=c(10,2), sticky="w")
	tcltk::tkgrid(entry, column=1, row=2, padx=10, pady=c(2,5), sticky="ew")
	tcltk::tkfocus(entry)
	
	# Done button
	doneButton <- tcltk::tkbutton(parent=topLevel, text="Done", width=10, command=doneCommand)
	tcltk::tkgrid(doneButton, column=1, row=3, padx=10, pady=c(5,10))
	
	# Wait
	tcltk::tkwait.window(topLevel)
	
	return(out)
}

