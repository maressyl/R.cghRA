\name{tk.value}
\alias{tk.value}

\title{
  Tk interface utilities
}

\description{
  This function prompt for a single value in a Tcl-tk interface.
}

\usage{
  tk.value(parent = NULL, type = c("character", "integer", "double"),
    title = "Enter a value", default = "", allowEmpty = FALSE)
}

\arguments{
  \item{parent}{
    Tcl-tk top-level to bind the popup window to.
}
  \item{type}{
    Single character value defining the type of the expected value.
}
  \item{title}{
    Single character value that will be displayed as the title of the popup window.
}
  \item{default}{
    Single value that will be used as default.
}
  \item{allowEmpty}{
    Single logical value, whether to raise an error if the user does not provide any value or not.
}
}

\value{
  Returns the entered value, casted to \code{type}.
}

\author{
  Sylvain Mareschal
}

