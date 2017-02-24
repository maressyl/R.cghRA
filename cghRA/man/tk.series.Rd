\name{tk.series}
\alias{tk.series}
\title{
  Interactive cghRA series processing
}

\description{
  This function provides a Tcl-Tk interface to perform series analysis on processed arrays and designs.
}

\usage{
  tk.series(globalTopLevel, localTopLevel)
}

\arguments{
  \item{globalTopLevel}{
    This argument should be filled only when embedding this Tcl-Tk interface in an other. It is the top level of the embedding interface, generally a call to \code{\link{tktoplevel}}.
}
  \item{localTopLevel}{
    This argument should be filled only when embedding this Tcl-Tk interface in an other. It is the local top level to use to build this interface, generally a \code{\link{tkframe}} or \code{\link{ttkframe}}.
}
}

\author{
  Sylvain Mareschal
}

\seealso{
  \code{\link{tk.cghRA}}
}
