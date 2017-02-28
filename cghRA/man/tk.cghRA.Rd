\name{tk.cghRA}
\alias{tk.cghRA}
\title{
  cghRA Tcl-Tk launcher
}

\description{
  This function produces a Tcl-Tk interface merging all the cghRA components installed.
}

\usage{
  tk.cghRA(blocking = FALSE, tkrplot.scale = 1)
}

\arguments{
  \item{blocking}{
    Single logical value, whether to wait for the interface window to be closed before unfreezing the R console. The \code{FALSE} default let you use R and the interface in parallel, the code{TRUE} is used essentially in the stand alone version.
}
  \item{tkrplot.scale}{
    Single numeric value to be passed to \code{\link{tk.modelize}}.
}
}

\author{
  Sylvain Mareschal
}

\seealso{
  \code{\link{tk.design}}, \code{\link{tk.process}}, \code{\link{tk.modelize}}, \code{\link{tk.annotate}}, \code{\link{tk.series}}, \code{\link[Rgb]{tk.convert}}, \code{\link[Rgb]{tk.browse}}
}
