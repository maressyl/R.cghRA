\name{drawableFromClass.cghRA.probes}
\alias{drawableFromClass.cghRA.probes}
\title{
  Extend Rgb compatibility to cghRA.probes
}

\description{
  This function is only defined to allow the selection of RDT files containing \code{\linkS4class{cghRA.probes}} in Rgb \code{drawable.list}s. It should not be called directly by users.
}

\usage{
  drawableFromClass.cghRA.probes(track, design, ...)
}

\arguments{
  \item{track}{
    The \code{\linkS4class{cghRA.probes}} object extracted from the currently parsed RDT file.
}
  \item{design}{
    Either a \code{\linkS4class{cghRA.design}} matching \code{track} or the path to a RDT file containing it. Alternatively a Tcl-tk dialog window will be summoned to select such a RDT file if \code{design} was not set in the \code{drawable.list$add()} call.
}
  \item{\dots}{
    Further arguments are silently ignored.
}
}

\value{
  A \code{\linkS4class{cghRA.array}} object binding \code{track} and \code{design}.
}

\author{
  Sylvain Mareschal
}

\seealso{
  \code{\link{cghRA.array}}
}

