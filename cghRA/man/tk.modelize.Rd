\name{tk.modelize}
\alias{tk.modelize}
\title{
  Interactive copy number modelization
}

\description{
  This function provides a Tcl-Tk interface to produce or adjust a CGH copy number model on single or multiple arrays.
}

\usage{
  tk.modelize(compress = "gzip", compression_level = 9, exclude = c("X", "Y", "Xp", "Xq",
    "Yp", "Yq"), globalTopLevel, localTopLevel, scaleFactor = NA, hscale = 2, vscale = 1)
}

\arguments{
  \item{compress}{
    To be passed to \code{cghRA-class} \code{toRdat} method.
}
  \item{compression_level}{
    To be passed to \code{cghRA-class} \code{toRdat} method.
}
  \item{exclude}{
    Vector, the chromosomes to exclude from the density computation and to plot with distinct symbols (use \code{NULL} to disable this feature). Sexual chromosomes should be excluded in heterogeneous DNA source, as their desequilibrium (2 'X' and no 'Y' in women) impact normal cells AND tumoral ones.
}
  \item{globalTopLevel}{
    This argument should be filled only when embedding this Tcl-Tk interface in an other. It is the top level of the embedding interface, generally a call to \code{\link{tktoplevel}}.
}
  \item{localTopLevel}{
    This argument should be filled only when embedding this Tcl-Tk interface in an other. It is the local top level to use to build this interface, generally a \code{\link{tkframe}} or \code{\link{ttkframe}}.
}
  \item{scaleFactor}{
    Single numeric value, scaling factor to convert pixel sizes in \code{tkrplot} scale units. \code{NA} will use the OS dependant default values.
}
  \item{hscale}{
    Single numeric value, to be passed to \code{tkrplot}.
}
  \item{vscale}{
    Single numeric value, to be passed to \code{tkrplot}.
}
}

\details{
  Currently two types of files are handled: \code{\link[cghRA:cghRA.regions-class]{cghRA.regions}} objects exported with \code{\link[Rgb]{saveRDT}} and custom tables of segments with an optional header line describing the model.
  
  Custom files are supposed to meet the following criteria:
  \itemize{
    \item{Filename extension must be ".txt".}
    \item{Table separated by tabulations, with dots as decimal separators.}
	\item{Each segment of the genome on a distinct row.}
	\item{A "chrom" column (preferably character) for segment chromosome location.}
	\item{"start" and "end" columns (1 based integers) for position on the chromosome.}
	\item{"probes" (integer) for probe amount in the segment.}
	\item{"logRatio" (numeric) for mean log-ratio of the segment.}
	\item{The first line can hold a model description, as returned by \code{\link{model.test}}. The line must begin with a "#" sign and describe values as "name=value" pairs separated by ", ".}
  }
}

\author{
  Sylvain Mareschal
}

\seealso{
  \code{\link{model.auto}}, \code{\link{model.test}}, \code{\link{tk.cghRA}}
}
