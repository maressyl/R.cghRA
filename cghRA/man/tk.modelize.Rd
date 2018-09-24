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
    "Yp", "Yq"), globalTopLevel, localTopLevel, png.res = 100,
    png.file = tempfile(fileext=".png"))
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
  \item{png.res}{
    Single integer value, the resolution of the plot in Pixels Per Inches. Passed to \code{\link{png}}, see the corresponding manual for further details.
}
  \item{png.file}{
    Single character value, the path to the PNG file that is displayed in the main window. The default behavior is to hide it in a temporary location, however you can define this argument to have an easier access to the images displayed in Rgb (the image will be replaced each time Rgb refresh its display).
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
