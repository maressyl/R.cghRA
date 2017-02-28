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
    "Yp", "Yq"), globalTopLevel, localTopLevel, render = c("auto", "png", "tkrplot"),
    tkrplot.scale = 1, png.res = 100, png.file = tempfile(fileext=".png"))
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
  \item{render}{
    Single character value from the ones listed, defining the rendering engine for the plot. "png" is recommended and the default on any platform supporting it (needs Tcl-tk version 8.6 or higher, already available on Linux and MacOS and on Windows with R version 3.4.0 or above), and consists in displaying an export to a PNG file. "tkrplot" is more limited and kept only for backward compatibility, it relies on the external package \code{tkrplot} and the Windows "metafile" format. "auto" (the dzfault) will select the best engine considering to the capabilities of your installation.
}
  \item{tkrplot.scale}{
    Single numeric value, defining a multiplying factor for plot size with the "tkrplot" engine. This argument is mainly provided to temper a bug with the "Font size multiplication factor" feature of last Windows operating system, and get plots filling the whole Tcl-tk window. As an example if you use a 150% factor for fonts and icons in Windows, you should define this argument to 1.5. This has no effect with the "png" engine used on other platforms.
}
  \item{png.res}{
    Single integer value, the resolution of the plot in Pixels Per Inches. Passed to \code{\link{png}}, see the corresponding manual for further details. This has no effect with the "tkrplot" engine used on Windows prior to R version 3.4.0.
}
  \item{png.file}{
    Single character value, the path to the PNG file that is displayed in the main window. The default behavior is to hide it in a temporary location, however you can define this argument to have an easier access to the images displayed in Rgb (the image will be replaced each time Rgb refresh its display). This has no effect with the "tkrplot" engine used on Windows prior to R version 3.4.0.
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
