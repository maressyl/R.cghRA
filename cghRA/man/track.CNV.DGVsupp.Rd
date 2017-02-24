\name{track.CNV.DGVsupp}
\alias{track.CNV.DGVsupp}
\title{
  DGV supporting variant parser
}

\description{
  This function constructs \code{\linkS4class{track.CNV}} objects from free annotation files provided by the Database of Genomic Variants.
  
  It is designed to parse \strong{supporting variants}, as opposed to \code{\link[Rgb]{track.CNV.DGV}} provided by Rgb which is designed to parse \strong{DGV Variants}.
}

\usage{
  track.CNV.DGVsupp(file, name = "DGV CNV (supporting variants)", quiet = FALSE, ...)
}

\arguments{
  \item{file}{
    Single character value, the path to the raw file to parse. See the 'References' section below.
}
  \item{name}{
    Single character value, the \code{name} field for the \code{\linkS4class{track.table}} object.
}
  \item{quiet}{
    Single logical value, whether to print diagnostic \code{\link{message}}s or not.
}
  \item{\dots}{
    Further arguments are passed to the class constructor, as a result most of the handled arguments are \code{\link{track.table}} arguments. Consider notably \code{.organism} and \code{.assembly} for track annotation.
}
}

\value{
  Return a \code{\linkS4class{track.CNV}} object.
}

\references{
  Example of raw file (human assembly 'hg19') : \url{http://dgv.tcag.ca/dgv/docs/GRCh37_hg19_supportingvariants_2014-10-16.txt}
}

\author{
  Sylvain Mareschal
}

\seealso{
  \code{\link{track.table-class}}, \code{\link{track.CNV-class}}, \code{\link{track.CNV.DGV}}
}

