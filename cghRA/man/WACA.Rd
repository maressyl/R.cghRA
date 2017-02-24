\name{WACA}
\alias{WACA}
\title{
  Waves aCGH Correction Algorithm
}

\description{
  This function applies the Waves aCGH Correction Algorithm to a series a logRatio (usually a complete series of probe logRatio from a single CGH array), using the probe-dependant biases computed by the \code{\link{bias}} function.
}

\usage{
  WACA(probeNames, probeLogRatios, bias, forceBiasOrdering = TRUE)
}

\arguments{
  \item{probeNames}{
    Character vector, the names of the probes to correct. All these names should be present in \code{bias} row.names.
}
  \item{probeLogRatios}{
    Numeric vector, the logRatios of the probes to correct.
}
  \item{bias}{
    A \code{data.frame}, as returned by the \code{\link{bias}} function.
}
  \item{forceBiasOrdering}{
    Single logical value, whether to force the \code{bias} data.frame ordering / subsetting / replication or not. \code{bias} must be ordered according to \code{probeNames} (that can contain duplicates), if they are not the former needs to be reordered. If they have different lengths, ordering is forced. If not, it is up to the user to assure they are or to set \code{forceBiasOrdering} to \code{TRUE} (the default value). It might be time-saving to order \code{bias} manually and set this parameter to \code{FALSE} when applying WACA on several arrays from the same design.
}
}

\value{
  Returns a numeric vector with the corrected logRatios, preserving the \code{probeNames} and \code{probeLogRatios} order.
}

\references{
  Lepretre F. et al. (2010) Waved aCGH: to smooth or not to smooth. Nucleic Acids Res. 2010 Apr;38(7):e94
}

\author{
  Sylvain Mareschal
}

\seealso{
  \code{\link{bias}}
}
