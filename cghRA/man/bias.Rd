\name{bias}
\alias{bias}
\title{
  WACA bias computation for a probe series
}

\description{
  This function computes the various probe-dependant biases used by the Waves aCGH Correction Algorithm (WACA), in order to apply this correction to CGH arrays using these probes.
}

\usage{
  bias(chromFiles, probeChrom, probeStarts, probeEnds,
    chromPattern = "^(.+)\\\\.[^\\\\.]+$", fragSites = c(AluI = "AG|CT", RsaI = "GT|AC"),
    digits = 6, verbose = 1)
}

\arguments{
  \item{chromFiles}{
    Character vector, paths to chromosome sequences (a single fasta file for each chromosome).
}
  \item{probeChrom}{
    \code{character} vector, for each probe its chromosome location.
}
  \item{probeStarts}{
    \code{integer} vector, for each probe its chromosome starting position (first base is 1, starting position is comprised in the probe).
}
  \item{probeEnds}{
    \code{integer} vector, for each probe its chromosome ending position (first base is 1, ending position is comprised in the probe).
}
  \item{chromPattern}{
    Single character value, a regular expression to be used for chromosome name extraction from \code{chromFiles}. It needs to capture a single value for replacement, default value will use the base names of the files without extension as chromosome names.
}
  \item{fragSites}{
    Named character vector describing the restriction enzymes used for the CGH experiment. Restriction sites are described in upper cases, with a pipe at the fragmentation position (see the default value for an example). Only A, C, G and T letters allowed.
}
  \item{digits}{
    Single integer value, to be passed to \code{\link{round}} for all bias values.
}
  \item{verbose}{
    Single integer value, whether to print diagnostic \code{message}s or not.
}
}

\value{
  Returns a double matrix, with probes in rows and the following columns :
    \item{wGC150}{GC frequency in a window of 150 kb on each side of the probe}
    \item{wGC500}{GC frequency in a window of 500 kb on each side of the probe}
    \item{wGCprobe}{GC frequency in the probe sequence}
    \item{wGCfrag}{GC frequency in the fragment holding the probe}
    \item{wFragSize}{Size (in bp) of the fragment holding the probe}
}

\references{
  Lepretre F. et al. (2010) Waved aCGH: to smooth or not to smooth. Nucleic Acids Res. 2010 Apr;38(7):e94
}

\author{
  Sylvain Mareschal
}

\seealso{
  \code{\link{WACA}}, \code{\link{localize}}
}
