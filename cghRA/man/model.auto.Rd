\name{model.auto}
\alias{model.auto}
\title{
  Automatic generation of copy number model
}

\description{
  This function computes a copy number model, as needed by \code{\link{model.apply}} to translate logRatios into copy numbers.
}

\usage{
  model.auto(segLogRatios, segChroms, segLengths = rep(1, length(segLogRatios)),
    from = 0.02, to = 0.5, by = 0.001, precision = 512, maxPeaks = 8, minWidth = 0.15,
    maxWidth = 0.9, minDensity = 0.001, peakFrom = -2, peakTo = 1.3, ploidy = 0,
    discreet = FALSE, method = c("stm", "sdd", "ptm"), exclude = c("X", "Y", "Xp", "Xq",
    "Yp", "Yq"))
}

\arguments{
  \item{segLogRatios}{
    Double vector, the log ratios of the CGH segments to modelize.
}
  \item{segChroms}{
    Vector, the chromosome holding the CGH segments to modelize.
}
  \item{segLengths}{
    Double vector, the lengths of the CGH segments to modelize. Amount of probes should be prefered if available, but nucleotide length or no length at all can also be used.
}
  \item{from}{
    Single double value, the minimal bandwidth to test for \code{\link{density}}.
}
  \item{to}{
    Single double value, the maximal bandwidth to test for \code{\link{density}}.
}
  \item{by}{
    Single double value, the precision of the bandwidths to test for \code{\link{density}}.
}
  \item{precision}{
    Single integer value, the amount of points to compute for \code{\link{density}}. As its help page suggests, values greater than 512 should be powers of 2.
}
  \item{maxPeaks}{
    Single integer value, the maximal amount of peaks in the density of distribution to consider a model as valid.
}
  \item{minWidth}{
    Single double value, minimal value allowed for the \code{width} model parameter (thus for tumoral cell proportion in the sample).
}
  \item{maxWidth}{
    Single double value, maximal value allowed for the \code{width} model parameter (thus for tumoral cell proportion in the sample).
}
  \item{minDensity}{
    Single double value, minimal density for a peak to be detected.
}
  \item{peakFrom}{
    Single double value, minimal logRatio for a peak to be detected. Use \code{NA} for no lower limit. Only 1, 2 and 3 copies peaks should be considered for a more precise model.
}
  \item{peakTo}{
    Single double value, maximal logRatio for a peak to be detected. Use \code{NA} for no upper limit. Only 1, 2 and 3 copies peaks should be considered for a more precise model.
}
  \item{ploidy}{
    Single numeric value, copy number supposed to be the most common within the analyzed genome.
}
  \item{discreet}{
    Single logical value, if \code{FALSE} a fail in modelization raises an error, if \code{TRUE} it returns a \code{NA} filled model.
}
  \item{method}{
    Single character value, the statistic to minimize ("stm" is default). See below for further details.
}
  \item{exclude}{
    Vector, the chromosomes to exclude from the density computation and to plot with distinct symbols (use \code{NULL} to disable this feature). Sexual chromosomes should be excluded in heterogeneous DNA source, as their desequilibrium (2 'X' and no 'Y' in women) impact normal cells AND tumoral ones.
}
}

\details{
  More details about the cghRA copy number model and modelization can be found in the vignette associated with this package, as well as in the related publication. Once the parameters of a model (width and center) are set, three scores can be computed to assess its fitness to the data :
  
  \strong{STM} is the "Segment To Model" score, computed at the segment level as the average of the residuals weighted by the segment size (in probe counts). Residuals are computed as the absolute difference between exact copy numbers (see the \code{\link{copies}} function) and their rounding, assuming that copy numbers should be integers and that decimal parts are noise in the model. This is the recommended score to use with cghRA.
  
  \strong{PTM} is the "Peak To Model" score, computed at the peak level as the average of the residuals. Residuals are computed as the absolute difference between exact copy numbers (see the \code{\link{copies}} function) and their rounding, assuming that copy numbers should be integers and that decimal parts are noise in the model.
  
  \strong{SDD} is the "Standard Deviation of peak Differences" score. As its name suggests, it is computed as the \code{sd} or differences between consecutive peaks, considering that good models should show very regularly spaced density peaks.  
}

\value{
  Returns a double vector, with the following values :
  \item{bw}{Bandwidth used for \code{\link{density}} computation.}
  \item{peaks}{Amount of peaks considered in the model.}
  \item{peakFrom}{See the \code{peakFrom} argument.}
  \item{peakTo}{See the \code{peakTo} argument.}
  \item{center}{Center parameter of the model.}
  \item{width}{Width paremeter of the model.}
  \item{ploidy}{Ploidy paremeter of the model, as provided.}
  \item{sdd}{Quality statistic, see 'Details'.}
  \item{ptm}{Quality statistic, see 'Details'.}
  \item{stm}{Quality statistic, see 'Details'.}
}

\author{
  Sylvain Mareschal
}

\seealso{
  \code{\link{model.test}}, \code{\link{model.apply}}
}

\examples{
  # Generating random segmentation results
  ## with 30% normal cells contamination
  ## with +10% for normal DNA labelling
  segLogRatios <- c(
    rnorm(
      sample(5:20, 1),
      mean = log((1*0.7 + 2*0.3)/(2*1.1), 2),   # One deletion
      sd = 0.08
    ),
    rnorm(
      sample(80:120, 1),
      mean = log(2/(2*1.1), 2),                 # No alteration
      sd = 0.08
    ),
    rnorm(
      sample(40:60, 1),
      mean = log((3*0.7 + 2*0.3)/(2*1.1), 2),   # One more copy
      sd = 0.08
    )
  )
  segLogRatios <- sample(segLogRatios)
  segLengths <- as.integer(3 + round(rchisq(length(segLogRatios), 1)*100))
  segEnds <- cumsum(segLengths)
  segStarts <- c(1L, head(segEnds, -1))
  segChroms <- rep("chr1", length(segEnds))
  
  # Generated genome
  genome <- data.frame(
    segChroms,
    segStarts,
    segEnds,
    segLogRatios,
    segLengths
  )
  print(genome)
  
  # Automatic modelization
  model <- model.auto(
    segLogRatios = segLogRatios,
    segChroms = segChroms,
    segLengths = segLengths
  )
  print(model)
}
