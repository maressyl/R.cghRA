\name{model.test}
\alias{model.test}
\title{
  Copy number model quality assessment
}

\description{
  This function provides various data to manually fit or upgrade a copy number model, as needed by \code{\link{model.apply}} to translate logRatios into copy numbers.
}

\usage{
  model.test(segLogRatios, segChroms, segLengths = rep(1, length(segLogRatios)),
    model = NA, center = model['center'], width = model['width'],
    ploidy = model['ploidy'],  bw = model['bw'], minDensity = 0.001,
    peakFrom = model['peakFrom'], peakTo = model['peakTo'], graph = TRUE,
    parameters = TRUE, returnPar = FALSE, xlim = c(0, 5), ylim = c(0, max(segLengths)),
    xlab = "Segment copy number", ylab = "Segment length", cex.seg = 0.4, cex.leg = 0.7,
    cex.l2r = 0.7, exclude = c("X", "Y", "Xp", "Xq", "Yp", "Yq"), title = NULL,
    panel = FALSE, klim = NULL, ...)
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
  \item{model}{
    A double vector, as returned by \code{\link{model.auto}} or \code{\link{model.test}}. Can be \code{NA} if parameters are provided via other arguments.
}
  \item{center}{
    Single double value, the center parameter to use in the model.
}
  \item{width}{
    Single double value, the width parameter to use in the model.
}
  \item{ploidy}{
    Single numeric value, copy number supposed to be the most common within the analyzed genome.
}
  \item{bw}{
    Single double value, the bandwidth parameter to use in the model.
}
  \item{minDensity}{
    Single double value, minimal density for a peak to be detected.
}
  \item{peakFrom}{
    Single double value, the peak logRatio lower limit parameter to use in the model.
}
  \item{peakTo}{
    Single double value, the peak logRatio upper limit parameter to use in the model.
}
  \item{graph}{
    Single logical value, whether to plot the density distribution of the segments with the modelized copy numbers or not.
}
  \item{parameters}{
    Single logical value, whether to add a legend to the plot with the parameters and statistics of the model or not.
}
  \item{returnPar}{
    Single logical value, whether to return the \code{par} content (for point identification in interactive plots) or the model statistics.
}
  \item{xlim}{
    Vector of two double values, the boundaries of the plot on the horizontal axis (in \code{\link{LCN}}).
}
  \item{ylim}{
    Vector of two double values, the boundaries of the plot on the vertical axis (in the same units than \code{segLengths}).
}
  \item{xlab}{
    Single character value, the title to print for the horizontal axis.
}
  \item{ylab}{
    Single character value, the title to print for the vertical axis.
}
  \item{cex.seg}{
    Single double value, the character expansion factor for points (segments) on the plot.
}
  \item{cex.leg}{
    Single double value, the character expansion factor for the plot legend.
}
  \item{cex.l2r}{
    Single double value, the character expansion factor for the log-ratio axis of the plot.
}
  \item{exclude}{
    Vector, the chromosomes to exclude from the density computation and to plot with distinct symbols (use \code{NULL} to disable this feature). Sexual chromosomes should be excluded in heterogeneous DNA source, as their desequilibrium (2 'X' and no 'Y' in women) impact normal cells AND tumoral ones.
}
  \item{title}{
    To be passed to \code{legend}, see there for allowed types (usually a single character value).
}
  \item{panel}{
    Single logical value, whether to plot a rotated minimalist graph or a classic one.
}
  \item{klim}{
    Double vector of two values, alternative definition of \code{xlim} in modelized copy numbers rather than LCN.
}
  \item{\dots}{
    Further graphical arguments to be passed to \code{\link{plot}}.
}
}

\value{
  When \code{returnPar} is \code{TRUE}, invisibly returns the \code{par} content, for point identification.
  
  When \code{returnPar} is \code{FALSE}, returns the same vector as \code{\link{model.auto}}, see its help page for further details.
}

\author{
  Sylvain Mareschal
}

\seealso{
  \code{\link{model.auto}}, \code{\link{model.apply}}
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
  autoModel <- model.auto(
    segLogRatios = segLogRatios,
    segChroms = segChroms,
    segLengths = segLengths
  )
  
  layout(matrix(1:2, ncol=1))
  
  # Show automatic model
  model.test(
    segLogRatios = segLogRatios,
    segChroms = segChroms,
    segLengths = segLengths,
	model = autoModel
  )
  
  # Standard model derived from the log ratios definition
  refModel <- model.test(
    segLogRatios = segLogRatios,
    segChroms = segChroms,
    segLengths = segLengths,
    center = 2,
    width = 1,
    bw = 0.1        # Arbitrary
  )
  
  # Differences in scores
  print(autoModel)
  print(refModel)
  
  layout(1)
}
