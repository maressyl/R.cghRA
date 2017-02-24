\name{model.apply}
\alias{model.apply}
\title{
  Computes copy number for a set of CGH segments
}

\description{
  This function translates log ratios of a set of segments into copy numbers, applying a copy number model as produced by \code{\link{model.auto}} or \code{\link{model.test}}.
  
  If \code{exact} is set set to \code{FALSE}, copy numbers are rounded and consecutive segments with the same copy number are merged.
}

\usage{
  model.apply(segStarts, segEnds, segChroms, segLogRatios, segLengths, model = NA,
    center = model['center'], width = model['width'], ploidy = model['ploidy'],
    exact = FALSE, merge = TRUE)
}

\arguments{
  \item{segStarts}{
    Numeric vector, the starting positions of the CGH segments to modelize.
}
  \item{segEnds}{
    Numeric vector, the endind positions of the CGH segments to modelize.
}
  \item{segChroms}{
    Vector, the chromosome holding the CGH segments to modelize.
}
  \item{segLogRatios}{
    Double vector, the log ratios of the CGH segments to modelize.
}
  \item{segLengths}{
    Numeric vector, the lengths of the CGH segments to modelize.
}
  \item{model}{
    A numeric vector, as returned by \code{\link{model.auto}} or \code{\link{model.test}}. Can be \code{NA} if parameters are provided via other arguments.
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
  \item{exact}{
    Single logical value, whether to return continue copy numbers (double) or discrete ones (integer).
}
  \item{merge}{
    Single logical value, whether to merge consecutive segments with the same copy number when \code{exact} is \code{FALSE}.
}
}

\value{
  Returns a \code{data.frame} describing the segments :
  \item{segStarts}{Extracted from the \code{segStarts} argument.}
  \item{segEnds}{Extracted from the \code{segEnds} argument.}
  \item{segChroms}{Extracted from the \code{segChroms} argument.}
  \item{segLogRatios}{Double, the theoretic log ratio of the segment, with 2 copies as reference.}
  \item{segCopies}{Numeric, the copy number of the segment.}
  \item{segLengths}{Extracted from the \code{segLengths} argument.}
}

\author{
  Sylvain Mareschal
}

\seealso{
  \code{\link{copies}}, \code{\link{model.auto}}, \code{\link{model.test}}
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
  
  # Profile simplification
  segments <- model.apply(
    segStarts,
    segEnds,
    segChroms,
    segLogRatios,
    segLengths,
    model = model,
    exact = FALSE,
    merge = TRUE
  )
  layout(matrix(1:2, ncol=1))
  plot(x=segStarts, y=segLogRatios, type="s", xlab="Position", ylab="Log Ratios")
  plot(x=segments$segStarts, y=segments$segCopies, type="s", xlab="Position", ylab="Copies")
  print(segments)
  
  layout(1)
}
