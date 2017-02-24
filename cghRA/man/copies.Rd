\name{copies}
\alias{LCN}
\alias{copies}
\title{
  LogRatio to copies conversion
}

\description{
  \code{copies} applies a model to a vector of logRatios, converting them into copy amounts.
  
  \code{LCN} is similar, but returns only Log-ratio related Copy Numbers, corresponding to a model with 0 as \code{center}, 1 as \code{width} and 2 as \code{ploidy}. See the references for further details on the models.
}

\usage{
  LCN(x, exact = TRUE)
  copies(x, model = NA, center = model['center'], width = model['width'],
    ploidy = model['ploidy'], exact = TRUE, from = c("logRatios", "LCN", "copies"))
}

\arguments{
  \item{x}{
    Numeric vector, the values to be converted (their nature depends on \code{from}).
}
  \item{model}{
    A numeric vector, as returned by \code{\link{model.auto}} or \code{\link{model.test}}. Can be \code{NA} if parameters are provided via other arguments.
}
  \item{center}{
    Single numeric value, the most common \code{\link{LCN}} within the analyzed genome.
}
  \item{width}{
    Single numeric value, \code{\link{LCN}} gaps between two consecutive real copy amounts.
}
  \item{ploidy}{
    Single numeric value, the real copy amount corresponding to \code{center} \code{\link{LCN}}. A few altered human genome should have a ploidy of 2, use 0 to compute relative copy numbers rather than absolute ones.
}
  \item{exact}{
    Single logical value, whether to \code{\link{round}} copy numbers or not.
}
  \item{from}{
    Single character value defining what computation apply to \code{x}. "logRatios" assumes \code{x} to be logRatios to be converted to copy numbers, applying a full model (\code{center}, \code{width}, \code{ploidy}). "LCN" assumes \code{x} to be Log-ratio related Copy Numbers, as returned by \code{\link{LCN}}, so only the \code{exact} argument is used. "copies" assumes \code{x} to be already modelized copy numbers to be turned back into logRatios, using \code{ploidy} as reference.
}
}

\value{
  A numeric vector the same length as \code{x}.
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
  model <- model.auto(
    segLogRatios = segLogRatios,
    segChroms = segChroms,
    segLengths = segLengths
  )
  
  # Relative copy numbers
  print(
    copies(
      segLogRatios,
      model = model,
      ploidy = 0,
	  exact = FALSE
    )
  )
  
  # Absolute copy number (assuming n=2)
  print(
    copies(
      segLogRatios,
      model = model,
      ploidy = 2,
      exact = FALSE
    )
  )
}
