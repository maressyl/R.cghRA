\name{STEPS}
\alias{STEPS}
\title{
  Selective Trends Evidenced by Penetrance Surge
}

\description{
  This function identifies and prioritize Selective Trends Evidenced by Penetrance Surge in a CGH array series. STEPS is an alternative to the Minimal Common Region (MCR) algorithms, with the aim to identify regions frequently amplified or deleted.
}

\usage{
  STEPS(segPenetrance, dpen = 2, vpen = 0.8, gpen = 0.3, threshold = NA,
    nested = c("merge", "flag", "none"), digits = 3, chromEnd = FALSE, quiet = FALSE)
}

\arguments{
  \item{segPenetrance}{
    A \code{data.frame}, as a single element from the \code{list} returned by the \code{\link{penetrance}} function.
}
  \item{dpen}{
    Single numeric value, penalty to apply to penetrance increases.
}
  \item{vpen}{
    Single numeric value, penalty to apply to penetrance differences between wide boundaries.
}
  \item{gpen}{
    Single numeric value, penalty to apply to genomic assymetry.
}
  \item{threshold}{
    Single numeric value, minimum STEPS score to filter results. 0 is the less stringent threshold to use, as negative scores correspond to assymetric STEPS (ascending only on a side). Higher values will return less results (focusing on the most significant ones), however scoring and boundaries of the results will not be impacted.
}
  \item{nested}{
    Single character value, defining how to deal with overlapping STEPS. "merge" will only keep for each set of overlapping STEPS the one with the highest score, "flag" will preserve all the STEPS but add a "nest" column with a distinct ID for each nest, and "none" won't do anything about this.
}
  \item{digits}{
    Single integer value, to be passed to \code{\link{round}} for score computations.
}
  \item{chromEnd}{
    Single logical value, whether to consider chromosome ending as a penetrance drop or not.
}
  \item{quiet}{
    Single logical value, whether to throw diagnosis messages or not.
}
}

\details{
  When a specific gene alteration induces a cell selection (like in tumors), it leads to different altered fragments from a patient to an other. All these fragments have a region in common : the region containing the selecting gene (the Minimal Common Region). Such patterns can be extracted from the penetrance, as they lead to 'stairway' patterns in specific locations.

  This function crawls along the penetrance from every available starting point, computing in both directions a score : a descending step grants the penetrance difference (in percents) while an ascending step penalizes by the penetrance difference multiplied by \code{penalty}. In each direction, the maximal score is used as boundary, and a total STEPS score for the starting point is computed as \code{2} \code{*} \code{(leftMax} \code{+} \code{rightMax)} \code{-} \code{abs(leftMax} \code{-} \code{rightMax)}.
  
  The greatest scores highlight symetric STEPS with high descending paths on both sides.
}

\value{
  Returns a subset of \code{segPenetrance} with the following additionnal columns :
  \item{score}{Numeric, the two-side score for the described starting point (see 'Details').}
  \item{leftBoundary}{Integer, position considered as the left boundary of the stairway pattern.}
  \item{leftScore}{Numeric, score for the left side of the STEPS (see 'Details').}
  \item{rightBoundary}{Integer, position considered as the right boundary of the stairway pattern.}
  \item{rightScore}{Numeric, score for the right side of the STEPS (see 'Details').}
}

\author{
  Sylvain Mareschal
}

\seealso{
  \code{\link{penetrance}}, \code{\link{SRA}}
}
