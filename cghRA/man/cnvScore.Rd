\name{cnvScore}
\alias{cnvScore}
\title{
  Polymorphism likelihood score for a genomic segment
}

\description{
  Computes for each genomic segment provided a score reflecting its likelihood to a polymorphism (CNV) dataset, as can be download from the Database of Genomic Variants.
}

\usage{
  cnvScore(sample.map, dgv.map, hangingThreshold = 0.8, minGain = 0.1, maxPaths = NA,
    trace = FALSE, expand = TRUE, quiet = TRUE)
}

\arguments{
  \item{sample.map}{
    A \code{\linkS4class{segmentMap}} object as returned by \code{\link{map2design}}, corresponding to the mapping of the segments to assess to a common design.
}
  \item{dgv.map}{
    A \code{\linkS4class{segmentMap}} object as returned by \code{\link{map2design}}, corresponding to the mapping of the true polymorphism (CNV) dataset to a common design.
}
  \item{hangingThreshold}{
    Single numeric value, segments to score must cover at least this proportion of union(CNV, segment) for a CNV to be considered. Decrease this value to allow poorly overlapping CNVs to (modestly) contribute to the final score, at the cost of longer computing time.
}
  \item{minGain}{
    Single numeric value, CNVs must add at least this value to the path's score to be retained. Increase this value to allow poorly overlapping CNVs to (modestly) contribute to the final score, at the cost of longer computing time.
}
  \item{maxPaths}{
    Single integer value, the maximal amount of paths to be computed for each segment (use \code{NA} to always compute all of them). Considering that most of the best paths are computed first and final score focus on them, an arbitrary value like 50 can be provided to decrease the computing time with marginal effects on the resulting scores.
}
  \item{trace}{
    Single logical value, whether to produce a trace of every path constructed or only the final CNV score. This is mainly provided for debugging purpose, and increase the computing time. \code{\link{trace2track}} provides graphical means to visualize these traces.
}
  \item{expand}{
    Single logical value, whether to return a vector of scores with one element for each row in \code{sample.map} (\code{FALSE}) or in the original mapped track (\code{TRUE}). As the mapping involves row compression (see \code{\link{map2design}}), producing a vector that can be directly used as a column in the original track needs an expansion step, that can be performed if requested via this argument.
}
  \item{quiet}{
    Single logical value, whether to print diagnostic \code{\link{message}}s or not.
}
}

\value{
  If \code{trace} is \code{FALSE}, returns a numeric vector storing the resulting CNV score. See \code{expand} for further details on this vector size.
  
  If \code{trace} is \code{TRUE}, returns a named \code{list} of two elements: "scores", that holds the numeric vector of scores (see above), and "traces", that described every path that has been built to compute the scores. The columns in "traces" are:
  \item{seg}{Range of the original track indexes corresponding to the assessed segment.}
  \item{seg.score}{Final CNV score for the assessed segment, all paths comprised.}
  \item{path.count}{How many times the CNV path described was built.}
  \item{path.jaccard}{Jaccard index between the assessed segment and the CNV path described.}
  \item{path.cnvCount}{How many CNVs are included in the CNV path described.}
  \item{path.cnvList}{Indexes in \code{dgv.map} of the CNVs retained in the CNV path described.}
  \item{path.score}{'path.jaccard' corrected for the amount of CNVs included in the CNV path described.}
}

\author{
  Sylvain Mareschal
}

\seealso{
  \code{\link{map2design}}, \code{\link{applyMap}}, \code{\link{trace2track}}
}

