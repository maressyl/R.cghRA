\name{fillGaps}
\alias{fillGaps}
\title{
  Fill gaps between consecutive segments
}

\description{
  This function enlarges segments on their upper boundary to fill gaps between consecutive segments.
  
  It may be crucial for \code{\link{penetrance}} computation, as they lead to small low steps in penetrance.
}

\usage{
  fillGaps(segTable, isOrdered = FALSE)
}

\arguments{
  \item{segTable}{
    A \code{data.frame} of segments, with at least "chrom" (character), "start" (integer) and "end" (integer) columns.
}
  \item{isOrdered}{
    Single logical value, whether \code{segTable} is already ordered by chromosome and starting position or not.
}
}

\value{
  Returns a \code{data.frame} similar to \code{segTable}.
}

\author{
  Sylvain Mareschal
}
