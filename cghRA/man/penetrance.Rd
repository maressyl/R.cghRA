\name{penetrance}
\alias{penetrance}
\title{
  Penetrance computation from a series of segments
}

\description{
  This function computes the penetrance of various states from a \code{\link{parallelize}}d series of segments.
  
  In each point of the genome, the penetrance is the proportion of the series arrays that show a specific alteration state.
}

\usage{
  penetrance(segParallel, states = list(deletion=c(-Inf, -0.5), gain=c(0.5, Inf)),
    na = c("fill", "keep", "false"), mergeOnValue = FALSE, bool = FALSE, quiet = FALSE)
}

\arguments{
  \item{segParallel}{
    A \code{data.frame}, as returned by \code{\link{parallelize}}.
}
  \item{states}{
    A named \code{list} of numerics defining the boundaries of each state. Each state may be defined by a single value (the only value in \code{segParallel} to link to the state) or by two boundaries (the lower boundary is part of the state, the upper one is not). \code{Inf} and \code{-Inf} can be used as boundaries.
}
  \item{na}{
    Single character value defining how to deal with \code{NA} segments : "fill" fills them when possible (chromosome ends and gaps for which the state is the same on each side), "keep" keeps all of them \code{NA} and "false" always considers them as "not in the state". When \code{NA} remains ("fill" or "keep"), the penetrance frequency is locally computed on non-\code{NA} samples.
}
  \item{mergeOnValue}{
    Single logical value, whether to merge consecutive regions with same penetrance value but distinct alterated sample list.
}
  \item{bool}{
    Single logical value, if \code{TRUE} the penetrance is not returned but logical matrixes of regions 'in state' are returned instead. This is a quite uncommon behavior, allowed essentially for code recycling by other packages, use \code{FALSE}.
}
  \item{quiet}{
    Single logical value, whether to throw diagnosis messages or not.
}
}

\value{
  If \code{bool} is \code{FALSE}, a \code{list} containing a distinct \code{data.frame} for each state, with the following columns :
  \item{chrom}{Character, the chromosomal location of the region described.}
  \item{start}{Integer, the lower coordinate of the region described.}
  \item{end}{Integer, the upper coordinate of the region described.}
  \item{value}{Numeric, the penetrance in the region described for the state described.}
}

\author{
  Sylvain Mareschal
}

\seealso{
  \code{\link{parallelize}}, \code{\link{STEPS}}
}
