\name{SRA & LRA}
\alias{xRA}
\alias{SRA}
\alias{LRA}
\title{
  Short/Long Recurrent Abnormalities detection
}

\description{
  These functions extract Short Reccurent Abnormalities (SRA) and Long Reccurent Abnormalities (LRA) from a CGH array series, as described by Lenz et al. (2008).
  
  The processing core \code{xRA} is common for both analysis, but is not intended to be called directly. Use the \code{SRA} and \code{LRA} wrappers instead.
}

\usage{
  xRA(segTables, value = "copies", states = list(deletion=c(-Inf,-0.5), gain=c(0.5,Inf)),
    sampleMin = 2, quiet = FALSE, lengthMax, lengthMin, gaps.width, gaps.ratio)
  SRA(...)
  LRA(...)
}

\arguments{
  \item{segTables}{
    An eventually named \code{list} of \code{data.frame}s to reshape. All the \code{data.frame}s must contain at least "chrom" (character), "start" (integer), "end" (integer) columns, and the column defined by \code{value}.
	
	Can also be a single \code{data.frame} containing all the segments, with a \code{.sampleIdentity} integer column.
}
  \item{value}{
    Single character value, the column name from which extract values that will fill the output cells.
}
  \item{states}{
    A named \code{list} of numerics defining the boundaries of each state. Each state may be defined by a single value (the only value in \code{segParallel} to link to the state) or by two boundaries (the lower boundary is part of the state, the upper one is not). \code{Inf} and \code{-Inf} can be used as boundaries.
}
  \item{sampleMin}{
    Single numeric value, minimal amount of samples in the 'overlapping group'. If lesser than 1, interpreted as a proportion of the sample count. Large values decrease processing time and SRA amounts.
}
  \item{quiet}{
    Single logical value, whether to print diagnostic \code{\link{message}}s or not.
}
  \item{lengthMax}{
    Single integer value, segments larger than this value will be filtered out (25 Mb for SRA, \code{NA} for LRA). Use \code{NA} to disabled length filtering.
}
  \item{lengthMin}{
    Single integer value, segments shorter than this value will be filtered out (\code{NA} for SRA, 15 Mb for LRA). Use \code{NA} to disabled length filtering.
}
  \item{gaps.width}{
    Single integer value, alterated segments separated by a gap shorter than this value will be merged (see also 'gaps.ratio'; 500 kb for SRA, 10 Mb for LRA). Use \code{NA} to disabled gap filling.
}
  \item{gaps.ratio}{
    Single numeric value, for a gap to be filled its two neighbors must be this value larger than it (see also 'gaps.width'; 1 for SRA, 1.5 for LRA). Use \code{NA} to disabled gap filling.
}
  \item{\dots}{
    The \code{SRA} and \code{LRA} functions are only wrappers to \code{xRA} with distinct \code{lengthMax}, \code{lengthMin}, \code{gaps.width} and \code{gaps.ratio} values, all other arguments are passed through to \code{xRA}.
}
}

\value{
  Returns a list with a \code{data.frame} for each state :
  \item{chrom}{Chromosomal location.}
  \item{inPeak}{Numeric, proportion of the sample series in the 'overlapping group'.}
  \item{overlap.start, overlap.end}{Integer, position on the chromosome for the highest peak of the SRA (region covered by the whole 'overlapping group').}
  \item{start, end}{Integer, position on the chromosome for the SRA itself (largest region covered by 2/3 of the 'overlapping group').}
  \item{extended.start, extended.end}{Integer, position on the chromosome for the extended SRA (largest region covered by 1/3 of the 'overlapping group').}
}

\references{
  Lenz G et al. "Molecular subtypes of diffuse large B-cell lymphoma arise by distinct genetic pathways". Proc Natl Acad Sci U S A. 2008 Sep 9;105(36):13520-5 (Supporting Information)
}

\author{
  Sylvain Mareschal
}

\note{
  For Long Reccurent Abnormalities, Lenz et al. suggest to filter out regions involved in abnormal chromosome arms. For technical reasons, this filter was \strong{NOT} implemented.
}

\seealso{
  \code{\link{STEPS}}
}
