\name{parallelize}
\alias{parallelize}
\title{
  Reshapes a list of segments
}

\description{
  This function reshapes a list of segment \code{data.frame}s (with chromosomal location and value) into a single \code{data.frame} containing a column for each element of the list (typically samples) and a the minimal amount of regions in rows.
}

\usage{
  parallelize(segTables, value = "logRatio", digits = 3, quiet = FALSE, chroms = NULL)
}

\arguments{
  \item{segTables}{
    An eventually named \code{list} of \code{data.frame}s to reshape. All the \code{data.frame}s must contain at least "chrom" (character), "start" (integer), "end" (integer) columns, and the column defined by \code{value}.
	
	Can also be a single \code{data.frame} containing all the segments, with a \code{.sampleIdentity} integer column.
}
  \item{value}{
    Single character value, the column name from which extract values that will fill the output cells.
}
  \item{digits}{
    Single integer value to be passed to \code{\link{round}} for each cell of the output (\code{NA} disables the rounding step).
}
  \item{quiet}{
    Single logical value, whether to throw diagnosis messages or not.
}
  \item{chroms}{
    Character vector, the names of chromosomes to restrain the analysis on (frequently autosomes). If \code{NULL}, all chromosomes in \code{segTable} will be used.
}
}

\value{
  Returns a \code{data.frame} with the following columns :
  \item{chrom}{Character, the chromosomal location of the region described.}
  \item{start}{Integer, the lower coordinate of the region described.}
  \item{end}{Integer, the upper coordinate of the region described.}
  \item{...}{For each element of \code{segTables} a column with the value extracted from the \code{value} column of the according \code{data.frame}.}
}

\author{
  Sylvain Mareschal
}

\seealso{
  \code{\link{penetrance}}
}
