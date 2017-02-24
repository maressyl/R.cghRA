\name{map2design}
\alias{map2design}
\alias{applyMap}
\title{
  Update a track coordinates to match a distinct CGH design
}

\description{
  Remapping a \code{\linkS4class{track.table}} object storing genomic segments to a specific CGH design consists of two steps :
  \itemize{
    \item{The production of a map, which defines the coordinates of each segment by the indexes of the first and last CGH probes included in it (\code{\link{map2design}}).}
    \item{The update of the genomic coordinates of the original track, using the map and the design (\code{\link{applyMap}}).}
  }
}

\usage{
  map2design(track, design, minProbes = 1, quiet = FALSE, warn = TRUE)
  applyMap(track, map, design)
}

\arguments{
  \item{track}{
    A \code{\linkS4class{track.table}}-inheriting object, storing one row for each genomic segment of interest in a CGH-like experiment.
}
  \item{design}{
    A \code{\linkS4class{track.table}}-inheriting object (preferably a \code{\linkS4class{cghRA.design}} object), storing one row for each probe in the design data is to be remapped on.
}
  \item{minProbes}{
    Single integer value, the amount of probes a segment in \code{track} must cover to be retained.
}
  \item{quiet}{
    Single logical value, whether to print diagnostic \code{\link{message}}s or not.
}
  \item{map}{
    An integer matrix defining the mapping of \code{track} to \code{design}, as produced by \code{\link{map2design}}.
}
  \item{warn}{
    Single logical value, to be passed to the \code{check} method of the newly created \code{\linkS4class{segmentMap}} object.
}
}

\value{
  \code{map2design} returns an integer matrix with 3 columns and row names. Columns "start" and "end" define the coordinates of a segment as probe indexes in \code{design}, and column "count" allow to group segments with the same remapped coordinates. Row names correspond to the index range of the corresponding segments in the original \code{track}.
  
  \code{applyMap} returns a copy of \code{track}, in which \code{start} and \code{end} coordinates have been updated to match coordinates of probes in \code{design}. Segments that do not overlap at least \code{minProbes} probe in \code{design} are excluded.
}

\author{
  Sylvain Mareschal
}

\seealso{
  \code{\link{cnvScore}}
}

