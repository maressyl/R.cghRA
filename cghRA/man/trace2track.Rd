\name{trace2track}
\alias{trace2track}
\title{
  Converts cnvScore traces to a drawable track
}

\description{
  This function converts the data.frame trace that can be produced by \code{\link{cnvScore}} into a \code{\linkS4class{track.table}} object that can be browsed using Rgb's functions \code{\link[Rgb]{tk.browse}}. and \code{\link[Rgb]{browsePlot}}.
}

\usage{
  trace2track(paths, dgv.map, dgv.track)
}

\arguments{
  \item{paths}{
    A \code{data.frame}, as produced by \code{\link{cnvScore}} with \code{trace=TRUE}.
}
  \item{dgv.map}{
    An integer matrix as returned by \code{\link{map2design}}, corresponding to the mapping of the polymorphism (CNV) dataset to a common design.
}
  \item{dgv.track}{
    A \code{\linkS4class{track.table}}-inheriting object, the original dataset used to produce \code{dgv.map}.
}
}

\value{
  Returns a copy of \code{dgv.track}, in which CNVs are grouped by paths labeled with the resulting score.
}

\author{
  Sylvain Mareschal
}

\seealso{
  \code{\link{cnvScore}}, \code{\link{map2design}}
}

