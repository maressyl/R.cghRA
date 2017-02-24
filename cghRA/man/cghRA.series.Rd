\name{cghRA.series}
\alias{cghRA.series}
\title{
  cghRA.series class constructor
}

\description{
  This function returns a new \code{\linkS4class{cghRA.series}} object. Elements may be added to the series via the \code{add} method in a second time.
}

\usage{
  cghRA.series(..., .name, warn = TRUE)
}

\arguments{
  \item{\dots}{
    Elements to include in the series, as a single \code{list} or multiple variables containing \code{\linkS4class{cghRA.regions}} objects. Alternatively, a character vector of RDT file paths can be provided.
}
  \item{.name}{
    Single character value, the name of the series.
}
  \item{warn}{
    Single logical value, to be passed to the \code{\linkS4class{cghRA.series}} \code{check} method.
}
}

\value{
  An object of class \code{\linkS4class{cghRA.series}}.
}

\author{
  Sylvain Mareschal
}

\seealso{
  \code{\link{cghRA.series-class}}
}
