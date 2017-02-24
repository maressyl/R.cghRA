\name{cghRA.array}
\alias{cghRA.array}
\title{
  cghRA.array class constructor
}

\description{
  This function returns a new \code{\linkS4class{cghRA.array}} object from various arguments.
}

\usage{
  cghRA.array(.design, .probes, .name, .parameters, warn = TRUE)
}

\arguments{
  \item{.design}{
    An object of class \code{\linkS4class{cghRA.design}}, as returned by the \code{\link{cghRA.design}} constructor.
  }
  \item{.probes}{
    An object of class \code{\linkS4class{cghRA.probes}}, as returned by the \code{\link{cghRA.probes}} constructor.
  }
  \item{.name}{
    Single character value, to fill the \code{name} field inherited from \code{\linkS4class{drawable}}.
  }
  \item{.parameters}{
    A \code{list} of drawing parameters, to fill the \code{parameters} field of the object.
  }
  \item{warn}{
    Single logical value, to be passed to the \code{\link{cghRA.array-class}} \code{check} method.
  }
}

\value{
  An object of class \code{\linkS4class{cghRA.array}}.
}

\author{
  Sylvain Mareschal
}

\seealso{
  \code{\link{cghRA.array-class}}
}
