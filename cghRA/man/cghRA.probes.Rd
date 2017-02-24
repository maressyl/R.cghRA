\name{cghRA.probes-constructor}
\alias{cghRA.probes}
\title{
  cghRA.probes class constructor
}

\description{
  This function returns a new \code{\linkS4class{cghRA.probes}} object from various arguments.
  
  Notice the \code{new()} alternative can be used to produce an empty object, setting only the fields not the content.
}

\usage{
  cghRA.probes(..., .name, warn = TRUE)
}

\arguments{
  \item{\dots}{
    Arguments to be passed through the inherited constructors up to \code{\link{refTable}}.
  }
  \item{.name}{
    Single character value, a custom name for the object (for printing purpose essentially).
  }
  \item{warn}{
    Single logical value, to be passed to the \code{\linkS4class{cghRA.probes}} \code{check} method.
  }
}

\value{
  An object of class \code{\linkS4class{cghRA.probes}}.
}

\author{
  Sylvain Mareschal
}

\seealso{
  \code{\link{cghRA.probes-class}}, \code{refTable-class}
  
  \code{\link{Agilent.probes}}
}
