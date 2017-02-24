\name{cghRA.design-constructor}
\alias{cghRA.design}
\title{
  cghRA.design class constructor
}

\description{
  This function returns a new \code{\linkS4class{cghRA.design}} object from various arguments.
  
  Notice the \code{new()} alternative can be used to produce an empty object, setting only the fields not the content.
}

\usage{
  cghRA.design(..., warn = TRUE)
}

\arguments{
  \item{\dots}{
    Arguments to be passed through the inherited constructors up to \code{\link{refTable}}.
}
  \item{warn}{
    Single logical value, to be passed to the \code{\linkS4class{cghRA.design}} \code{check} method.
}
}

\value{
  An object of class \code{\linkS4class{cghRA.design}}.
}

\author{
  Sylvain Mareschal
}

\seealso{
  \code{\link{cghRA.design-class}}, \code{track.table-class}, \code{refTable-class}
  
  \code{\link{Agilent.design}}
}
