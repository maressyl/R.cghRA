\name{cghRA.regions-constructor}
\alias{cghRA.regions}
\title{
  cghRA.regions class constructor
}

\description{
  This function returns a new \code{\linkS4class{cghRA.regions}} object from various arguments.
  
  Notice the \code{new()} alternative can be used to produce an empty object, setting only the fields not the content.
}

\usage{
  cghRA.regions(..., .model, warn = TRUE)
}

\arguments{
  \item{\dots}{
    Arguments to be passed through the inherited constructors up to \code{\link{refTable}}.
}
  \item{.model}{
    Numeric vector, to fill the \code{model} field of the object.
}
  \item{warn}{
    Single logical value, to be passed to the \code{\linkS4class{cghRA.regions}} \code{check} method.
}
}

\value{
  An object of class \code{\linkS4class{cghRA.regions}}.
}

\author{
  Sylvain Mareschal
}

\seealso{
  \code{\link{cghRA.regions-class}}, \code{track.table-class}, \code{refTable-class}
}
