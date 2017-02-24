\name{segmentMap-class}
\Rdversion{1.1}
\docType{class}
\alias{segmentMap-class}

\title{Class \code{"segmentMap"}}
\description{
  Efficient storage of a large collection of genomic intervals, located using probe IDs from a specific array design rather than genomic coordinates. Objects of this class are essentially intended to be produced by the \code{\link{map2design}} function, and used by the \code{\link{cnvScore}} function.
}

\author{
  Sylvain Mareschal
}

\seealso{
  \code{\link{map2design}}, \code{\link{cnvScore}}
}

\keyword{classes}
\section{Extends}{
  All reference classes extend and inherit methods from \code{\linkS4class{envRefClass}}.
}
\section{Fields}{
  \describe{
    \item{\code{designName}:}{Single character value, the content of the \code{name} field of the \code{\linkS4class{cghRA.design}} object used to produce the object.}
    \item{\code{designSize}:}{Single integer value, the row count in the \code{\linkS4class{cghRA.design}} object used to produce the object.}
    \item{\code{map}:}{Integer matrix with one row for each distinct genomic interval in the mapped \code{\linkS4class{track.table}} object. The columns are \code{start} and \code{end}, the indexes of the first and last design elements in the interval and \code{count}, the amount of such intervals in the mapped object. Row names of this matrix list the indexes of the corresponding mapped object intervals.}
    \item{\code{trackName}:}{Single character value, the content of the \code{name} field of the mapped \code{\linkS4class{track.table}} object.}
    \item{\code{trackSize}:}{Single integer value, the row count in the mapped \code{\linkS4class{track.table}} object.}
  }
}
\section{Methods}{
  \describe{
    \item{\code{check(warn = )}:}{
      Raises an error if the object is not valid, else returns TRUE
    }
    \item{\code{initialize(map = , trackName = , trackSize = , designName = , designSize = , ...)}:}{}
  }
  The following methods are inherited (from the corresponding class):
  \itemize{
    \item callSuper (\linkS4class{envRefClass})
    \item copy (\linkS4class{envRefClass})
    \item export (\linkS4class{envRefClass})
    \item field (\linkS4class{envRefClass})
    \item getClass (\linkS4class{envRefClass})
    \item getRefClass (\linkS4class{envRefClass})
    \item import (\linkS4class{envRefClass})
    \item initFields (\linkS4class{envRefClass})
    \item show (\linkS4class{envRefClass}, overloaded)
    \item trace (\linkS4class{envRefClass})
    \item untrace (\linkS4class{envRefClass})
    \item usingMethods (\linkS4class{envRefClass})
  }
}
