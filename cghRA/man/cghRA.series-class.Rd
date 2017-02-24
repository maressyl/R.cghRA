\name{cghRA.series-class}
\Rdversion{1.1}
\docType{class}
\alias{cghRA.series-class}

\title{Class \code{"cghRA.series"}}
\description{
  Objects from this class are collections of \code{\linkS4class{cghRA.regions}} objects, and provide various methods for CGH series analysis.
  
  Objects from this class should always be produced by the \code{\link{cghRA.series}} constructor.
}

\author{
  Sylvain Mareschal
}

\seealso{
  \code{\link{cghRA.series}}, \code{\linkS4class{cghRA.regions}}
}

\keyword{classes}
\section{Extends}{
  All reference classes extend and inherit methods from \code{\linkS4class{envRefClass}}.
}
\section{Fields}{
  \describe{
    \item{\code{arrays}:}{A possibly named \code{list} of \code{\linkS4class{cghRA.regions}} objects.}
    \item{\code{name}:}{Single character value, the custom name of the series.}
  }
}
\section{Methods}{
  \describe{
    \item{\code{add(object)}:}{
      Add an object to the series
    }
    \item{\code{applyMethod(.method, ..., .simplify = , .quiet = )}:}{
      Calls a method on each array of the series\cr
      - \strong{.method}     : single character value, the method to be called.\cr
      - \strong{...}         : arguments to be passed to the method.\cr
      - \strong{.simplify}   : same behavior as sapply() 'simplify' argument.\cr
      - \strong{.quiet}      : single logical value, whether to print iterations or not.
    }
    \item{\code{check(warn = )}:}{
      Raises an error if the object is not valid, else returns TRUE
    }
    \item{\code{get(arrayName)}:}{
      Returns an element from the series
    }
    \item{\code{getArrayNames()}:}{
      Returns a vector of array names
    }
    \item{\code{initialize(name = , arrays = , ...)}:}{}
    \item{\code{last()}:}{
      Refers to the last array added in the series
    }
    \item{\code{LRA(value = , tracks = , ...)}:}{
      Apply the LRA() function to list Long Recurrent Abnormalities (Lenz et al, PNAS 2008).\cr
      - \strong{value}    : single character value, the name of the column to use as copy number estimate ('copies' or 'logRatio').\cr
      - \strong{tracks}   : single logical value, whether to convert output to track.table class or not.
    }
    \item{\code{parallelize(value = , quiet = , tracks = , ...)}:}{
      Apply the parallelize() function to build a summary matrix of the series.\cr
      - \strong{tracks}   : single logical value, whether to convert output to track.table class or not.
    }
    \item{\code{penetrance(tracks = , ...)}:}{
      Apply the penetrance() function to compute the proportion of altered samples for each genomic position.\cr
      - \strong{tracks}   : single logical value, whether to convert output to track.table class or not.
    }
    \item{\code{pool(tracks = , value = , group = , states = , others = , quiet = )}:}{
      Collect and pool all alterated segments from the various samples of the series.\cr
      - \strong{tracks}   : single logical value, whether to convert output to track.table class or not.\cr
      - \strong{value}    : column on which apply a filtering.\cr
      - \strong{group}    : single logical value, whether to visually group segments per samples or not (valid only for tracks=TRUE).\cr
      - \strong{states}   : list of states, see penetrance help page for details. If 'states' is not empty, segments without state will be filtered out.\cr
      - \strong{others}   : character vector, names of other columns to keep.\cr
      - \strong{quiet}    : single logical value, whether to throw diagnosis messages or not.
    }
    \item{\code{SRA(value = , tracks = , ...)}:}{
      Apply the SRA() function to list Short Recurrent Abnormalities (Lenz et al, PNAS 2008).\cr
      - \strong{value}    : single character value, the name of the column to use as copy number estimate ('copies' or 'logRatio').\cr
      - \strong{tracks}   : single logical value, whether to convert output to track.table class or not.
    }
    \item{\code{STEPS(tracks = , ...)}:}{
      Apply the STEPS() function to prioritize commonly altered regions.\cr
      - \strong{tracks}   : single logical value, whether to convert output to track.table class or not.
    }
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
