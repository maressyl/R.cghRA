\name{cghRA.probes-class}
\Rdversion{1.1}
\docType{class}
\alias{cghRA.probes-class}

\title{Class \code{"cghRA.probes"}}
\description{
  This class is part of the \code{\linkS4class{cghRA.array}} class, designed to store all probe-related values of a single CGH array.
  
  Objects from this class can be produced by the \code{\link{cghRA.array}} constructor or by the \code{\link{process}} function, its interfaces \code{\link{tk.process}} and \code{\link{tk.cghRA}} or their sub-functions.
}

\author{
  Sylvain Mareschal
}

\seealso{
  \code{\link{cghRA.array-class}}, \code{\link{refTable-class}}, \code{\link{tk.process}}
}

\keyword{classes}
\section{Extends}{
  Class \code{\linkS4class{refTable}}, directly.

  All reference classes extend and inherit methods from \code{\linkS4class{envRefClass}}.
}
\section{Fields}{
  \describe{
    \item{\code{name}:}{Custom name for the object, as a \code{character} vector of length 1.}
  }
  The following fields are inherited (from the corresponding class):
  \itemize{
    \item colCount (\linkS4class{refTable})
    \item colIterator (\linkS4class{refTable})
    \item colNames (\linkS4class{refTable})
    \item colReferences (\linkS4class{refTable})
    \item rowCount (\linkS4class{refTable})
    \item rowNamed (\linkS4class{refTable})
    \item rowNames (\linkS4class{refTable})
    \item values (\linkS4class{refTable})
  }
}
\section{Methods}{
  The following methods are inherited (from the corresponding class):
  \itemize{
    \item addColumn (\linkS4class{refTable})
    \item addDataFrame (\linkS4class{refTable})
    \item addEmptyRows (\linkS4class{refTable})
    \item addList (\linkS4class{refTable})
    \item addVectors (\linkS4class{refTable})
    \item callSuper (\linkS4class{envRefClass})
    \item check (\linkS4class{refTable}, overloaded)
    \item coerce (\linkS4class{refTable})
    \item colOrder (\linkS4class{refTable})
    \item copy (\linkS4class{refTable})
    \item delColumns (\linkS4class{refTable})
    \item erase (\linkS4class{refTable})
    \item export (\linkS4class{envRefClass})
    \item extract (\linkS4class{refTable})
    \item field (\linkS4class{envRefClass})
    \item fill (\linkS4class{refTable})
    \item getClass (\linkS4class{envRefClass})
    \item getColCount (\linkS4class{refTable})
    \item getColNames (\linkS4class{refTable})
    \item getLevels (\linkS4class{refTable})
    \item getRefClass (\linkS4class{envRefClass})
    \item getRowCount (\linkS4class{refTable})
    \item getRowNames (\linkS4class{refTable})
    \item import (\linkS4class{envRefClass})
    \item indexes (\linkS4class{refTable})
    \item initFields (\linkS4class{envRefClass})
    \item initialize (\linkS4class{refTable}, overloaded)
    \item metaFields (\linkS4class{refTable})
    \item rowOrder (\linkS4class{refTable})
    \item setColNames (\linkS4class{refTable})
    \item setLevels (\linkS4class{refTable})
    \item setRowNames (\linkS4class{refTable})
    \item show (\linkS4class{refTable}, overloaded)
    \item trace (\linkS4class{envRefClass})
    \item types (\linkS4class{refTable})
    \item untrace (\linkS4class{envRefClass})
    \item usingMethods (\linkS4class{envRefClass})
  }
}
