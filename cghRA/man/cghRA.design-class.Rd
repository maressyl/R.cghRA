\name{cghRA.design-class}
\Rdversion{1.1}
\docType{class}
\alias{cghRA.design-class}

\title{Class \code{"cghRA.design"}}
\description{
    This class is part of the \code{\linkS4class{cghRA.array}} class. A single object from this class is used to store informations about probes for series of arrays sharing the same CGH design, in order to store only array-specific values in the array variables.
	
	Objects from this class can be produced by the \code{\link{cghRA.design}}, \code{\link{Agilent.design}} and \code{\link{custom.design}} constructors. Alternatively they can be produced by the interactive function \code{\link{tk.design}}, included in \code{\link{tk.cghRA}}.
}

\author{
  Sylvain Mareschal
}

\seealso{
  \code{\link{cghRA.design}}, \code{\link{Agilent.design}}, \code{\link{custom.design}}, \code{\link{tk.design}}
  
  \code{\link{cghRA.array-class}}, \code{\link{refTable-class}}
}

\keyword{classes}
\section{Extends}{
  Class \code{\linkS4class{track.table}}, directly.\cr
  Class \code{\linkS4class{refTable}}, by class \code{\linkS4class{track.table}}, distance 2.\cr
  Class \code{\linkS4class{crossable}}, by class \code{\linkS4class{track.table}}, distance 2.\cr
  Class \code{\linkS4class{sliceable}}, by class \code{\linkS4class{track.table}}, distance 3.\cr
  Class \code{\linkS4class{drawable}}, by class \code{\linkS4class{track.table}}, distance 4.

  All reference classes extend and inherit methods from \code{\linkS4class{envRefClass}}.
}
\section{Fields}{
  The following fields are inherited (from the corresponding class):
  \itemize{
    \item assembly (\linkS4class{track.table})
    \item checktrack (\linkS4class{track.table})
    \item colCount (\linkS4class{refTable})
    \item colIterator (\linkS4class{refTable})
    \item colNames (\linkS4class{refTable})
    \item colReferences (\linkS4class{refTable})
    \item index (\linkS4class{track.table})
    \item name (\linkS4class{drawable})
    \item organism (\linkS4class{track.table})
    \item parameters (\linkS4class{drawable})
    \item rowCount (\linkS4class{refTable})
    \item rowNamed (\linkS4class{refTable})
    \item rowNames (\linkS4class{refTable})
    \item sizetrack (\linkS4class{track.table})
    \item subtrack (\linkS4class{track.table})
    \item values (\linkS4class{refTable})
  }
}
\section{Methods}{
  \describe{
    \item{\code{bias(...)}:}{
      Computes the Waves aCGH Correction Algorithm (Lepretre et al. 2009) bias for the current design.\cr
      - \strong{...}   : arguments to be passed to the bias() function (except from 'probeChrom', 'probeStarts' and 'probeEnds').
    }
    \item{\code{remap(...)}:}{
      Recomputes the coordinates of the probes from the probes and genome sequences. Forces 'chrom' to factor, keeping levels if available.\cr
      - \strong{...}   : arguments to be passed to the localize() function.
    }
  }
  The following methods are inherited (from the corresponding class):
  \itemize{
    \item addArms (\linkS4class{track.table})
    \item addColumn (\linkS4class{track.table})
    \item addDataFrame (\linkS4class{refTable})
    \item addEmptyRows (\linkS4class{refTable})
    \item addList (\linkS4class{track.table})
    \item addVectors (\linkS4class{refTable})
    \item buildCalls (\linkS4class{track.table})
    \item buildGroupPosition (\linkS4class{track.table})
    \item buildGroupSize (\linkS4class{track.table})
    \item buildIndex (\linkS4class{track.table})
    \item callParams (\linkS4class{drawable})
    \item callSuper (\linkS4class{envRefClass})
    \item check (\linkS4class{track.table}, overloaded)
    \item chromosomes (\linkS4class{track.table})
    \item coerce (\linkS4class{track.table})
    \item colOrder (\linkS4class{refTable})
    \item copy (\linkS4class{refTable})
    \item cross (\linkS4class{crossable})
    \item defaultParams (\linkS4class{track.table}, overloaded)
    \item delColumns (\linkS4class{track.table})
    \item draw (\linkS4class{sliceable})
    \item erase (\linkS4class{refTable})
    \item eraseArms (\linkS4class{track.table})
    \item export (\linkS4class{envRefClass})
    \item extract (\linkS4class{refTable})
    \item field (\linkS4class{envRefClass})
    \item fill (\linkS4class{track.table})
    \item fix.param (\linkS4class{drawable})
    \item getChromEnd (\linkS4class{track.table})
    \item getClass (\linkS4class{envRefClass})
    \item getColCount (\linkS4class{refTable})
    \item getColNames (\linkS4class{refTable})
    \item getLevels (\linkS4class{refTable})
    \item getName (\linkS4class{drawable})
    \item getParam (\linkS4class{drawable})
    \item getRefClass (\linkS4class{envRefClass})
    \item getRowCount (\linkS4class{refTable})
    \item getRowNames (\linkS4class{refTable})
    \item import (\linkS4class{envRefClass})
    \item indexes (\linkS4class{refTable})
    \item initFields (\linkS4class{envRefClass})
    \item initialize (\linkS4class{track.table}, overloaded)
    \item isArmed (\linkS4class{track.table})
    \item metaFields (\linkS4class{track.table})
    \item rowOrder (\linkS4class{track.table})
    \item segMerge (\linkS4class{track.table})
    \item segOverlap (\linkS4class{track.table})
    \item setColNames (\linkS4class{track.table})
    \item setLevels (\linkS4class{track.table})
    \item setName (\linkS4class{drawable})
    \item setParam (\linkS4class{drawable})
    \item setRowNames (\linkS4class{refTable})
    \item show (\linkS4class{track.table}, overloaded)
    \item size (\linkS4class{track.table})
    \item slice (\linkS4class{track.table})
    \item trace (\linkS4class{envRefClass})
    \item types (\linkS4class{refTable})
    \item untrace (\linkS4class{envRefClass})
    \item usingMethods (\linkS4class{envRefClass})
  }
}
