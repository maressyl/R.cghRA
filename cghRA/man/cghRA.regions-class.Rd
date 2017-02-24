\name{cghRA.regions-class}
\Rdversion{1.1}
\docType{class}
\alias{cghRA.regions-class}

\title{Class \code{"cghRA.regions"}}
\description{
  This class is intended to store a list of genomic segments produced by a segmentation algorithm, with a mean log-ratio for each segment.
  
  Objects from this class are intended to be produced by the \code{DNAcopy} method of the \code{\linkS4class{cghRA.array}} class, or the \code{\link{cghRA.regions}} constructor. Producing such objects is part of the \code{\link{process}} function and its interfaced version \code{\link{tk.process}}, found in \code{\link{tk.cghRA}}.
}

\author{
  Sylvain Mareschal
}

\seealso{
  \code{\link{cghRA.array-class}}, \code{\link{process}}, \code{\link{tk.process}}, \code{\link{refTable-class}}
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
  \describe{
    \item{\code{model}:}{Numeric vector, storing the parameters and fitness of a copy-number model. See \code{\link{model.auto}} for further details on the components.}
    \item{\code{modelizeCall}:}{The R \code{call} which produced the stored copy-number model.}
    \item{\code{segmentCall}:}{The R \code{call} which produced the segments stored in the object.}
  }
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
    \item{\code{fillGaps(...)}:}{
      Apply the fillGaps() function to extend regions in order to fill inter-segment gaps.
    }
    \item{\code{karyotype(bandTrack, value = , thresholds = , precision = )}:}{
      Returns a karyotype formula of altered regions.\cr
      - \strong{bandTrack}    : a track.table object, as returned by track.UCSC_bands().\cr
      - \strong{value}        : column to use to select altered regions.\cr
      - \strong{thresholds}   : length 2 numeric vector defining altered values.\cr
      - \strong{precision}    : single integer value from 1 to 4, amount of digits to consider in banding.
    }
    \item{\code{model.apply(...)}:}{
      Call the model.apply() function to produce a cghRA.copies object with predicted copy number for each region.
    }
    \item{\code{model.auto(save = , ...)}:}{
      Call the model.auto() function to automatically fit a copy-number prediction model.\cr
      - \strong{save}   : single logical value, whether to save the model or only return it
    }
    \item{\code{modelized()}:}{
      Does the object embed a complete model or not
    }
    \item{\code{model.test(...)}:}{
      Call the model.test() function to plot the current copy-number model.
    }
    \item{\code{proportions(chrom = , value = , states = , mode = )}:}{
      Returns the proportion of the chromosomes in given states (in bp involved).\cr
      - \strong{chrom}    : character vector, chromosome location of the regions to query. Consider track.table$eraseArms() to focus on chromosome arms.\cr
      - \strong{value}    : single character value, name of the column to use for state assignation.\cr
      - \strong{states}   : list of states, see penetrance help page for details.
    }
    \item{\code{status(chrom, start, end, value = , na = , fuzzy = , states = )}:}{
      Returns the copy states in various windows, mimicing penetrance behavior.\cr
      - \strong{chrom}    : character vector, chromosome location of the regions to query.\cr
      - \strong{start}    : integer vector, starting position on the chromosome for the regions to query.\cr
      - \strong{end}      : integer vector, ending position on the chromosome for the regions to query.\cr
      - \strong{value}    : single character value, name of the column to use for state assignation.\cr
      - \strong{na}       : single character value, see penetrance() help page for details ('false' is not handled).\cr
      - \strong{fuzzy}    : single logical value, whether to assign the state when some sub-regions are out or not.\cr
      - \strong{states}   : list of states, see penetrance help page for details.
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
