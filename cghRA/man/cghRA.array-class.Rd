\name{cghRA.array-class}
\Rdversion{1.1}
\docType{class}
\alias{cghRA.array-class}

\title{Class \code{"cghRA.array"}}
\description{
  This class is the main component of the cghRA object-oriented package. Each CGH array must be stored in a distinct cghRA.array object.
  
  Objects from this class should always be produced by the \code{\link{cghRA.array}} constructor.
  
  This class is a hub, it provides methods to apply various CGH analysis tools in a straight-forward way.
  
  The Reference Class system is used notably to share designs objects between arrays, as several arrays may have values for the same probes.
}

\author{
  Sylvain Mareschal
}

\seealso{
  \code{\link{cghRA.array}}
  
  \code{\link{cghRA.series-class}}, \code{\link{cghRA.design-class}}, \code{\link{cghRA.probes-class}}, \code{\link{cghRA.regions-class}}
}

\keyword{classes}
\section{Extends}{
  Class \code{\linkS4class{crossable}}, directly.\cr
  Class \code{\linkS4class{sliceable}}, by class \code{\linkS4class{crossable}}, distance 2.\cr
  Class \code{\linkS4class{drawable}}, by class \code{\linkS4class{crossable}}, distance 3.

  All reference classes extend and inherit methods from \code{\linkS4class{envRefClass}}.
}
\section{Fields}{
  \describe{
    \item{\code{assembly}:}{Single \code{character} value, the assembly version for the coordinates stored in the object. Must have length 1, should not be \code{NA}.}
    \item{\code{design}:}{Object of class \code{\linkS4class{cghRA.design}}}
    \item{\code{organism}:}{Single \code{character} value, the name of the organism whose data is stored in the object. Must have length 1, should not be \code{NA}.}
    \item{\code{probes}:}{Object of class \code{\linkS4class{cghRA.probes}}.}
  }
  The following fields are inherited (from the corresponding class):
  \itemize{
    \item name (\linkS4class{drawable})
    \item parameters (\linkS4class{drawable})
  }
}
\section{Methods}{
  \describe{
    \item{\code{as.CNA()}:}{
      Returns a CNA object (DNAcopy) with the object content.
    }
    \item{\code{as.profileCGH(chrom = , quiet = )}:}{
      Returns a profileCGH object (GLAD) with the object content.\cr
      - \strong{chrom}   : single character value defining how to deal with chromosome names :\cr
                'merged' forces chromosome arms to be merged (as chromosome arms are not handled)\cr
                'levels' converts chromosome to integers (can be deceiving for factors)\cr
      - \strong{quiet}   : single logical value, whether to warn for factor to integer conversion or not.
    }
    \item{\code{DLRS(method = , na.rm = )}:}{
      Computes the Derivative Log Ratio Spread from the probes.\cr
      - \strong{method}   : 'agilent' or 'original', implying distinct formulas.
    }
    \item{\code{DNAcopy(smooth = , ...)}:}{
      Apply the Circular Binary Segmentation, as implemented in DNAcopy, and return a cghRA.regions object.\cr
      - \strong{smooth}      : a list of arguments to be passed to smooth.CNA(), TRUE to use the default parameters or FALSE to skip smoothing.\cr
      - \strong{...}         : arguments to be passed to segment().
    }
    \item{\code{extract(i = , j = )}:}{
      Extracts values from 'probes' and 'design' into a data.frame.\cr
      - \strong{i}   : row selection, see the R5Table method for further details.\cr
      - \strong{j}   : column selection, see the R5Table method for further details.
    }
    \item{\code{GADA(...)}:}{
      Apply the Genome Alteration Detection Analysis, as implemented in GADA, and return a cghRA.regions object.\cr
      - \strong{smooth}   : a list of arguments to be passed to smooth.CNA(), TRUE to use the default parameters or FALSE to skip smoothing.\cr
      - \strong{...}      : arguments to be passed to segment().
    }
    \item{\code{GLAD(chrom = , quiet = , output = , ...)}:}{
      Apply the Gain and Loss Analysis of Dna, as implemented in GLAD, and return a cghRA.regions object.\cr
      - \strong{chrom, quiet}   : to be passed to the as.profileCGH method.\cr
      - \strong{output}         : single character value defining the returned value :\cr
                         'regions' returns a cghRA.regions object with the segmented genome\cr
                         'raw' returns the glad() output\cr
                         'both' adds a 'cghRA.regions' element to the glad() output list to return both\cr
      - \strong{...}            : arguments to be passed to glad().
    }
    \item{\code{MAplot(pch = , cex = , xlab = , ylab = , ...)}:}{
      MA plot of all the probes.\cr
      - \strong{...}   : arguments to be passed to plot().
    }
    \item{\code{maskByFlag(flags = , pattern = , multiple = , na = )}:}{
      Replaces logRatios of flagged probes by NA.\cr
      - \strong{flags}      : character vector, the columns to coerce as boolean and use as flags.\cr
      - \strong{pattern}    : single logical value, whether to consider 'flags' as regular expressions or fixed values.\cr
      - \strong{multiple}   : mask a probe when 'all' its flag columns are TRUE or when 'any' is.
    }
    \item{\code{replicates(fun = , na.rm = , ...)}:}{
      Apply 'fun' to replicated probes (same name), masking all members but one.\cr
      - \strong{fun}   : single character value, the function to apply.\cr
      - \strong{...}   : to be passed to 'fun'.
    }
    \item{\code{spatial(filename = , palSize = , palEnds = , ...)}:}{
      Produces a spatial representation of the logRatios, to identify spatial biases.\cr
      - \strong{filename}   : single character value, the path to the PNG output.\cr
      - \strong{palSize}    : single integer value, the amount of color levels for logRatios. Should be lesser or equal to 254 to produce small PNG files.\cr
      - \strong{palEnds}    : character vector to be passed to colorRampPalette() for palette generation.
    }
    \item{\code{WACA()}:}{
      Apply the Waves aCGH Correction Algorithm (Lepretre et al. 2009) to the array logRatios.
    }
  }
  The following methods are inherited (from the corresponding class):
  \itemize{
    \item callParams (\linkS4class{drawable})
    \item callSuper (\linkS4class{envRefClass})
    \item check (\linkS4class{drawable}, overloaded)
    \item chromosomes (\linkS4class{drawable}, overloaded)
    \item copy (\linkS4class{envRefClass})
    \item cross (\linkS4class{crossable})
    \item defaultParams (\linkS4class{sliceable}, overloaded)
    \item draw (\linkS4class{sliceable})
    \item export (\linkS4class{envRefClass})
    \item field (\linkS4class{envRefClass})
    \item fix.param (\linkS4class{drawable})
    \item getChromEnd (\linkS4class{sliceable}, overloaded)
    \item getClass (\linkS4class{envRefClass})
    \item getName (\linkS4class{drawable})
    \item getParam (\linkS4class{drawable})
    \item getRefClass (\linkS4class{envRefClass})
    \item import (\linkS4class{envRefClass})
    \item initFields (\linkS4class{envRefClass})
    \item initialize (\linkS4class{drawable}, overloaded)
    \item setName (\linkS4class{drawable})
    \item setParam (\linkS4class{drawable})
    \item show (\linkS4class{sliceable}, overloaded)
    \item slice (\linkS4class{sliceable}, overloaded)
    \item trace (\linkS4class{envRefClass})
    \item untrace (\linkS4class{envRefClass})
    \item usingMethods (\linkS4class{envRefClass})
  }
}
