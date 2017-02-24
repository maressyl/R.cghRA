\name{parseKaryo}
\alias{parseKaryo}
\title{
  Parses a karyotype-like formula
}

\description{
  This function produces a \code{\linkS4class{cghRA.regions}} object from a simplified karyotype formula, associating copy numbers to numeric coordinates.
}

\usage{
  parseKaryo(formula, bandTrack, name = as.character(NA), design = NULL,
    alteratedOnly = TRUE)
}

\arguments{
  \item{formula}{
    Single character value, the formula to be parsed. See 'Examples'.
}
  \item{bandTrack}{
    A \code{track.table} object with cytoband definition, as returned by the \code{track.UCSC_bands} function from the \code{Rgb} package.
}
  \item{name}{
    Single character value, to be used as name for the produced object.
}
  \item{design}{
    A \code{\linkS4class{cghRA.design}} object, or \code{NULL}. If provided, a \code{\linkS4class{cghRA.copies}} object will be produced, using \code{design} to compute probe content of each region. Else, a \code{track.table} object will be returned.
}
  \item{alteratedOnly}{
    Single logical value, if \code{TRUE} normal clones (2n without alteration) will not be averaged with alterated clones for the final copy amount computation. If all clones are normals, a normal genome will be returned anyway.
}
}

\value{
  Returns a \code{list} with two elements : "clones" and "copies".
  
  "clones" is a summary of the clones found in the formula as an integer value, with mitosis counts as values and ploidy as names.
  
  "copies" is a \code{track.table}-inheriting object with genomic regions of distinct copy numbers. If \code{design} is provided, the object is a \code{\linkS4class{cghRA.copies}} object, else a \code{track.table} object.
}

\author{
  Sylvain Mareschal
}

\seealso{
  \code{\linkS4class{cghRA.copies}}
}

\examples{
  \dontrun{
    karyo <- paste(
      "111<5n>,6(1qt-p11),4(1p11-pt),4(2),8(3),4(4),6(5),6(6pt-q22),6(6q26-qt),",
      "2(6q22-q26),6(7pt-q31),3(7q31-qt),6(9),4(10),4(11),4(12),6(13),4(14),",
      "4(15pt-q22),2(15q22-qt),2(16),4(17),6(18),4(19),4(21),4(22) [6] ; 46<2n> [7]",
      collapse = ""
    )
    parseKaryo(karyo, bandTrack)
  }
}

