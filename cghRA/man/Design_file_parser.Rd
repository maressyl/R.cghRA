\name{Design file parser}
\alias{Design file parser}
\alias{Agilent.design}
\alias{custom.design}
\title{
  Design file parser
}

\description{
  These functions are examples of design file parsers, as can be used directly or by \code{\link{tk.design}} to produce a \code{\linkS4class{cghRA.probes}} object from a CGH design file.
}

\usage{
  Agilent.design(file, name = NULL, organism = as.character(NA),
    assembly = as.character(NA), chromosomes = NULL, ...)
  custom.design(file, name = NULL, organism = as.character(NA),
    assembly = as.character(NA), chromosomes = NULL, ...)
}

\arguments{
  \item{file}{
    Single character value, path to the file to extract the design from (Agilent TDT design file for \code{Agilent.design}, CSV file as described below for \code{custom.design}).
}
  \item{name}{
    Single character value, the name of the design. \code{NULL} will generate an automatic design name with the array dimensions (e.g. "Agilent 125 x 50").
}
  \item{organism}{
    Single character value, the name of the organism studied by the current design.
}
  \item{assembly}{
    Single character value, the genome assembly version for probe coordinates.
}
  \item{chromosomes}{
    Character vector, the ordered list of the chromosome names for the design organism. If \code{NULL} the \code{factor} levels of the chrom column will be extracted, if not \code{chromosomes} will be used as levels to coerce the chrom column to \code{factor}.
}
  \item{\dots}{
    Further arguments are ignored by \code{Agilent.design} and \code{custom.design}, but can be used by other design file parsers.
  }
}

\details{
  As the package was developed with Agilent arrays, only the corresponding parser and a generic one are currently provided. Parsing design files from other brands can be achieved providing a custom design file parser suiting the manufacturer file format. Common brand file parsers may be added in the future, if you developed one (or need one to be developed) and wish it to be added to the package, please contact the package maintainer.
  
  "Custom" files must be CSV files, using tabulations as column separators, periods as decimal separators and a first row naming columns. No comment line is allowed, and cell content protection (quoting) can be performed using double-quotes. The mandatory columns are "chrom" (character), "start" (integer) and "end" (integer), describing the genomic location of each probe in the design. Additionally it is recommended to provide "strand" ("+", "-" or \code{NA}), "id" (an integer ID that will be used to match probes between design and data files), "name" (character), "row" and "col" (integers, the physical position of the probe on the slide). Further columns will be stored as provided.
}

\value{
  An object of class \code{\linkS4class{cghRA.design}}.
}

\author{
  Sylvain Mareschal
}

\seealso{
  \code{\link{cghRA.design-class}}, \code{\link{tk.design}}
}
