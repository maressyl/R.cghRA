\name{Probe file parser}
\alias{Probe file parser}
\alias{Agilent.probes}
\alias{custom.probes}
\title{
  Probe file parser
}

\description{
  These functions are examples of probe file parsers, as requested by \code{\link{process}} to produce a \code{\linkS4class{cghRA.probes}} object from a CGH array data file.
}

\usage{
  Agilent.probes(
    file,
    columns = c(
      rFin = "rProcessedSignal",
      gFin = "gProcessedSignal",
      flag_rIsSaturated = "rIsSaturated",
      flag_gIsSaturated = "gIsSaturated",
      flag_rIsFeatNonUnifOL = "rIsFeatNonUnifOL",
      flag_gIsFeatNonUnifOL = "gIsFeatNonUnifOL",
      flag_rIsBGNonUnifOL = "rIsBGNonUnifOL",
      flag_gIsBGNonUnifOL = "gIsBGNonUnifOL",
      flag_rIsFeatPopnOL = "rIsFeatPopnOL",
      flag_gIsFeatPopnOL = "gIsFeatPopnOL",
      flag_rIsBGPopnOL = "rIsBGPopnOL",
      flag_gIsBGPopnOL = "gIsBGPopnOL"
    ),
    ...
  )
  custom.probes(file, columns = NULL, ...)
}

\arguments{
  \item{file}{
    Single character value, path to the file to extract the design from (Agilent Feature Extraction file).
}
  \item{columns}{
    Character vector defining the columns to extract, the names are the names to use in the \code{\linkS4class{cghRA.probes}} object while the values are the names used in the Feature Extraction file.
}
  \item{\dots}{
    Further arguments are ignored by \code{Agilent.probes} and \code{custom.probes}, but can be used by other probe file parsers.
}
}

\details{
  As the package was developped with Agilent arrays, only the corresponding parser and a generic one are currently provided. Parsing arrays from other brands can be achieved providing a custom probe file parser suiting the manufacturer file format. Common brand file parsers may be added in the future, if you developped one (or need one to be developped) and wish it to be added to the package, please contact the package maintainer.
  
  As this function will be exported for parallel computing, dependencies need to be explicit : packages need \code{\link{library}} calls (even the core ones) or usage of \code{\link{::}} operators and sub-functions should be declared inside the parser body.
  
  "Custom" files must be CSV files, using tabulations as column separators, periods as decimal separators and a first row naming columns. No comment line is allowed, and cell content protection (quoting) can be performed using double-quotes. The mandatory columns are "id" (an integer ID that will be used to match probes between design and data files) and "logRatio" (numeric). Additionally one can provide boolean columns starting with "flag_", to be used as probe filters by \code{\link{process.mask}} during the array processing. Further columns will be stored as provided.
}

\value{
  An object of class \code{\linkS4class{cghRA.probes}}.
}

\author{
  Sylvain Mareschal
}

\seealso{
  \code{\link{cghRA.probes-class}}
}
