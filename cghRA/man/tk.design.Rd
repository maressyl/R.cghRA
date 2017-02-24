\name{tk.design}
\alias{tk.design}
\title{
  Interactive cghRA design processing
}

\description{
  This function provides a Tcl-Tk interface to import a CGH array design file into a \code{\linkS4class{cghRA.design}} object and apply various cghRA tools on it.
}

\usage{
  tk.design(organism = "Human", assembly = "GRCh37",
    chromosomes = "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y",
    chromFiles = "", restrictionSites = "AluI=AG|CT, RsaI=GT|AC", globalTopLevel,
    localTopLevel)
}

\arguments{
  \item{organism}{
    Single character value, default value for the Organism field.
}
  \item{assembly}{
    Single character value, default value for the Assembly field.
}
  \item{chromosomes}{
    Single character value, default value for the Chromosomes field.
}
  \item{chromFiles}{
    Character vector, default chromosome files.
}
  \item{restrictionSites}{
    Single character value, default value for the Restriction sites field.
}
  \item{globalTopLevel}{
    This argument should be filled only when embedding this Tcl-Tk interface in an other. It is the top level of the embedding interface, generally a call to \code{\link{tktoplevel}}.
}
  \item{localTopLevel}{
    This argument should be filled only when embedding this Tcl-Tk interface in an other. It is the local top level to use to build this interface, generally a \code{\link{tkframe}} or \code{\link{ttkframe}}.
}
}

\author{
  Sylvain Mareschal
}

\seealso{
  \code{\link{tk.cghRA}}, \code{\linkS4class{cghRA.design}}, \code{\link{Agilent.design}}, \code{\link{custom.design}}
}
