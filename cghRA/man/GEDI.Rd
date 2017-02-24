\name{GEDI}
\alias{GEDI}
\title{
  Gene Expression and Dosage Integrator
}

\description{
  This function implements the "Gene Expression and Dosage Integrator" CGH / transcriptome correlation, as described by Lenz et al.
}

\usage{
  GEDI(cgh, cgh.chrom, cgh.start, cgh.end, cgh.genes, expr, expr.genes,
    permutations = 1000, type = c("amplifications", "deletions"), quiet = FALSE)
}

\arguments{
  \item{cgh}{
    Logical matrix, with regions in rows and samples in columns. Alterated samples for a given region are to be \code{TRUE}, germline \code{FALSE} and other \code{NA}.
}
  \item{cgh.chrom}{
    Character vector, the chromosome location of the regions described in \code{cgh}.
}
  \item{cgh.start}{
    Integer vector, the starting position on the chromosome for the regions described in \code{cgh}.
}
  \item{cgh.end}{
    Integer vector, the ending position on the chromosome for the regions described in \code{cgh}.
}
  \item{cgh.genes}{
    Character vector, the names of the genes in each region described in \code{cgh}, separated by ", ". See the \code{cross} method of the \code{sliceable} class (in \code{Rgb} package) for an easy way to produce this, in combination with \code{track.NCBI_genes}.
}
  \item{expr}{
    Numeric matrix of gene expressions, with probesets in rows and samples in columns.
}
  \item{expr.genes}{
    Character vector, the names of the genes associated with each probeset described in \code{expr}, separated by ", ". Notice probesets associated with multiple genes will not be used, as they are not specific.
}
  \item{permutations}{
    Single integer value, the amount of permutations to use for score computation. Time consumption and score accuracy increases with this value.
}
  \item{type}{
    Single character value, describing the type of alterations studied (as the alternative hypothesis for the t-test depends on it).
}
  \item{quiet}{
    Single logical value, when \code{FALSE} a message will be sent for each region processing, in order to evaluate the processing time.
}
}

\value{
  Returns a \code{list} with the following elements :
  \item{gediScore}{Numeric vector with for each \code{cgh} row the proportion of permutated scores lesser than the observed one. The algorithm authors consider an association to be present if this score is greater than 0.9.}
  \item{gediGenes}{Character vector with for each \code{cgh} row the list of the genes used for the score computation (intersection of \code{cgh.genes} and \code{expr.genes} for the considered region).}
}

\references{
  Lenz G et al. "Molecular subtypes of diffuse large B-cell lymphoma arise by distinct genetic pathways".
  Proc Natl Acad Sci U S A. 2008 Sep 9;105(36):13520-5 (Supporting Information)
}

\author{
  Sylvain Mareschal
}
