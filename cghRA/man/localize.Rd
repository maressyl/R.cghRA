\name{localize}
\alias{localize}
\alias{blatInstall}
\title{
  Localize CGH probes in a genome
}

\description{
  \code{localize} returns genomic coordinates (chromosome, strand, starting position, ending position) of a set of probes into a given genome. It relies on the external \code{Blast-Like Alignment Tool} to perform fuzzy both-strands matching, and provides various filters suitable to CGH probes.
  
  \code{blatInstall} needs to be executed once after the R package installation in order to use \code{localize}.
}

\usage{
  blatInstall(blat, cygwin)

  localize(probeFile, chromFiles, chromPattern = "^(.+)\\\\.[^\\\\.]+$",
    blatArgs = character(0), rawOutput = FALSE, noMulti = TRUE, noOverlap = TRUE,
    noPartial = TRUE, verbose = 2)
}

\arguments{
  \item{blat}{
    Single character value, path to the BLAT executable file to use for localization.
}
  \item{cygwin}{
    Single character value, path to the cygwin1.dll file that might be needed to run BLAT on Windows.
}
  \item{probeFile}{
    Single character value, path to a multi-fasta file describing the probes to compute the bias for. FASTA comments are used as probe names, and should be unique.
}
  \item{chromFiles}{
    Character vector, paths to chromosome sequences (a single fasta file for each chromosome).
}
  \item{chromPattern}{
    Single character value, a regular expression to be used for chromosome name extraction from \code{chromFiles}. It needs to capture a single value for replacement, default value will use the base names of the files without extension as chromosome names.
}
  \item{blatArgs}{
    Character vector, arguments to be passed to BLAT ("name=value" or "-flag"). See the BLAT documentation in 'References' for further details.
}
  \item{rawOutput}{
    Single logical value, whether to return the merged BLAT output or the processed one (see 'Value'). Notice raw output is not filtered.
}
  \item{noMulti}{
    Single logical value, whether to filter out probes located in multiple genomic positions or not. Ignored if \code{rawOutput}.
}
  \item{noOverlap}{
    Single logical value, whether to filter out overlapping probes or not (when two overlapping probes are detected, both are discarded). Ignored if \code{rawOutput}.
}
  \item{noPartial}{
    Single logical value, whether to filter out partial matches or not (they will still be used by other filters, to disable them completely consider using different BLAT arguments). Ignored if \code{rawOutput}.
}
  \item{verbose}{
    Single numeric value, the level of verbosity (0, 1 or 2).
}
}

\value{
  If \code{rawOutput}, \code{localize} returns the tabular section of merged \code{psLayout 3} file returned by BLAT (see the BLAT documentation in 'References' for further details).
  
  Else returns a \code{data.frame} with a row for each probe that was found and not filtered, ordered by \code{chrom}, \code{start} then \code{name} :
  \item{name}{Character, the probe names, as defined by comments in \code{probeFile}.}
  \item{chrom}{Character, the chromosomal location of the probe, as defined by the \code{chromNames} corresponding to the code{chromFiles} in which the probe matched.}
  \item{strand}{Character, "+" for a forward match, "-" for a reverse complement match.}
  \item{start}{Integer, the lower position of the probe in the chromosome. See 'Coordinate system'.}
  \item{end}{Integer, the upper position of the probe in the chromosome. See 'Coordinate system'.}
  \item{insertions}{Integer, amount of nucleotides inserted in the probe when refering to the chromosome sequence.}
  \item{deletions}{Integer, amount of nucleotides deleted in the probe when refering to the chromosome sequence.}
  \item{mismatches}{Integer, amount of mismatching nucleotides between probe and chromosome sequence.}
  \item{freeEnds}{Integer, amount of nucleotides at probe extremities ignored in the alignment.}
}

\section{Coordinate system}{
  When \code{rawOutput} is \code{FALSE}, coordinates begin at 1, both boundaries are comprised in the sequence and length can be computed as \code{end - start + 1} (\code{Biostrings} behavior).
  
  When \code{rawOutput}, refer to BLAT specifications (See 'References').
  
  In both cases, backward matches (strand = "-") are expressed in forward coordinates (start < end) (BLAT behavior).
}

\section{BLAT installation}{
  BLAT relies on a single executable file, so installation is straight-forward.
  
  Download the executable file or compile it for your computer architecture, then simply use the \code{blatInstall} function to copy it to the proper package folder for further uses. Precompiled executables for various systems can be found on the author website (see 'References'), as part of the BlatSuite (only 'blat.exe' or 'blat' is needed).
  
  \subsection{Windows specificities}{
    Running BLAT on Windows needs Cygwin. You can install Cygwin entirely on your system (see 'References'), or download the "cygwin1.dll" file and provide it to \code{blatInstall}, as it is the only Cygwin component needed. DLL is a common format for informatic viruses, so be sure of the website you download this file from. You can safely (no guarantee !) download it from the official website (see 'References') mirrors, they generally keep compressed archives in /release/cygwin in which you can find the DLL (in /usr/bin).
  }
}

\references{
  BLAT is an open-source software freely available for academic, nonprofit and personal use. See the FAQ for further details.
  \href{http://genome.ucsc.edu/FAQ/FAQblat.html}{FAQ}, \href{http://genome.ucsc.edu/goldenPath/help/blatSpec.html}{specifications}, \href{http://hgwdev.cse.ucsc.edu/~kent/src/}{source code} and \href{http://hgwdev.cse.ucsc.edu/~kent/exe/}{executables}
  
  Cygwin is a free and open-source software under GNU General Public Licencing.
  \href{http://www.cygwin.com/}{Official website}
}

\author{
  Sylvain Mareschal
}

\seealso{
  \code{\link{bias}}
}
