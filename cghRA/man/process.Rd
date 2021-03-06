\name{process}
\alias{tk.process}
\alias{process}
\alias{process.applyModel}
\alias{process.cnvScore}
\alias{process.core}
\alias{process.default}
\alias{process.export}
\alias{process.fill}
\alias{process.filter}
\alias{process.fittest}
\alias{process.fixLast}
\alias{process.log}
\alias{process.mask}
\alias{process.modelize}
\alias{process.parse}
\alias{process.probes}
\alias{process.regions}
\alias{process.replicates}
\alias{process.segment}
\alias{process.spatial}
\alias{process.waca}
\title{
  cghRA array processing
}

\description{
  These functions implement the cghRA workflow, as a sequence of \code{process} subfunction calls. Each of them rely on \code{\linkS4class{cghRA.array}} and \code{\linkS4class{cghRA.regions}} methods, so custom processing can be easily achieved using them directly if the \code{steps} argument is not flexible enough to your purpose.
  
  Custom steps can be added as well on the model of existing ones, defining a function called \code{process.NAME} and adding "NAME" to the \code{steps} vector during the call to \code{process}. Step functions need to handle at least an \code{input} parameter which will be returned directly by the previous step, thus forming a pipeline.
  
  The \code{tk.process} function is a wrapper for \code{process}, built around a Tcl-Tk interface for more user-friendliness.
  
  The \code{process} function is a multi-core command line interface that will dispatch its arguments to individual \code{process.core} calls, and should be the prefered entry point even on single core computers. \code{process.log} is a wrapper to \code{process.core} which captures warnings and errors into a log file.
  
  The \code{process.default} function is a common way for \code{process} and \code{tk.process} to obtain default values for complex arguments like 'segmentArgs' and 'modelizeArgs'. It can be used to obtain the profiles proposed by \code{tk.process} in \code{process}.
}

\usage{
  process(inputs, logFile = "process.log", cluster = NA, output = FALSE, ...)
  process.log(..., logFile)
  process.core(input, inputName, output, steps = c("parse", "mask", "replicates", "waca",
    "export", "spatial", "segment", "fill", "fixLast", "modelize", "export", "fittest",
    "export", "applyModel", "export"), ...)
  
  process.applyModel(input, ...)
  process.cnvScore(input, design, dgv.map, cnvScoreCol = "cnvScore", ...)
  process.export(input, outDirectory, ...)
  process.fill(input, ...)
  process.filter(input, filter = NULL, ...)
  process.fittest(input, ...)
  process.fixLast(input, design, ...)
  process.mask(input, ...)
  process.modelize(input, modelizeArgs = process.default("modelizeArgs"), ...)
  process.parse(input, design, probeParser = Agilent.probes, probeArgs = list(), ...)
  process.probes(input, design, ...)
  process.regions(input, ...)
  process.replicates(input, replicateFun = stats::median, ...)
  process.segment(input, segmentArgs = process.default("segmentArgs"), ...)
  process.spatial(input, outDirectory, ...)
  process.waca(input, ...)
  
  tk.process(globalTopLevel, localTopLevel)
  process.default(argName, profileName)
}

\arguments{
  \item{inputs}{
    List of \code{input} to dispatch to each node (preferably named). The default workflow expects it to be a character vector naming raw data files to be parsed.
}
  \item{logFile}{
    Single character value, the path to the log file to produce with messages, warnings and errors. If the file already exists, it will be emptied first. The behavior when \code{logFile} is set to \code{NA} or "" depends on \code{cluster}: if \code{cluster} is \code{FALSE} (unparallelized mode), messages and errors will be passed to the R console rather than logged in a file, if \code{cluster} is anything else they will be silently ignored.
}
  \item{cluster}{
    Arguments to be passed to \code{\link{makeCluster}} as a list, for parallel processing (requires the optionnal \code{parallel} package). Remote machines are not handled properly in the current version of \code{process}, you should limit to "spec" defining how many processors can be used on the local machine as an integer value. The \code{FALSE} value requires an unparallelized mode, slower but more suitable for error tracking. The \code{NA} default value tries to detect the CPU count on the local machine if \code{parallel} is installed, else switches to unparallelized mode.
}
  \item{output}{
    Single logical value, whether to return the output of the last step or not (see the \code{Value} section)).
}
  \item{\dots}{
    Further arguments to be passed to \code{process} sub-functions, depending on the \code{steps} choosen (see below). The default workflow expects at least \code{design} and \code{outDirectory} to be provided.
}
  \item{input}{
    A single input to process on one node. The default workflow expects it to be a single character value naming a raw data file to be parsed.
}
  \item{inputName}{
    Single character value, the name of the input currently processed (for logging only).
}
  \item{steps}{
    Ordered character vector, naming the processing steps to apply. Custom steps can be named as well, as long as a function named "process.[step]" exists in the global environment. Each step will take as input the output of the previous step, the first step taking the value of the \code{input} argument as input.
}
  \item{probeParser}{
    The function to parse \code{probeFiles} into \code{\linkS4class{cghRA.probes}} objects, such as \code{\link{Agilent.probes}} for Agilent FeatureExtraction arrays.
}
  \item{probeArgs}{
    A list of arguments to pass to \code{probeParser} (apart from 'file' which is always provided).
}
  \item{design}{
    Single character vector, the path and name of the RDT design file, as produced by \code{\link{tk.design}}.
}
  \item{replicateFun}{
    The function to apply to replicate groups, if the "replicate" step is to be applied. This function must use a vector of numeric values (logRatios) as input, and return a single representative value (typically \code{median} or \code{mean}).
}
  \item{outDirectory}{
    Single character value, the directory in which produce the output files.
}
  \item{segmentArgs}{
    Character vector, the arguments to be passed to the \code{DNAcopy} method of the \code{\linkS4class{cghRA.array}} class. Arguments are defined as a character string that will be parsed, multiple values define multiple segmentation profiles to apply sequentially.
}
  \item{modelizeArgs}{
    Single character value, the arguments to be passed to the \code{model.auto} method of the \code{\linkS4class{cghRA.array}} class. Arguments are defined as a character string that will be parsed.
}
  \item{dgv.map}{
    Single character value, the path and name of a RDS file containing a \code{\linkS4class{segmentMap}} object to be passed to \code{\link{cnvScore}}.
}
  \item{cnvScoreCol}{
    Single character value, the name of the column that should be produced (or replaced) with the newly computed CNV score.
}
  \item{filter}{
    An R expression referring to columns as variables and returning a subsetting vector (either logical or numeric).
}
  \item{argName}{
    Single character value, 'segmentArgs' or 'modelizeArgs', the argument to get the default value for. If missing, the list of profiles and arguments handled is returned.
}
  \item{profileName}{
    Single character value, altering the default values returned. If missing, the default profile is returned.
}
  \item{globalTopLevel}{
    This argument should be filled only when embedding this Tcl-Tk interface in an other. It is the top level of the embedding interface, generally a call to \code{\link{tktoplevel}}.
}
  \item{localTopLevel}{
    This argument should be filled only when embedding this Tcl-Tk interface in an other. It is the local top level to use to build this interface, generally a \code{\link{tkframe}} or \code{\link{ttkframe}}.
}
}

\section{Processing steps}{
  The complete workflow involves the following steps :
  \describe{
    \item{parse}{Read a raw data file and return a \code{cghRA.array} object.}
    \item{probes}{Read a \code{cghRA.probes} object stored in a RDT file and return a \code{cghRA.array} object.}
    \item{regions}{Reads one or many \code{cghRA.regions} file(s) stored in RDT file(s).}
    \item{mask}{Discard flagged probes (saturated, high background ...) in a \code{cghRA.array} object. Any \code{TRUE} value in a column whose name begins with "flag_" is enough to discard a probe (turn its logRatio into \code{NA}. See the cghRA.array$maskByFlag() method for further details.}
    \item{replicates}{Replace replicated probe groups (same "name") by a single representative value (all logRatios are turned to \code{NA} except from the first one which will hold the representative value). See the cghRA.array$replicates() method for further details.}
    \item{waca}{Apply the WACA algorithm to the logRatios. See the cghRA.array$WACA() method for further details.}
    \item{spatial}{Produce a PNG file to visually check spatial biases. See the cghRA.array$spatial() method for further details.}
    \item{segment}{Compute regions with similar logRatios along the genome, using the CBS algorithm. See the cghRA.array$DNAcopy() method for further details.}
    \item{fill}{Extend segments to the right to join consecutive segments. See the cghRA.regions$fillGaps() method for further details.}
    \item{modelize}{Fit a copy number model to segments, in order to convert logRatios to true copy numbers. If \code{segmentArgs} contains multiple values, each segmentation profile will lead to distinct "copies" and "regions" files numbered according to its position in \code{segmentArgs}. See the cghRA.regions$model.auto() method for further details.}
    \item{applyModel}{Convert a modelized \code{cghRA.regions} objects into \code{cghRA.copies}.}
    \item{fittest}{If multiple segmentation profiles have been used, select the fittest model ("copies" and "regions" files duplicated without number). For further details on the STM score used for fittest model selection, see the \code{model.auto} function of the cghRA.copies package.}
    \item{clean}{Erase "copies" and "regions" files of the different segmentation profiles tested, as "fittest" should have saved the best.}
  }
  
  The following steps are currently optional :
  \describe{
    \item{cnvScore}{Compute for each segment of a \code{track.table} object its likelihood to be a polymorphism (see \code{\link{cnvScore}}).}
    \item{filter}{Remove rows from a \code{refTable} object using an expression based on column names.}
  }
}

\value{
  \code{process.default} : if \code{argName} is provided it returns the default value for the queried argument, else a list of profiles available for each handled argument. When many profiles are handled, the first value in the list is the default one (returned when \code{profileName} is missing).
  
  \code{process} : if \code{output} is \code{TRUE}, the function invisibly returns the output of the last step of the pipeline (in most cases a list with one element or a sub-list for each sample). If \code{output} is \code{FALSE}, the function invisibly returns \code{inputs} untouched.
}

\author{
  Sylvain Mareschal
}

\seealso{
  \code{\link{tk.design}}, \code{\linkS4class{cghRA.array}}
}

