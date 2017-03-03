\name{cghRA-news}
\title{cghRA package NEWS}
\section{Changes in version 1.6.0 (2017-03-03)}{
  \itemize{
    \item{Manual page cleanup.}
    \item{Vignette addition.}
    \item{Non-base package functions are now properly identified, to comply with R CMD check.}
    \item{C entry points are now registered, to comply with R CMD check.}
    \item{'states' argument default values harmonized across the package.}
    \item{cghRA.copies$newSTEPS() replaces cghRA.copies$STEPS().}
    \item{Better rendering for plots in tk.modelize() for Linux and/or R >= 3.4.0.}
    \item{Better automatic sizing of the plot in tk.browse().}
    \item{cghRA.copies$STEPS() no longer crashes if 'nested' is not 'none' with an empty region table.}
    \item{tk.series() rebuilt, moving some elements to the new tk.annotate().}
    \item{White image at tk.modelize() start.}
    \item{segmentMap class added for cnvScore maps.}
    \item{cghMCR algorithm no longer available in cghRA.series.}
    \item{GADA segmentation no longer available in cghRA.array.}
  }
}
\section{Changes in version 1.5.1 (2016-09-14)}{
  \itemize{
    \item{cghRA.copies$segMerge() removed as conflicting with Rgb::track.table.}
    \item{process.replicates() no longer crashes when no duplicated probe is present.}
  }
}
\section{Changes in version 1.5.0 (2016-01-19)}{
  \itemize{
    \item{cnvScore(), map2design(), applyMap() and trace2track() added.}
    \item{"dgv.track" dataset and track.CNV.DGVsupp() added.}
    \item{chromEnd argument added to STEPS(), default behavior changed from TRUE to FALSE.}
    \item{newSTEPS() added, as a function and cghRA.series method.}
    \item{cghRA.array class code reformated.}
    \item{cghRA.array$chromosomes() added.}
  }
}
\section{Changes in version 1.4.1 (2015-09-11)}{
  \itemize{
    \item{NAMESPACE and dependencies updated to comply with new R policies.}
    \item{cghRA.series$show() fixed for empty series.}
    \item{cghRA.series$pool() now accepts a 'group' argument.}
    \item{cghRA.series$pool() now complies with Rgb::track.CNV types.}
    \item{'cluster' argument from process() modified, see the manual page.}
    \item{cghRA.series$SRA(), $LRA() and $STEPS() now produce tracks based on Rgb::draw.steps().}
  }
}
\section{Changes in version 1.4.0 (2015-08-05)}{
  \itemize{
    \item{process() rewritten for more flexibility.}
    \item{cghRA.copies, cghRA.regions and cghRA.probes packages merged into cghRA.}
    \item{segMerge() and UTILS.tk removed, as already defined in package Rgb.}
    \item{IRanges package dependency in localize() replaced by Rgb.}
  }
}
\section{Changes in version 1.3.2 (2015-07-22)}{
  \itemize{
    \item{process() can now be launched without log file.}
    \item{DNAcopy bug temporarly fixed in cghRA.array$DNAcopy().}
    \item{"CBS default" added to process.default() for 'segmentArgs'.}
  }
}
\section{Changes in version 1.3.1 (2014-06-05)}{
  \itemize{
    \item{'quiet' argument added to cghRA.series$pool().}
    \item{Agilent.design() now relies on cghRA.design() constructor.}
    \item{'colorFun' defaults for cghRA.regions and cghRA.copies redefined to comply with Rgb 1.0.0.}
    \item{cghRA.series$applyMethod() no longer crashes with methods returning NULL.}
    \item{cghRA.series$pool() now relies on Rgb::track.CNV class.}
  }
}
\section{Changes in version 1.3.0 (2014-01-28)}{
  \itemize{
    \item{track.lone occurences renammed in track.table, for Rgb compatibility.}
    \item{cghRA.copies now relies on its model's ploidy for default drawing origin.}
    \item{cghRA.copies no longer enforces a new default parameter for colorFun.}
    \item{Error in cghRA.series$pool() fixed.}
  }
}
\section{Changes in version 1.2.5 (2013-09-02)}{
  \itemize{
    \item{cghRA.probes now inherits from refTable, no more from track.part.}
    \item{cghRA.probes now have its own 'name' field, no more inherited.}
    \item{cghRA.array now inherits an independant 'name' field from the drawable class, no more linked to its 'probes' 'name' field.}
  }
}
\section{Changes in version 1.2.4 (2013-08-30)}{
  \itemize{
    \item{cghRA.array$initialize() and cghRA.array$defaultParams() updated to current Rgb standards (0.5.0).}
    \item{cghRA.array$maskByFlag() no longer raises error about data class.}
    \item{cghRA.array$replicates() no longer raises error about data class.}
    \item{parseKaryo() argument 'alteratedOnly' added, with a new default (previous behavior was like 'FALSE').}
    \item{parseKaryo() now returns a list with an additional clone summary.}
  }
}
\section{Changes in version 1.2.3 (2013-07-12)}{
  \itemize{
    \item{parseKaryo() added.}
    \item{cghRA.series$MCR() renammed into cghRA.series$STEPS().}
    \item{'mode' argument added to cghRA.regions$proportions().}
    \item{defParam() overloads replaced by defaultParams() overloads.}
    \item{Dependency 'cghRA.track' becomes 'Rgb'.}
  }
}
\section{Changes in version 1.2.2 (2013-05-13)}{
  \itemize{
    \item{Tk windows names and titles fixed on Gnome 3.}
    \item{Icons replaced by 16x16 and 32x32 GIF images.}
    \item{Tk style 'clam' is the new default on unix platforms.}
    \item{Tk file selection windows become transients.}
    \item{process.one() timers in logs fixed on unix platforms.}
    \item{X11 warnings silenced in tk.file() to fix 'X11 BadDrawable'.}
    \item{'karyotype' method moved from cghRA.array class to cghRA.regions.}
    \item{'mergeArms' and 'splitArms' renamed into 'eraseArms' and 'addArms' to avoid confusion.}
    \item{cghRA.series constructor now accepts RDT file path vectors.}
  }
}
\section{Changes in version 1.2.1 (2013-02-25)}{
  \itemize{
    \item{Agilent.design() now uses R NA values for 'NA.[number]' probe names.}
    \item{'chrom' argument related systematic error corrected in cghRA.array$GLAD().}
  }
}
\section{Changes in version 1.2.0 (2013-01-04)}{
  \itemize{
    \item{process() 'do...' arguments merged in a single 'steps' argument.}
    \item{process() 'output' argument replaced by 'outputDir'.}
    \item{process() 'fun' argument replaced by 'replicateFun'.}
    \item{process() 'replicateFun' type changed from function name to function.}
    \item{'Spatial' checkbox added to tk.process().}
    \item{process.one() no longer takes its arguments in the global environnement.}
  }
}
\section{Changes in version 1.1.1 (2012-12-19)}{
  \itemize{
    \item{Profiles added in tk.process().}
    \item{process.default() added.}
  }
}
\section{Changes in version 1.1.0 (2012-12-12)}{
  \itemize{
    \item{cghRA.copies class added.}
    \item{cghRA.regions$model.apply() no longer alter itself, but return a cghRA.copies object.}
    \item{getParameters() overloads replaced by defParam() overloads for the package classes.}
    \item{'maxElements' and 'label' default parameters added to cghRA.design class.}
    \item{Drawing parameters pre-set in cghRA.array$pool() track output.}
    \item{cghRA.array$pool() unified to use 'states' and 'values' rather than 'keep'.}
    \item{'Pool' track production added to tk.series().}
    \item{cghRA.array also implements getChromEnd() for optimisation.}
    \item{cghRA.regions$proportions() added.}
  }
}
\section{Changes in version 1.0.11 (2012-12-12)}{
  \itemize{
    \item{'segMerge' method added to cghRA.regions class.}
    \item{TK file choosing windows now remind last directory browsed.}
    \item{require() calls turned into library() ones to raise errors.}
    \item{GLAD segmentation and calling added to cghRA.array class.}
  }
}
\section{Changes in version 1.0.10 (2012-11-23)}{
  \itemize{
    \item{draw() overloads replaced by getParameters() overloads for the package classes.}
  }
}
\section{Changes in version 1.0.9 (2012-11-09)}{
  \itemize{
    \item{cghRA.regions$status() rebuilt from ground to mimic cghRA.regions::penetrance behavior.}
    \item{cghRA.series$pool() added.}
  }
}
\section{Changes in version 1.0.8 (2012-10-31)}{
  \itemize{
    \item{process() now accepts other Probe file parsers than Agilent.probes.}
    \item{process() 'column' argument removed, 'arrayFiles' replaced by 'probeFiles', 'probeParser' and 'probeArgs' added.}
    \item{process() and tk.process() arguments 'filePattern' and 'fileReplace' removed.}
    \item{process() now ignore 'spatial' step if 'row' and 'col' are not available.}
    \item{tk.process() 'probeParser' argument added.}
    \item{tk.process() 'fillGaps' argument is no longer checked by default, as the MCR algorithm can handle it.}
    \item{Agilent.probes() 'name' argument removed.}
  }
}
\section{Changes in version 1.0.7 (2012-10-26)}{
  \itemize{
    \item{Segmentation 'min.width' default value switched back to DNAcopy default (2 presently) in tk.process().}
    \item{'Keep all segmentations' behavior added to tk.process(), process() and process.one().}
    \item{'exclude' argument definition added to tk.process().}
    \item{tk.design() know produces '*.design.rdt' files, for smarter file selection.}
  }
}
\section{Changes in version 1.0.6 (2012-10-23)}{
  \itemize{
    \item{'index' field and buildIndex methods added to track.lone class, to fulfill cghRA.track::subtrack 1.0.5 updates.}
    \item{tk.series() MCR now produces a 'boundariesProbeCount' column.}
    \item{'tracks' argument added to parallelize, penetrance, MCR, SRA and LRA methods for cghRA.series.}
    \item{'Penalty' replaced by 'Allow nested' in tk.series().}
    \item{cghRA.series$hasState method replaced by cghRA.regions$status.}
    \item{cghRA.series$applyMethod .mode argument added, to control output type.}
    \item{cghRA.series$applyMethod nrow relies on 'arrays' names rather than 'name' for faster computation.}
  }
}
\section{Changes in version 1.0.5 (2012-10-12)}{
  \itemize{
    \item{Cytoband file selection added in tk.design().}
    \item{'...' was not passed through by cghRA.regions$fillGaps.}
  }
}
\section{Changes in version 1.0.4 (2012-10-05)}{
  \itemize{
    \item{'utils', 'stats', 'graphics' and 'grDevices' were not loaded on cluster members in process.one().}
    \item{Extra vertical space deplaced from top to bottom in tk.series().}
    \item{tk.cghRA() window name changed to 'cghRA - Data processor'.}
  }
}
\section{Changes in version 1.0.3 (2012-09-26)}{
  \itemize{
    \item{'blocking' argument added to tk.cghRA(), for stand alone version.}
  }
}
\section{Changes in version 1.0.2 (2012-09-21)}{
  \itemize{
    \item{'strand' column added to pass cghRA.regions$check() in tk.series().}
  }
}
\section{Changes in version 1.0.1 (2012-09-21)}{
  \itemize{
    \item{Error message reshaping in tk.design().}
    \item{Error message reshaping in process.one().}
    \item{'doFittest' was not exported from process() to process.one().}
    \item{'No fittest segmentation found' warning added to process.one().}
    \item{'fittest' section added to process.one() log file.}
    \item{Default 'nperm' in tk.process() changed to 25/alpha.}
  }
}
