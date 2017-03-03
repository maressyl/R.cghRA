# CNV track from the new version of DGV (supporting variants)
# Author : Sylvain Mareschal <mareschal@ovsa.fr>
# Example file : http://dgvbeta.tcag.ca/dgv/app/downloads (Release Version)

track.CNV.DGVsupp = function(
		file,
		name = "DGV CNV (supporting variants)",
		quiet = FALSE,
		...
		)
	{
	# Parse DGV supporting variants
	if(!quiet) message("Parsing data file ...")
	col.names <-  c("variantaccession", "chr",       "start",   "end",     "varianttype", "variantsubtype", "reference", "pubmedid", "method",    "platform",  "mergedvariants", "supportingvariants", "mergedorsample" ,"frequency", "samplesize", "observedgains", "observedlosses", "cohortdescription", "genes",     "samples")
	colClasses <- c("character",        "character", "integer", "integer", "character",   "character",      "character", "integer",  "character", "character", "character",      "character",          "character",      "double",    "integer",    "integer",       "integer",        "character",         "character", "character")
	tab <- utils::read.table(file, header=TRUE, sep="\t", quote="\"", dec=".", col.names=col.names, colClasses=colClasses, comment.char="")
	if(!quiet) message(nrow(tab), " variants read")

	# Filter on type
	if(!quiet) message("Filtering on type ...")
	tab <- tab[ tab$varianttype == "CNV" ,]
	tab <- tab[ tab$observedgains > 0 | tab$observedlosses > 0 ,]

	# Filter out merged variants (cf PMID:16327809, technique concluding at series level)
	if(!quiet) message("Filtering on samples ...")
	tab <- tab[ tab$mergedorsample == "S" ,]

	# Filter on cohort size
	if(!quiet) message("Filtering on cohort size ...")
	tab <- tab[ tab$samplesize > 20 ,]

	# Filter on chromosome
	if(!quiet) message("Filtering on chromosome ...")
	tab <- tab[ tab$chr %in% c(1:22, "X", "Y") ,]

	### # Methods
	### methods <- unique(unlist(strsplit(unique(tab$method), split=",")))
	### retained <- c(
	### 	"Oligo aCGH", "SNP array", "ROMA",                          # Micro-arrays
	### 	"qPCR", "MLPA", "Read-depth_analysis",                      # Quantitative PCR (confirmation usually)
	### 	"Sequence_alignment", "Sequencing", "Paired-end_mapping",   # Sequencing ("Sequencing" seems to refer to NGS)
	### 	"PCR", "MassSpec",                                          # PCR product size (mass spectrometry PMID:17122850) (confirmation usually)
	### 	"Composite_approach"                                        # 42M SNP array anyway PMID:19812545
	### )
	### excluded <- c(
	### 	"BAC aCGH", "FISH", "Karyotyping",   # Low resolution
	### 	"OEA_assembly",                      # One-End Anchored fosmid clones (PMID:20440878), low resolution
	### 	"Optical_mapping"                    # Unusual technique (cell culture, restriction site based alignment), PMID:20534489
	### )
	### # "SNP_genotyping_analysis" ???
	### # "MCD_analysis" ???
	### # "Digital array" ??? -> co-occurs with "Oligo aCGH"
	### # "Merging" ??? -> co-occurs with "Oligo aCGH"

	# Filter on methods
	if(!quiet) message("Filtering on method ...")
	required <- c("Oligo aCGH", "SNP array", "ROMA", "Sequencing")
	exclude <- c("BAC aCGH", "FISH", "Karyotyping")
	splitMethod <- strsplit(tab$method, split=",", fixed=TRUE)
	selected <- sapply(X=splitMethod, FUN=function(x, req, exc){ any(!is.na(match(x, req))) && all(is.na(match(x, exc))) }, req=required, exc=exclude)
	tab <- tab[ selected ,]

	if(!quiet) message(nrow(tab), " variants kept")

	# Reformat
	if(!quiet) message("Reformatting ...")
	dgv <- track.CNV(
		name = tab$variantaccession,
		chrom = factor(tab$chr, levels=c(1:22, "X", "Y")),
		start = tab$start,
		end = tab$end,
		strand = factor(rep("+", nrow(tab)), levels=c("-","+")),
		type = factor(tab$variantsubtype),
		PMID = tab$pubmedid,
		parent = tab$mergedvariants,
		.name = name,
		...
	)

	return(dgv)
}

