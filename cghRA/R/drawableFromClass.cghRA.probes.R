# Add cghRA.probes compatibility to Rgb::drawable.list$add()
# Author : Sylvain Mareschal <mareschal@ovsa.fr>

drawableFromClass.cghRA.probes <- function(track, design, ...) {
	# Corresponding design
	if(missing(design)) {
		design <- tk.file(
			title = "Choose the corresponding design file",
			typeNames = c("Standard design file", "RDT file", "Custom R drawable object"),
			typeExt = c("*.design.rdt", "*.rdt", "*.rds"),
			multiple = FALSE,
			mandatory = TRUE,
			type = "open"
		)	
	}
	
	# Parse design
	if(is.character(design) && length(design) == 1L) design <- readRDT(design)
	
	# Check elements
	if(!is(design, "cghRA.design")) stop("'design' must be a 'cghRA.design' object or the path to a RDT file")
	if(!is(track, "cghRA.probes"))  stop("'track' must be a 'cghRA.probes' object")
	
	# Merge as a cghRA.array
	arr <- cghRA.array(.design=design, .probes=track)
	
	return(arr)
}

