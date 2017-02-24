# Moves BLAT files to the package location for straight BLAT usage
# http://genome.ucsc.edu/goldenPath/help/blatSpec.html
# Author : Sylvain Mareschal <mareschal@ovsa.fr>

blatInstall = function(
		blat,
		cygwin
		)
	{
	# Directory
	package = find.package("cghRA", quiet=TRUE, verbose=FALSE)
	dir.create(sprintf("%s/blat", package), showWarnings=FALSE)
	
	# Cygwin DLL
	message("Cygwin DLL : ", appendLF=FALSE)
	if(.Platform$OS.type == "windows") {
		# Windows
		if(missing(cygwin)) {
			message("WARNING   (missing and required if Cygwin is not installed)")
		} else {
			if(basename(cygwin) != "cygwin1.dll") {
				message("ERROR     (incorrect file name, should be 'cygwin1.dll')")
			} else if(!file.exists(cygwin)) {
				message("ERROR     (file not found)")
			} else if(!suppressWarnings(file.copy(cygwin, sprintf("%s/blat/cygwin1.dll", package), overwrite=TRUE))) {
				message("ERROR     (copy failure at file system level)")
			} else {
				message("OK        (valid and copied)")
			}
		}
	} else {
		# UNIX
		if(missing(cygwin)) message("OK        (missing and not required)")
		else                message("OK        (provided and ignored)")
	}
	
	# BLAT exe
	message("BLAT exe   : ", appendLF=FALSE)
	if(missing(blat)) {
		message("ERROR     (not provided)")
	} else if(!file.exists(blat)) {
		message("ERROR     (file not found)")
	} else if(!suppressWarnings(file.copy(blat, sprintf("%s/blat/blat.exe", package), overwrite=TRUE))) {
		message("ERROR     (copy failure at file system level)")
	} else {
		message("OK        (valid and copied)")
	}
}
