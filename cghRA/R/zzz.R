.onLoad <- function(libname, pkgname) {
	library.dynam(pkgname, pkgname, libname)
}
.onUnload <- function(libpath) {
	library.dynam.unload("cghRA", libpath)
}
