.onLoad <- function(libname, pkgname) {
  library.dynam("dfertility", pkgname, libname)
}

.onUnload <- function(libpath) {
  library.dynam.unload("dfertility", libpath)
}
