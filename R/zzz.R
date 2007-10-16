.onAttach <- function(libname, pkgname) {
    RFver <- read.dcf(file=system.file("DESCRIPTION", package=pkgname), 
                      fields="Version")
    message(paste(pkgname, RFver))
    message("Type rfNews() to see new features/changes/bug fixes.")
}
