.onAttach <- function(libname, pkgname) {
    RFver <- read.dcf(file=system.file("DESCRIPTION", package=pkgname), 
                      fields="Version")
    cat(paste(pkgname, RFver, "\n"))
    cat("Type rfNews() to see new features/changes/bug fixes.\n")
}
