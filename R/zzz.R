.onAttach <- function(libname, pkgname) {
    desc <- packageDescription(pkgname)
    cat(paste(desc$Package, desc$Version), "\n")
    cat("Type rfNews() to see new features/changes/bug fixes.\n")
}
