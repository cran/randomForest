.onAttach <- function(libname, pkgname) {
    RFver <- if (as.numeric(R.version$major) < 2 && 
                 as.numeric(R.version$minor) < 9.0)
      package.description("randomForest")["Version"] else
    packageDescription("randomForest")$Version
    cat(paste("randomForest", RFver), "\n")
    cat("Type rfNews() to see new features/changes/bug fixes.\n")
}
