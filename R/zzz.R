.First.lib <- function(lib, pkg)
{
    library.dynam("randomForest", pkg, lib)
    require(mva)
}
