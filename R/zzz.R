.First.lib <- function(lib, pkg)
{
  library.dynam("randomForest", pkg, lib)
  Rvers <- sapply(R.Version()[c("major", "minor")], as.numeric)
  if (Rvers["major"] > 1 || Rvers["minor"] >= 9) {
    require(stats)
  } else {
    require(mva)
  }
}
