".First.lib" <-
function(lib, pkg)
{
  library.dynam("evd", package = pkg, lib.loc = lib)
  return(invisible(0))
}

