".First.lib" <-
function(lib, pkg) {
    library.dynam("mvpart", pkg, lib)
#    cat(" mvpart package loaded: extends rpart to include\n") 
#    cat(" multivariate and distance-based partitioning\n\n")
}

