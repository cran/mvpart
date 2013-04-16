"gdist" <-
function(x, method = "bray", keepdiag = FALSE , full = FALSE, sq = FALSE)
{
    METHODS <- c("manhattan", "euclidean", "canberra", "bray", "kulczynski", "gower", 
    "maximum", "binary", "chisq", "chord", "beta0", "beta1", "beta2")
    method <- pmatch(method, METHODS)
    if(is.na(method))
        stop("invalid distance method")
    N <- nrow(x <- as.matrix(x))
    if (method == 6) x <- scaler(x,col=c("min0","max1"))     
    if(method == 9) {
        rr <- apply(x, 1, sum)
        cc <- apply(x, 2, sum)
        x <- diag(1/sqrt(rr)) %*% x %*% diag(1/sqrt(cc))
        method <- 2
    }
    else if(method == 10) {
        mns <- sqrt(apply(x^2, 1, sum))
        x <- x/(mns * sqrt(2))
        method <- 2
    }
	else if(method > 10) method <- method - 2
    d <- .C("gdistance",
        x = as.double(x),
        nr = N,
        nc = ncol(x),
        d = double((N * (N - 1))/2),
        keepdiag = as.integer(FALSE),
        method = as.integer(method),
        PACKAGE="mvpart")$d
    attr(d, "Size") <- N
    class(d) <- "dist"
    if (full) d <- distfull(d)
    if (sq) d <- d^2
    d
}

