"scaler" <-
function(x, col = c("mean1", "max1", "min0", "ssq1", "range01", "zsc", "pa", "rank")[2],
    row = c("mean1", "max1", "min0", "ssq1", "range01", "zsc", "pa", "rank")[1])
{
fun = function (x, method, MARGIN) 
    {
    x <- as.matrix(x)
    switch(method, mean1 = {
       tmp <- apply(x, MARGIN, mean, na.rm = TRUE)
        x <- sweep(x, MARGIN, tmp, "/")
    }, max1 = {
        tmp <- apply(x, MARGIN, max, na.rm = TRUE)
        x <- sweep(x, MARGIN, tmp, "/")
    }, min0 = {
        tmp <- apply(x, MARGIN, min, na.rm = TRUE)
        x <- sweep(x, MARGIN, tmp, "-")
    }, ssq1 = {
        tmp <- apply(x^2, MARGIN, sum, na.rm = TRUE)
        tmp <- sqrt(tmp)
        x <- sweep(x, MARGIN, tmp, "/")
    }, range01 = {
        tmp <- apply(x, MARGIN, min, na.rm = TRUE)
        ran <- apply(x, MARGIN, max, na.rm = TRUE)
        ran <- ran - tmp
        x <- sweep(x, MARGIN, tmp, "-")
        x <- sweep(x, MARGIN, ran, "/")
    }, zsc = {
        if (MARGIN == 1) 
            x <- t(scale(t(x)))
        else x <- scale(x)
    }, pa = {
        tmp <- dim(x)
        nam <- dimnames(x)
        x <- as.numeric(x > 0)
        dim(x) <- tmp
        dimnames(x) <- nam
    }, rank = {
       x <- apply(x, MARGIN, rank)
       if (MARGIN == 1) 
       x <- t(x)
    })
    x
}
    METHODS <- c("mean1", "max1", "min0", "ssq1", "range01", "zsc", "pa", "rank")
    if (is.null(col) & is.null(row)) 
        cat("Scalings are",METHODS,"\n")
    if (!is.null(col)) {
    for (i in 1:length(col)){
        method <- match.arg(col[i], METHODS)
        x <- fun(x, method = method, MARGIN = 2)
    }
    }
    if (!is.null(row)) {
        for (i in 1:length(row)){
        method <- match.arg(row[i], METHODS)
        x <- fun(x, method = method, MARGIN = 1)
    }
}
x
}
