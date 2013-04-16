"plot.rpart" <- function (x, uniform = FALSE, branch = 1, compress = FALSE, nspace, 
    margin = 0.0, minbranch = 0.3, bar = 0.03, ...) 
{
    if (!inherits(x, "rpart")) 
        stop("Not an rpart object")
    if (!is.null(x$frame$splits)) 
        x <- rpconvert(x)
    if (compress & missing(nspace)) 
        nspace <- branch
    if (!compress) 
        nspace <- -1
    dev <- dev.cur()

##  added 15/04/13 
    parms <- list(uniform = uniform, branch = branch, nspace = nspace,
                 minbranch = minbranch)
    temp <- rpartco(x,parms)
    xx <- temp$x
    yy <- temp$y
    temp1 <- range(xx) + diff(range(xx)) * c(-margin, margin)
    temp2 <- range(yy) + diff(range(yy)) * c(-margin, margin)
    plot(temp1, temp2, type = "n", axes = FALSE, xlab = "", ylab = "", 
        ...)
    node <- as.numeric(row.names(x$frame))
    temp <- rpart.branch(xx, yy, node, branch)
    if (branch > 0) 
        lines(c(xx[1], xx[1]), c(yy[1], yy[1] + bar*diff(range(yy))), ...)
    lines(c(temp$x), c(temp$y))
    invisible(list(x = xx, y = yy))
}

