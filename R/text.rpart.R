"text.rpart" <-
function (x, splits = TRUE, which = 4, label = "yval", FUN = text, 
    all.leaves = FALSE, pretty = NULL, digits = getOption("digits") - 2,
    tadj = 0.65, stats = TRUE, use.n = FALSE, bars = TRUE, 
    legend = FALSE, xadj = 1, yadj = 1, bord = FALSE, big.pts = FALSE,
    uniform = FALSE, branch = 1, nspace = -1, minbranch = 0.3, ...) 
{
    if (!inherits(x, "rpart")) 
        stop("Not legitimate rpart")
   if (!is.null(x$frame$splits)) 
        x <- rpconvert(x)
    frame <- x$frame
    col <- names(frame)
    method <- x$method
    ylevels <- attr(x, "ylevels")
    if (!is.null(ylevels <- attr(x, "ylevels"))) 
        col <- c(col, ylevels)
    if (is.na(match(label, col))) 
        stop("Label must be a column label of the frame component of the tree")
    cxy <- par("cxy")
    if (!is.null(srt <- list(...)$srt) && srt == 90) 
        cxy <- rev(cxy)
    parms <- list(uniform = uniform, branch = branch, nspace = nspace,
                 minbranch = minbranch)
    xy <- rpartco(x,parms)
    node <- as.numeric(row.names(x$frame))
    is.left <- (node%%2 == 0)
    node.left <- node[is.left]
    parent <- match(node.left/2, node)
    bars <- bars & is.matrix(frame$yval2)
    text.adj <- ifelse(bars, yadj * diff(range(xy$y))/12, 0)
    if (splits) {
        left.child <- match(2 * node, node)
        right.child <- match(node * 2 + 1, node)
        rows <- labels(x, pretty = pretty)
        if (which == 1) 
            FUN(xy$x, xy$y + tadj * cxy[2], rows[left.child], 
                ...)
        else {
            if (which == 2 | which == 4) 
                FUN(xy$x, xy$y + tadj * cxy[2], rows[left.child], 
                  pos = 2, ...)
            if (which == 3 | which == 4) 
                FUN(xy$x, xy$y + tadj * cxy[2], rows[right.child], 
                  pos = 4, ...)
        }
    }
    leaves <- if (all.leaves) 
        rep(TRUE, nrow(frame))
    else frame$var == "<leaf>"
    if (stats) {
        if (is.null(frame$yval2)) 
            stat <- x$functions$text(yval = frame$yval[leaves], 
                dev = frame$dev[leaves], wt = frame$wt[leaves], 
                ylevel = ylevels, digits = digits, n = frame$n[leaves], 
                use.n = use.n)
        else stat <- x$functions$text(yval = frame$yval2[leaves, 
            ], dev = frame$dev[leaves], wt = frame$wt[leaves], 
            ylevel = ylevels, digits = digits, n = frame$n[leaves], 
            use.n = use.n)
        FUN(xy$x[leaves], xy$y[leaves] - tadj * cxy[2] - text.adj, 
            stat, adj = 0.5, ...)
    }
    if (bars) {
        bar.vals <- x$functions$bar(yval2 = frame$yval2)
        sub.barplot(xy$x, xy$y, bar.vals, leaves, xadj = xadj, 
            yadj = yadj, bord = bord, line = TRUE, col = c("lightblue", 
                "blue", "darkblue"))
        rx <- range(xy$x)
        ry <- range(xy$y)
        if (!is.null(ylevels)) 
            bar.labs <- ylevels
        else bar.labs <- dimnames(x$y)[[2]]
     if (legend & !is.null(bar.labs)) 
            legend(min(xy$x) - 0.1 * rx, max(xy$y) + 0.05 * ry, bar.labs, 
            col = c("lightblue", "blue", "darkblue"), pch = 15, bty = "n", ...)
    }
    if (big.pts) 
        points(xy$x[leaves], xy$y[leaves], pch = 16, cex = 3 * 
            par()$cex, col = 2:(sum(leaves) + 1))
    invisible()
}
