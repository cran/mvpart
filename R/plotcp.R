"plotcp" <-
function (x, xvse = 1, minline = TRUE , lty = 3, col = 1, upper = c("size",
    "splits", "none"), tab, resub.err = TRUE , adj.df = FALSE , ...)
{
    if (!inherits(x, "rpart"))
        stop("Not legitimate rpart object")
    upper <- match.arg(upper)
    p.rpart <- x$cptable
    if (xv <- (ncol(p.rpart) > 3))
    {
    xstd <- p.rpart[, 5]
    xerror <- p.rpart[, 4]
    }
    error <- p.rpart[, 3]
    nsplit <- p.rpart[, 2]
    ns <- seq(along = nsplit)
    cp0 <- p.rpart[, 1]
    cp <- sqrt(cp0 * c(Inf, cp0[-length(cp0)]))
    if (xv) {
        ylo <- min(c(xerror - xstd, error)) - 0.05
        yhi <- max(c(xerror + xstd, error)) + 0.05
    }
    else {
        ylo <- min(error) - 0.05
        yhi <- max(error) + 0.05
    }
    ylim <- c(ylo, yhi)
    plot(ns, error, axes = FALSE , xlab = "cp", ylab = "X-val Relative Error",
        ylim = ylim, type = "n", ...)
    if (xv) {
    inpt <- (xerror == min(xerror))
    points(ns[inpt], xerror[inpt], col = "red", pch = 16, cex = 2)
    inpt <- min(ns[xerror < min(xerror + xvse * xstd)])
    points(ns[inpt], xerror[inpt], col = "orange", pch = 16,
        cex = 2)
    points(ns, xerror, type = "b", col = "blue", ...)
    segments(ns, xerror - xstd, ns, xerror + xstd, col = "blue",
        ...)
    }
    if (resub.err)
        points(ns, error, type = "b", lty = 1, col = "darkgreen",
            ...)
    box()
    axis(2, ...)
    axis(1, at = ns, labels = as.character(signif(cp, 2)), ...)
    if (!missing(tab)) {
        xp <- as.numeric(names(tab))
        segments(ns[match(xp, nsplit + 1)], yhi, ns[match(xp,
            nsplit + 1)], yhi - 0.5 * (tab/sum(tab)) * (yhi -
            ylo), col = col + 1, lwd = 2, ...)
    }
    switch(upper, size = {
        axis(3, at = ns, labels = as.character(nsplit + 1), ...)
        mtext("Size of tree", side = 3, line = 3, cex = par()$cex,
            ...)
    }, splits = {
        axis(3, at = ns, labels = as.character(nsplit), ...)
        mtext("Number of splits", side = 3, line = 3, ...)
    }, )
    if (xv) {
    minpos <- min(seq(along = xerror)[xerror == min(xerror)])
    if (minline) {
        abline(h = (xerror + xvse * xstd)[minpos], lty = 1, col = col,
            xpd = FALSE )
        text(ns[2], (xerror + 0.5 * xvse * xstd)[minpos], paste("Min +",
            xvse, "SE"), col = col, ...)
    }
    }
    invisible()
}

