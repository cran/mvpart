"trclcomp" <- function (x, method = "com", km = TRUE, mrt = TRUE)
{
    if (class(x) != "rpart")
        stop("Rpart object needed!")
    if (x$method != "mrt")
        stop("Multivariate tree needed!")
    pruner <- function(x) {
        cps <- x$cptable[,1]
        nps <- length(cps)
        groups <- matrix(0,nrow=length(x$where),ncol=length(cps))
        groups[,1] <- x$where
        for (i in 1:nps) {
        new.x <- prune.rpart(x, cp=cps[i])
        cat(length(unique(new.x$where)),"")
        groups[,i] <- new.x$where
    }
    cat("\n")
    groups
    }
    cpt <- x$cptable
    size <- cpt[, 2] + 1
    nr <- nrow(cpt)
    mrt.err <- cpt[, 3]
    mrt.clust.err <- clust.err <- rep(1, nr)
    sst <- sum(scale(x$y, scale = FALSE)^2)
    n <- nrow(x$y)
    d <- dist(x$y)
    if (any(is.na(d))) {
        cat("Warning -- NA distances in cluster -- replacing by 0\n")
        d[is.na(d)] <- 0
    }
    hclout <- hclust(d, method = method)
    grp.mrt <- pruner(x)
    for (i in 2:nr) {
        grp.clust <- factor(cutree(hclout, k = size[i]))
        cents <- t(sapply(split(as.data.frame(x$y), grp.clust), colMeans))
        grp.clust <- factor(kmeans(x$y, centers = cents)$cluster)
        cents.mrt <- t(sapply(split(as.data.frame(x$y), grp.mrt[, i]), colMeans))
        grp.mrt.clust <- factor(kmeans(x$y, centers = cents.mrt)$cluster)
        clust.err[i] <- sum(resid(lm(x$y ~ factor(grp.clust), singular.ok = TRUE))^2)/sst
        mrt.clust.err[i] <- sum(resid(lm(x$y ~ factor(grp.mrt.clust), singular.ok = TRUE))^2)/sst
    }
    minerr <- min(c(mrt.err, mrt.clust.err, clust.err))
    plot(size, mrt.err, type = "n", ylim = c(minerr, 1), xlab = "Size",
        ylab = "Resubstition Error")
    points(size, mrt.err, type = "o", col = 2, pch = 16)
    points(size, mrt.clust.err, type = "o", col = 3, pch = 16)
    points(size, clust.err, type = "o", col = 4, pch = 16)
    legend(mean(size), 1, c("MRT", "MRT-Cluster", "Cluster"), col = c(2:4),
        lty = 1, bty = "n")
    title("Comparison of tree and cluster errors across size")
    cat("MRT error                : ", signif(mrt.err, 3), "\n")
    cat("MRT.Cluster error        : ", signif(mrt.clust.err, 3), "\n")
    cat("Cluster error            : ", signif(clust.err, 3), "\n")
    invisible(list(mrt.err = mrt.err, mrt.clust.err = mrt.clust.err, clust.err = clust.err))
}
