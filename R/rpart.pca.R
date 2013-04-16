"rpart.pca" <-
function (tree, pts = TRUE, plt.allx = TRUE, speclabs = TRUE, 
    specvecs = TRUE, wgt.ave = FALSE, add.tree = TRUE, cv1 = 1, 
    cv2 = 2, chulls = TRUE, interact = FALSE, ...) 
{
    if (tree$method != "mrt") 
        stop("Only for multivariate trees !! \n")
    if (nrow(tree$frame) < 4) 
        stop("Only 2 terminal nodes -- PCA not done !! \n")
    old.par <- par(mar = rep(2, 4), xpd = TRUE)
    on.exit(par(old.par))
    frame <- tree$frame
    ncf <- ncol(frame)
    data <- tree$y
    ny <- ncol(data)
    treegrps <- tree$where
    specs <- dimnames(data)[[2]]
    leaves <- frame$var == "<leaf>"
    n <- length(leaves)
    ln <- sum(leaves)
    lnot <- sum(!leaves)
    key <- dimnames(frame)[[1]]
    node <- as.numeric(key)
    even.node <- node[even <- node%%2 == 0]
    num <- length(specs)
    node.means <- as.matrix(frame[, ncf])
    tnode.means <- node.means[leaves, ]
    dimnames(node.means) <- list(key, specs)
    mat <- amat <- node.means - node.means[rep(1, n), ]
    mat <- mat[leaves, ]
    temp <- mat[rep(1:ln, frame[leaves, 2]), ]
    z <- svd(temp)
    maxd <- sum(z$d > 1e-06)
    d <- diag(z$d[1:maxd])
    xall <- z$u[, 1:maxd, drop = FALSE] %*% d
    x <- amat %*% (z$v)[, 1:maxd, drop = FALSE]
    xlv <- x[leaves, ]
    if (!wgt.ave) 
        y <- z$v[, 1:maxd, drop = FALSE]
    else {
        specvecs <- FALSE
        rc <- apply(tnode.means * frame$n[leaves], 2, sum)
        wgt <- diag(1/rc) %*% t(tnode.means * frame$n[leaves])
        y <- wgt %*% xlv
    }
    label <- 4:(3 + num)
    dstat <- signif(frame[leaves, "yval2"], digits = options()$digits)
    ln <- dim(dstat)[1]
    stat <- vector("character", length = ln)
    for (i in 1:ln) stat[i] <- paste(dstat[i, ], collapse = ", ")
    ymax <- max(dstat)
    ymin <- min(0, min(dstat))
    treegrps <- as.numeric(factor(treegrps))
    xx <- (scale(as.matrix(data), center = TRUE, scale = FALSE) %*% 
        z$v)[, 1:maxd, drop = FALSE]
    xrb <- rbind(x, xx)
    if (plt.allx) {
        mxx <- sqrt(apply(xrb[, c(cv1, cv2)]^2, 1, sum))
    }
    else mxx <- sqrt(apply(x[, c(cv1, cv2)]^2, 1, sum))
    cvar <- round((100 * z$d[1:maxd]^2)/sum(z$d[1:maxd]^2), digits = 2)
    cvar2 <- round(diag(cor(xall, xx[order(tree$where), ]))[1:maxd], 
        3)
    dlabs <- paste("   Dim ", c(1:maxd), " ", cvar, "% : [", 
        cvar2, "]")
    myy <- sqrt(apply(y[, c(cv1, cv2)]^2, 1, sum))
    sc <- ifelse(wgt.ave, 1, max(mxx)/max(myy))
    repeat {
        plot(c(sc * y[, cv1], xx[, cv1]), c(sc * y[, cv2], xx[, 
            cv2]), axes = FALSE, xlab = "", ylab = "", type = "n", 
            asp = 1)
        cxy <- par("cxy")
        sze <- par()$fin/par()$din
        adj <- ifelse(pts, cxy[2] * sze[2], 0)
        if (specvecs) 
            segments(sc * y[, cv1], sc * y[, cv2], rep(0, nrow(y)), 
                rep(0, nrow(y)), col = "gray", lty = 1)
        mtext(dlabs[cv1], side = 1, las = 0, adj = 0, line = 0, 
            cex = 0.85 * par()$cex)
        mtext(dlabs[cv2], side = 2, las = 0, adj = 0, line = 0, 
            cex = 0.85 * par()$cex)
        if (add.tree) {
            pp <- match(c(even.node, even.node + 1), node)
            nn <- length(even.node)
            from <- pp[1:nn]
            to <- pp[(nn + 1):(2 * nn)]
            segments(x[from, cv1], x[from, cv2], x[to, cv1], 
                x[to, cv2])
        }
        if (chulls) {
            unitg <- sort(unique(treegrps))
            for (i in 1:length(unitg)) {
                hpts <- chull(xx[unitg[i] == treegrps,c(cv1,cv2)])
                hpts <- c(hpts,hpts[1])
                lines(xx[unitg[i] == treegrps, c(cv1,cv2)][hpts,], col = i + 1)        
                }
        }
        if (plt.allx) {
            unitg <- sort(unique(treegrps))
            for (i in 1:length(unitg)) 
                points(xx[unitg[i] == treegrps, cv1], xx[unitg[i] == treegrps, cv2], 
                pch = 21, col = 1, bg = i + 1, cex = 1.2*par()$cex)
        }
        if (pts) {
            lvnode <- sort(node[leaves])
            for (i in 1:length(lvnode)) 
                points(xlv[, cv1][lvnode[i] == lvnode], xlv[, cv2][lvnode[i] == lvnode],
                pch = 21, cex = 2 * par()$cex, col = 1, bg = i + 1)
        }

        if (speclabs) 
            text(sc * y[, cv1], sc * y[, cv2] + 0.5 * adj * specvecs * 
                (y[, cv2] > 0), specs, col = "black", cex = par()$cex)
        points(0, 0, pch = 3, cex = par()$cex * 2.5, col = 1)
        if (interact) {
            z <- locator(1)
            if (length(z$x)) {
                if (z$x > 0 & z$y < 0) 
                  if (cv1 < maxd) 
                    cv1 <- cv1 + 1
                  else cv1 <- 1
                else if (z$x < 0 & z$y > 0) 
                  if (cv2 < maxd) 
                    cv2 <- cv2 + 1
                  else cv2 <- 2
                else if (z$x < 0 & z$y < 0) {
                  cv1 <- 1
                  cv2 <- 2
                }
            }
            else (break)
        }
        else (break)
    }
    invisible(list(y = sc * y, xlv = xlv, xx = xx))
}
