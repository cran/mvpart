"mvpart" <- function (form, data, minauto = TRUE, size, xv = c("1se", "min", 
    "pick", "none"), xval = 10, xvmult = 0, xvse = 1, snip = FALSE, 
    plot.add = TRUE, text.add = TRUE, digits = 3, margin = 0, 
    uniform = FALSE, which = 4, pretty = TRUE, use.n = TRUE, 
    all.leaves = FALSE, bars = TRUE, legend, bord = FALSE, 
    xadj = 1, yadj = 1, prn = FALSE, branch = 1, rsq = FALSE, 
    big.pts = FALSE, pca = FALSE, interact.pca = FALSE, 
    wgt.ave.pca = FALSE, keep.y = TRUE, ...) 
{
    call <- match.call()

	number <-
		function (x) {
		    match(x, sort(unique(x)))
	}
    cv.var <- function(x, cv = 10)  {
        x <- match(x, sort(unique(x)))
        luni <- length(unique(x))
        if(luni >= cv) {
        grps <- ceiling((cv * cumsum(table(x)))/length(x))
        x <- number(grps[x])
        }
        x
    }
    if (length(xval) > 1) {
        if (xvmult > 1) 
            xvalvar <- xval
        xval <- cv.var(xval)
    }
    choice <- c("1se", "min", "pick", "none")
    xv <- choice[pmatch(xv[1], choice)]
    if (!missing(size) || xv == "none") 
        xval <- 0
    if (minauto) {
        n <- nrow(data)
        minsplit <- ceiling(log2(n))
        minbucket <- ceiling(minsplit/3)
    }
    z <- rpart(form, data = data, ...)
    if (all(z$where==1)) {
    cat("No splits possible -- try decreasing cp\n")
    return(z)
    }
    old.par <- par(mar = c(6, 4, 4, 4) + 0.1, xpd = NA, cex = par()$cex)
    on.exit(par(old.par))

    if (!is.null(z)) {
        xval <- z$control$xval
        if (xvmult > 1) {
            zresse <- zres <- matrix(NA, nrow = nrow(z$cptable), 
                ncol = xvmult)
            zres[, 1] <- z$cptable[, 4]
            zresse[, 1] <- z$cptable[, 5]
            cat("X-Val rep : 1")
            for (i in 2:xvmult) {
                if (length(xval) == nrow(data)) 
                  xval <- cv.var(xvalvar)
                ztemp <- rpart(form, data = data, ...)$cptable[, 
                  4:5]
                zres[, i] <- ztemp[, 1]
                zresse[, i] <- ztemp[, 2]
                cat(" ", i)
                NULL
            }
            cat("\n")
            z$cptable[, 4] <- apply(zres, 1, mean)
            z$cptable[, 5] <- apply(zresse, 1, mean)
            tabmins <- apply(zres, 2, function(x, nc, sizes) {
                sizes[x == min(x)][1]
            }, nc = nrow(zres), sizes = z$cptable[, 2] + 1)
            cat("Minimum tree sizes\n")
            print(table(tabmins))
        }
        if (missing(size)) {
        if (xv == "pick") {
            if (xvmult <= 1) 
                plotcp(z, xvse, pch = 16, col = 2)
            else plotcp(z, xvse, pch = 16, col = 2, tab = table(tabmins))
            size.loc <- locator(1)
            if (!is.null(size.loc)) {
                splt <- round(size.loc$x)
                if (splt < 2) 
                  splt <- 2
                else if (splt > length(z$cptable[, 1])) 
                  splt <- length(z$cptable[, 1])
                cpp <- z$cptable[, 1][splt]
                z <- prune.rpart(z, cp = cpp)
            }
        }
        else if ((xv == "1se" | xv == "min") && (xval[1] != 0)) {
            xerror <- z$cptable[, 4]
            xstd <- z$cptable[, 5]
            if (xv == "min") 
                splt <- min(seq(along = xerror)[xerror == min(xerror)])
            else splt <- min(seq(along = xerror)[xerror <= min(xerror) + 
                xvse * xstd])
            if (!is.na(splt)) {
                if (splt == 1) 
                  splt <- 2
                cpp <- z$cptable[, 1][splt]
                z <- prune.rpart(z, cp = cpp)
            }
            else {
                (cat("No pruning possible : size 2 tree produced ?? \n"))
                use.size <- TRUE
                size <- 2
            }
            }
        }
        else {
            if (size <= 2) 
                cpp <- z$cptable[2, 1]
            else if (size >= max(z$cptable[, 2] + 1)) 
                cpp <- z$cptable[dim(z$cptable)[1], 1]
            else cpp <- z$cptable[, 1][min(abs(size - z$cptable[, 
                2] - 1)) == abs(size - z$cptable[, 2] - 1)][1]
            z <- prune.rpart(z, cp = cpp)
        }
        if (snip) {
            plot(z)
            z <- snip.rpart(z)
        }
        if (rsq && xval != 0 && z$method != "class") {
            par(mfrow = c(1, 2))
            rsq.rpart(z)
            locator(1)
            par(mfrow = c(1, 1))
        }
        if (plot.add) {
            plot.rpart(z, uniform = uniform, branch = branch, 
                margin = margin)
            if (text.add) 
                text.rpart(z, digits = digits, xadj = xadj, yadj = yadj, 
                  which = which, pretty = pretty, use.n = use.n, bars = bars,
                  legend = ifelse(missing(legend),(z$method=="mrt") 
                  && (ncol(z$frame$yval2)<20),legend),
                  all.leaves = all.leaves, bord = bord, big.pts = big.pts | pca)
            len <- dim(z$cptable)[1]
            foot <- paste("Error : ", signif(z$cptable[len, 3], 
                digits))
            if (xval[1] > 0) 
                foot <- paste(foot, "  CV Error : ", signif(z$cptable[len, 
                  4], digits), "  SE : ", signif(z$cptable[len, 
                  5], digits))
            mtext(foot, side = 1, line = 3.5, cex = par()$cex)
            n <- dim(data)[1]
            if (z$method == "class") {
                nex <- max(table(z$y))
                foot2 <- paste("Missclass rates : Null = ", signif(1 - 
                  nex/n, digits), " : Model = ", signif((1 - 
                  nex/n) * z$cptable[len, 3], digits))
                if (xval[1] > 0) 
                  foot2 <- paste(foot2, " : CV = ", signif((1 - 
                    nex/n) * z$cptable[len, 4], digits))
                mtext(foot2, side = 1, line = 4.5, cex = par()$cex)
            }
        }
        if (prn) 
            printcp(z)
        if (pca) {
            locator(1)
            rpart.pca(z, interact = interact.pca, wgt.ave = wgt.ave.pca)
        }
    }
    else {
        plot(c(-1, 1), c(-1, 1), axes = FALSE, type = "n", xlab = "", 
            ylab = "")
        text(c(0), c(0), "No splits could be formed", col = 2)
        cat("No splits could be formed\n")
    }
    if (!is.null(z)) {
        if (!keep.y) z$y <- NULL
        z$call <- call
        invisible(z)
    }
}
