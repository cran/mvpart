"cmds.diss" <-
function (data, k = ncol(data), x.use = FALSE, zero.chk = TRUE,
     plt = FALSE, plot.subset = FALSE , plot.subn = 5, ...) 
{
    if (x.use) {
        xdists <- xdiss(data, ...)
        xds <- cmdscale(xdists, k = k)
        colnames(xds) <- paste("s",1:ncol(xds),sep="")
            }
    else {
        xdists <- gdist(data, ...)
        xds <- cmdscale(xdists, k = k)
        colnames(xds) <- paste("s",1:ncol(xds),sep="")

           }
    if (zero.chk) {
        drop.cols <- apply(xds, 2, function(x) (all(is.nan(x)) || 
            all(x == 0) || all(is.na(x))))
        if (any(drop.cols)) {
            cat(sum(drop.cols), " columns with NAs or all zeros dropped \n")
            xds <- xds[, !drop.cols]
        }
    }
    if (plt) {
        n <- nrow(data)
        if (n < 30 || !plot.subset) 
        plot(xdists, dxds <- dist(xds), xlim = c(0, max(xdists)), 
            ylim = c(0, max(dxds)), xlab = "Dists", ylab = "CMDS Dists", pch=1)
        else { 
        samp <- sample(n*n, floor(750 + n * plot.subn))
        dxds <- dist(xds)
        plot(xdists[samp], dxds[samp], xlim = c(0, md <- max(xdists)), 
            ylim = c(0, max(dxds)), xlab = "Dists", ylab = "CMDS Dists", pch=1)
        }    
        abline(c(0, 1), col = 2, xpd = FALSE)
        mtext("Pairwise distances vs CMD scaled pairwise distances", 
            3, line = 1.5)
        mtext(paste("R2 =", signif(cor(xdists, dxds)^2, 4), sep = ""), 
            3, line = -1.5)
        locator(1)
    }
    xds
}
