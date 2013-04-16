"sub.barplot" <-
function (x, y, z, keep = rep(TRUE , length(x)), row.scale = FALSE , xadj = 1,
    yadj = 1, bord = TRUE , line = TRUE , col = col)
{
    par(xpd = TRUE )
    drawbar <- function(x, y, z, line.adj = 0, xwid, ywid, border = TRUE ,
        line = TRUE , colbar = 10:12) {
        xx <- c(x - xwid/2, x + xwid/2, x + xwid/2, x - xwid/2)
        yy <- c(y - ywid, y - ywid, y, y)
        nbar <- length(z)
        xbwid <- xwid/nbar
        for (i in 1:nbar) {
            xb <- x - xwid/2 + xbwid * c(i - 1, i, i, i - 1)
            yb <- y - ywid + z[i] * ywid * c(0, 0, 1, 1) - line.adj *
                ywid
            polygon(xb, yb, col = colbar[i])
        }
        if (border)
            polygon(xx, yy, col = 1, density = 0)
        if (line)
            lines(xx[1:2], yy[1:2] - line.adj * ywid)
    }

    xrnge <- range(x)
    yrnge <- range(y)
    n <- length(x)
    xdiv <- max(sum(keep), 6)
    xwid <- (xadj * diff(xrnge))/xdiv
    ywid <- (yadj * diff(yrnge))/12
    x <- x[keep]
    y <- y[keep]
    z <- z[keep, ]
    nkeep <- sum(keep)
    nz <- ncol(z)
    if (length(col) < nz)
        col <- rep(col, nz, length = nz)
    if (any(z < 0)) {
        z <- z/diff(range(z))
        ladj <- min(z)
    }
    else {
        if (row.scale)
            z <- z/apply(z, 1, sum)
        else z <- z/max(z)
        ladj <- 0
    }
    for (i in 1:nkeep) drawbar(x[i], y[i], z[i, ], line.adj = ladj,
        xwid, ywid, border = bord, line = line, colbar = col)
}
