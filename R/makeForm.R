"makeform" <-
function (data, ycol, xcol, zcol, FUN, maxy = 20, extra) 
{
    ny <- length(ycol)
    nx <- ifelse(!missing(xcol), length(xcol), 0)
    nz <- ifelse(!missing(zcol), length(zcol), 0)
    dnames <- colnames(data)
    if (ny > maxy) {
        yy <- deparse(substitute(ycol))
        ty <- paste("as.matrix(", deparse(substitute(data)), 
            "[,", yy, "])", collapse = "", sep = "")
        if (!missing(FUN)) 
            ty <- paste(deparse(substitute(FUN)), "(", ty, ")", 
                sep = "")
        if (nx > 1) 
            tx <- paste(dnames[xcol], collapse = "+")
        else if (nx == 1) 
            tx <- dnames[xcol]
        else tx <- "1"
    }
    else {
        if (ny > 1) {
            ty <- paste(dnames[ycol], collapse = ",")
            ty <- paste("cbind(", ty, ")", collapse = "", sep = "")
        }
        else if (ny == 1) 
            ty <- dnames[ycol]
        else ty <- ""
        if (!missing(FUN)) 
            ty <- paste(deparse(substitute(FUN)), "(", ty, ")", 
                sep = "")
        if (nx > 1) 
            tx <- paste(dnames[xcol], collapse = "+")
        else if (nx == 1) 
            tx <- dnames[xcol]
        else tx <- "1"
    }
    if (missing(extra)) 
            form <- paste(ty, "~", tx, collapse = "", sep = "")
    else form <- paste(ty, "~", extra, tx, collapse = "", sep = "")

    if (nz > 1) 
        tz <- paste(dnames[zcol], collapse = "+")
    else if (nz == 1) 
        tz <- dnames[zcol]

    if (nz > 0) 
        form <- paste(form, "+ Condition(", tz, ")", sep = "")

    formula(form)
}
