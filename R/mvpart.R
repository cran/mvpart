.packageName <- "mvpart"

# SCCS  @(#)formatg.s	1.3 06/06/01
# format a set of numbers using C's "g" format
#  It is applied on an element by element basis, which is more
#  appropriate for rpart output than the standard Splus format()
#  command.
# For instance if x=(123, 1.23, .00123)
#	  format(x) = "123.00000", "1.23000", "0.00123"
#  but formatg does not add all of those zeros to the first two numbers
#
formatg <- function(x, digits= unlist(options('digits')),
		         format= paste("%.", digits, "g", sep='')) {
    if (!is.numeric(x)) stop("x must be a numeric vector")

    n <- length(x)
    #
    # the resultant strings could be up to 8 characters longer,
    #   assume that digits =4,  -0.dddde+104 is a worst case, where
    #   dddd are the 4 significant digits.
    dummy  <- paste(rep(" ", digits+8), collapse='')
    temp <- .C("formatg", as.integer(n),
	                  as.double(x),
                          rep(format,n),
                          out= rep(dummy, n), NAOK=TRUE, PACKAGE="mvpart")$out
    if (is.matrix(x)) matrix(temp, nrow=nrow(x))
    else temp
    }

sub.barplot <- function (x, y, z, keep = rep(TRUE , length(x)), row.scale = FALSE , xadj = 1, 
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
            polygon(xx, yy, col = 1, dens = 0)
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


eqscplt <- function(x, y, tol = 0.15, sym = FALSE , ...)
{
        if(is.matrix(x)) {
                y <- x[, 2]
                x <- x[, 1]
        }
        if(is.list(x)) {
                y <- x$y
                x <- x$x
        }
        oldpin <- par("pin")
        if(!sym) {
                xlim <- range(x, na.rm = TRUE )
                midx <- 0.5 * (xlim[2] + xlim[1])
                xlim <- midx + (1 + tol) * 0.5 * c(-1, 1) * (xlim[2] - xlim[1])
                ylim <- range(y, na.rm = TRUE )
                midy <- 0.5 * (ylim[2] + ylim[1])
                ylim <- midy + (1 + tol) * 0.5 * c(-1, 1) * (ylim[2] - ylim[1])
        }
        else {
                xlim <- c(-1, 1) * max(abs(x), na.rm = TRUE ) * (1 + tol)
                ylim <- c(-1, 1) * max(abs(y), na.rm = TRUE ) * (1 + tol)
                midx <- midy <- 0
        }
        xr <- oldpin[1]/(xlim[2] - xlim[1])
        yr <- oldpin[2]/(ylim[2] - ylim[1])
        if(yr > xr) {
                ylim <- midy + (yr * c(-1, 1) * (ylim[2] - ylim[1]))/(2 * xr)
        }
        else {
                xlim <- midx + (xr * c(-1, 1) * (xlim[2] - xlim[1]))/(2 * yr)
        }
        plot(x, y, xlim = xlim, ylim = ylim, ...)
        list(xlim = xlim, ylim = ylim)
}

mvpart <- function (form, data, minauto = TRUE, size, xv = c("1se", "min", 
    "pick", "none"), xval = 10, xvmult = 0, xvse = 1, snip = FALSE, 
    plot.add = TRUE , text.add = TRUE , digits = 3, margin = 0, uniform = FALSE, 
    which = 1, pretty = TRUE, use.n = TRUE, all = FALSE, bord = FALSE, 
    xadj = 1, yadj = 1, prn = FALSE, branch = 1, rsq = FALSE, 
    big.pts = FALSE, pca = FALSE, interact.pca = FALSE, wgt.ave.pca = FALSE, 
    ...) 
{
    call <- match.call()
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
    use.size <- FALSE
    z <- rpart(form, data = data, xval = xval, ...)
    old.par <- par(mar = c(6, 4, 4, 4) + 0.1)
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
        if (use.size || !missing(size)) {
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
                  which = which, pretty = pretty, use.n = use.n, all = all, 
                  bord = bord, big.pts = big.pts | pca)
            len <- dim(z$cptable)[1]
            foot <- paste("Error : ", signif(z$cptable[len, 3], 
                digits))
            if (xval[1] > 0) 
                foot <- paste(foot, "  CV Error : ", signif(z$cptable[len, 
                  4], digits), "  SE : ", signif(z$cptable[len, 
                  5], digits))
            mtext(foot, side = 1, line = 3.5, cex = par()$cex)
            n <- dim(data)[1]
            nex <- max(table(z$y))
            if (z$method == "class") {
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
        plot(c(-1, 1), c(-1, 1), axes = FALSE , type = "n", xlab = "", 
            ylab = "")
        text(c(0), c(0), "No splits could be formed", col = 2)
        cat("No splits could be formed\n")
    }
    if (!is.null(z)) {
        z$call <- call
        invisible(z)
    }
}

# SCCS @(#)labels.rpart.s	1.5 07/17/01
# Make the nice labels used by print and summary
#   digits = obvious
#   minlength = 0 = don't abbrev factors
#               1 = use single letters
#               2+= the same arg as the "abbreviate" function
#   collapse = an oddly named argument
#              FALSE = return a matrix with two columns, containing the labels of
#                    the left and right descendants of each node
#              TRUE = return a vector of 1 column, with the label of the parent
#   pretty: for historical compatability
#               0   -> minlength =0
#              NULL -> minlength =1
#               TRUE   -> minlength =4
#   ... = other args for abbreviate()
#
labels.rpart <- function(object, digits=4, minlength=1, pretty,
			      collapse=TRUE, ...) {
    if (missing(minlength) && !missing(pretty)) {
	if (is.null(pretty)) minlength <-1
	else if (is.logical(pretty)) {
	    if (pretty) minlength <- 4
	    else        minlength <- 0
	    }
	else minlength <- 0
	}

    ff <- object$frame
    n  <- nrow(ff)
    if (n==1) return("root")  #special case of no splits

    is.leaf <- (ff$var == "<leaf>")
    whichrow <- !is.leaf
    vnames <- ff$var[whichrow]  #the variable names for the primary splits

    index <- cumsum(c(1, ff$ncompete + ff$nsurrogate + 1*(!is.leaf)))
    irow  <- index[c(whichrow, FALSE)]     #we only care about the primary split
    ncat  <- object$splits[irow, 2]

    # Now to work: first create labels for the left and right splits,
    #  but not for leaves of course
    #
    lsplit <- rsplit <- vector(mode='character', length= length(irow))

    if (any(ncat <2)) {  # any continuous vars ?
	jrow <- irow[ncat <2]
	cutpoint <- formatg(object$splits[jrow,4], digits)
	temp1 <- (ifelse(ncat<0, "< ", ">="))[ncat <2]
	temp2 <- (ifelse(ncat<0, ">=", "< "))[ncat <2]
	lsplit[ncat<2] <- paste(temp1, cutpoint, sep='')
	rsplit[ncat<2] <- paste(temp2, cutpoint, sep='')
	}

    if (any(ncat >1)) { # any categorical variables ?
	xlevels <- attr(object, 'xlevels')
	#
	# jrow will be the row numbers of factors within lsplit and rsplit
	# crow the row number in "csplit"
	# and cindex the index on the "xlevels" list
	#
	jrow <- (seq(along=ncat))[ncat>1]
	crow <- object$splits[irow[ncat>1],4]    #row number in csplit
	cindex <- (match(vnames, names(xlevels)))[ncat >1]

	# Now, abbreviate the levels
	if (minlength ==1) {
	    if (any(ncat>52))
		warning(paste("More than 52 levels in a predicting factor,",
			      "truncated for printout"))
	    xlevels <- lapply(xlevels,
			       function(z) {
				   k <- length(z)
				   k <- pmin(1:k, 52)
				   c(letters, LETTERS)[k]
				   })
	    }
	else if (minlength >1)
	    xlevels <- lapply(xlevels, abbreviate, minlength, ...)

	# Now tuck in the labels
	# I'll let some other clever person vectorize this
	for (i in 1:(length(jrow))) {
	    j <- jrow[i]
	    splits <- object$csplit[crow[i],]
	    # splits will contain 1=left, 2=right, 3= neither
	    ltemp <- (1:length(splits))[splits== 1]
	    rtemp <- (1:length(splits))[splits== 3]
	    if (minlength==1) {
		lsplit[j] <- paste((xlevels[[cindex[i]]])[ltemp], collapse='')
		rsplit[j] <- paste((xlevels[[cindex[i]]])[rtemp], collapse='')
		}
	    else {
		lsplit[j] <-paste((xlevels[[cindex[i]]])[ltemp], collapse=',')
		rsplit[j] <-paste((xlevels[[cindex[i]]])[rtemp], collapse=',')
		}
	    }
	}

    if (!collapse) {  #called by no routines that I know of
	ltemp <- rtemp <- rep("<leaf>", n)
	ltemp[whichrow] <- lsplit
	rtemp[whichrow] <- rsplit
	return(cbind(ltemp, rtemp))
	}

    lsplit <- paste(ifelse(ncat<2, "", "="), lsplit, sep='')
    rsplit <- paste(ifelse(ncat<2, "", "="), rsplit, sep='')

    #
    # Now match them up to node numbers
    #   The output will have one label per row of object$frame, each
    #   corresponding the the line segement joining this node to its parent
    #
    varname <- (as.character(vnames))
    node <- as.numeric(row.names(ff))
    parent <- match(node %/% 2, node[whichrow])
    odd <- (as.logical(node %%2))

    labels <- vector('character', length=n)
    labels[odd] <- paste(varname[parent[odd]], rsplit[parent[odd]], sep="")
    labels[!odd]<- paste(varname[parent[!odd]],lsplit[parent[!odd]], sep="")
    labels[1] <- "root"
    labels
    }
# SCCS 02/18/97 @(#)meanvar.rpart.s	1.2

meanvar.rpart <- function(tree, xlab = "ave(y)", ylab = "ave(deviance)", ...)

{
	if(!inherits(tree, "rpart"))
		stop("Not legitimate rpart object")
	if(!tree$method=='anova')
		stop("Plot not useful for classification or poisson trees")
	frame <- tree$frame
	frame <- frame[frame$var == "<leaf>",  ]
	x <- frame$yval
	y <- frame$dev/frame$n
	label <- row.names(frame)
	plot(x, y, xlab = xlab, ylab = ylab, type = "n", ...)
	text(x, y, label)
	invisible(list(x = x, y = y, label = label))
}

meanvar <- function(tree,...) UseMethod('meanvar')
# sccs @(#)model.frame.rpart.s	1.3 01/21/97
model.frame.rpart <- function(formula, ...)
{
	m <- formula$model
	if(!is.null(m))
		return(m)
	oc <- formula$call
	if(substring(deparse(oc[[1]]), 1, 7) == "predict") {
		m <- eval(oc$newdata)
		if(is.null(attr(m, "terms"))) {
			object <- eval(oc$object)
			m <- model.frame(object$terms, m, na.rpart)
		}
		return(m)
	}
	while(deparse(oc[[1]]) != "rpart") oc <- eval(oc[[2]])$call
	oc$subset <- names(formula$where)
	oc$method <- formula$method
	eval(oc)
}

#SCCS  @(#)na.rpart.s	1.5 12/13/99
na.rpart <- function(x){
    Terms <- attr(x, 'terms')
    if(!is.null(Terms)) yvar <- attr(Terms, "response") else yvar <- 0
    if (yvar==0) {
	xmiss <- is.na(x)
	keep <-  (xmiss %*% rep(1,ncol(xmiss))) < ncol(xmiss)
	}
    else {
	xmiss <- is.na(x[-yvar])
	ymiss <- is.na(x[[yvar]])
	if (is.matrix(ymiss))
	    keep <- ((xmiss %*% rep(1,ncol(xmiss))) < ncol(xmiss)) &
		    ((ymiss %*% rep(1,ncol(ymiss))) == 0 )
	else
	    keep <- ((xmiss %*% rep(1,ncol(xmiss))) < ncol(xmiss)) & !ymiss
	}
    if (all(keep)) x
    else {
	temp <- seq(keep)[!keep]
	names(temp) <- row.names(x)[!keep]
	#the methods for this group are all the same as for na.omit
	class(temp) <- c("na.rpart", "omit")
	structure(x[keep,], na.action=temp)
	}
    }
## submitted by Anantha Prasad 1/26/98

path.rpart <- function(tree, nodes, pretty = 0, print.it = TRUE)
{
        if(!inherits(tree, "rpart"))
                stop("Not legitimate tree")
        splits <- labels.rpart(tree, pretty = pretty)
        frame <- tree$frame
        n <- row.names(frame)
        node <- as.numeric(n)
        which <- descendants(node)      #ancestors are columns
        path <- list()
        if(missing(nodes)) {
                xy <- rpartco(tree)
                while(length(i <- identify(xy, n = 1, plot = FALSE)) > 0) {
                        path[[n[i]]] <- path.i <- splits[which[, i]]
                        if(print.it) {
                                cat("\n", "node number:", n[i], "\n")
                                cat(paste("  ", path.i), sep = "\n")
                        }
                }
        }
        else {

                if(length(nodes <- node.match(nodes, node)) == 0)
                        return(invisible())
                for(i in nodes)
                       { path[[n[i]]] <- path.i <- splits[which[, i]]
			if(print.it) {
                                cat("\n", "node number:", n[i], "\n")
                                cat(paste("  ", path.i), sep = "\n")
                                }
                       }
        }
        invisible(path)
}





#SCCS @(#)plot.rpart.s	1.8 06/08/01

plot.rpart <- function (x, uniform = FALSE, branch = 1, compress = FALSE, nspace, 
    margin = 0.0, minbranch = 0.3, bar = 0.03, ms.fudge = TRUE,...) 
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
    if (dev == 1) 
        dev <- 2
    assign(paste(".rpart.parms", dev, sep = "."), list(uniform = uniform, 
        branch = branch, nspace = nspace, minbranch = minbranch), 
        envir = .GlobalEnv)
    temp <- rpartco(x)
    xx <- temp$x
    yy <- temp$y
    temp1 <- range(xx) + diff(range(xx)) * c(-margin, margin)
    temp2 <- range(yy) + diff(range(yy)) * c(-margin, margin)
    plot(temp1, temp2, type = "n", axes = FALSE, xlab = "", ylab = "", 
        ...)
    node <- as.numeric(row.names(x$frame))
    temp <- rpart.branch(xx, yy, node, branch)
#
#  fudge for MS PPT "closing boxes on ungroup"  GD 11/02
#
    if(ms.fudge)
        temp$y[4,  ] <- temp$y[4,  ] + 0.01*diff(range(temp$y[4,]))
    if (branch > 0) 
        lines(c(xx[1], xx[1]), c(yy[1], yy[1] + bar*diff(range(yy))), ...)
    lines(c(temp$x), c(temp$y))
    invisible(list(x = xx, y = yy))
}


# SCCS @(#)plotcp.s	1.1 02/08/98
# Contributed by B.D. Ripley 97/07/17
#

plotcp <- function (x, xvse = 1, minline = TRUE , lty = 3, col = 1, upper = c("size", 
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
    axis(1, at = ns, lab = as.character(signif(cp, 2)), ...)
    if (!missing(tab)) {
        xp <- as.numeric(names(tab))
        segments(ns[match(xp, nsplit + 1)], yhi, ns[match(xp, 
            nsplit + 1)], yhi - 0.5 * (tab/sum(tab)) * (yhi - 
            ylo), col = col + 1, lwd = 2, ...)
    }
    switch(upper, size = {
        axis(3, at = ns, lab = as.character(nsplit + 1), ...)
        mtext("Size of tree", side = 3, line = 3, cex = par()$cex, 
            ...)
    }, splits = {
        axis(3, at = ns, lab = as.character(nsplit), ...)
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


# SCCS 05/11/01 @(#)post.rpart.s	1.13
#
post.rpart <- function(tree, title.,
		       filename=paste(deparse(substitute(tree)),".ps",sep=""),
		       digits=getOption("digits") - 3, pretty=TRUE,
		       use.n=TRUE,  horizontal=TRUE, ...)
{
    if(filename !=""){
	postscript(file = filename, horizontal=horizontal, ...)
	par(mar=c(2,2,4,2)+.1)
	on.exit(dev.off())
	}
    else {
	oldpar <- par(mar=c(2,2,4,2)+.1)
	on.exit(invisible(par(oldpar)))
	}

    plot(tree, uniform=TRUE, branch=.2, compress=TRUE, margin=.1)
    text(tree, all=TRUE, use.n=use.n, digits=digits, pretty=pretty)
    method <- tree$method

    if(missing(title.)) {
        temp  <- attr(tree$terms,'variables')[2]
        title(paste("Endpoint =",temp),cex=.8)
    } else if (title. !="") title(title.,cex=.8)
}

## SCCS @(#)post.s	1.3 02/27/98
post <- function(tree, ...) UseMethod("post")

# SCCS @(#)pred.rpart.s	1.3 09/03/97
#
# Do Rpart predictions given a tree and a matrix of predictors
pred.rpart <- function(fit, x) {

    frame <- fit$frame
    if(nrow(frame) == 1) { # handle root-only tree separately
        temp <- rep(1, nrow(x))
    } else {
        nc <- frame[, c('ncompete', 'nsurrogate')]
        frame$index <- 1 + c(0, cumsum((frame$var != "<leaf>") +
                                       nc[[1]] + nc[[2]]))[-(nrow(frame)+1)]
        frame$index[frame$var == "<leaf>"] <- 0
        vnum <- match(dimnames(fit$split)[[1]], dimnames(x)[[2]])
        if (any(is.na(vnum)))
            stop("Tree has variables not found in new data")
        temp <- .C("pred_rpart",
                        as.integer(dim(x)),
                        as.integer(dim(frame)[1]),
                        as.integer(dim(fit$splits)),
                        as.integer(if(is.null(fit$csplit)) rep(0,2)
                                   else dim(fit$csplit)),
                        as.integer(row.names(frame)),
                        as.integer(unlist(frame[,
                               c('n', 'ncompete', 'nsurrogate', 'index')])),
                        as.integer(vnum),
                        as.double(fit$splits),
                        as.integer(fit$csplit -2),
                        as.integer((fit$control)$usesurrogate),
                        as.double(x),
                        as.integer(is.na(x)),
                        where = integer(dim(x)[1]),
                        NAOK = TRUE, PACKAGE = "mvpart")
        temp <- temp$where
    }
    names(temp) <- rownames(x)
    temp
}
## SCCS @(#)predict.rpart.s	1.11 06/03/01
predict.rpart <-
function(object, newdata = list(),
	 type = c("vector", "prob", "class", "matrix"), ...) {
    if(!inherits(object, "rpart"))
	stop("Not legitimate tree")
    mtype <- missing(type)
    type <- match.arg(type)
    if(missing(newdata))
	where <- object$where
    else {
	if(is.null(attr(newdata, "terms"))) {
	    Terms <- delete.response(object$terms)
	    act <- (object$call)$na.action
	    if (is.null(act)) act<- na.rpart
	    newdata <- model.frame(Terms, newdata, na.action = act,
                                      xlev=attr(object, "xlevels"))
        }
	where <- pred.rpart(object, rpart.matrix(newdata))
    }
    frame <- object$frame
    method <- object$method
    ylevels <- attr(object, "ylevels")
    nclass <- length(ylevels)
    if(mtype && nclass > 0) type <- "prob"
    if(mtype && method=="mrt") type <- "matrix"
    if(type == "vector" || (type=="matrix" && is.null(frame$yval2))) {
	pred <- frame$yval[where]
	names(pred) <- names(where)
    }
    else if (type == "matrix") {
	pred <- frame$yval2[where,]
	dimnames(pred) <- list(names(where), NULL)
    }
    else if(type == "class" && nclass > 0) {
	pred <- factor(ylevels[frame$yval[where]], levels=ylevels)
	names(pred) <- names(where)
    }
    else if (type == "prob" && nclass > 0) {
	pred <- frame$yval2[where, 1 + nclass + 1:nclass]
	dimnames(pred) <- list(names(where), ylevels)
    }
    else stop("Invalid prediction for rpart object")

    # Expand out the missing values in the result
    # But only if operating on the original dataset
    if (missing(newdata) && !is.null(object$na.action))
        pred <- naresid(object$na.action, pred)
    pred
}

#SCCS  @(#)print.rpart.s	1.15 06/06/01
print.rpart <- function(x, minlength=0, spaces=2, cp,
               digits=getOption("digits"), ...) {
    if(!inherits(x, "rpart")) stop("Not legitimate rpart object")
    if (!is.null(x$frame$splits)) x <- rpconvert(x)  #help for old objects

    if (!missing(cp)) x <- prune.rpart(x, cp=cp)
    frame <- x$frame
    ylevel <- attr(x, "ylevels")
    node <- as.numeric(row.names(frame))
    depth <- tree.depth(node)
    indent <- paste(rep(" ", spaces * 32), collapse = "")
    #32 is the maximal depth
    if(length(node) > 1) {
        indent <- substring(indent, 1, spaces * seq(depth))
        indent <- paste(c("", indent[depth]), format(node), ")", sep = "")
    }
    else indent <- paste(format(node), ")", sep = "")

    tfun <- (x$functions)$print
    if (!is.null(tfun)) {
	if (is.null(frame$yval2))
		yval <- tfun(frame$yval,  ylevel, digits)
	else    yval <- tfun(frame$yval2,  ylevel, digits)
	}
    else yval <- format(signif(frame$yval, digits = digits))
    term <- rep(" ", length(depth))
    term[frame$var == "<leaf>"] <- "*"
    z <- labels(x, digits=digits, minlength=minlength, ...)
    n <- frame$n
    z <- paste(indent, z, n, format(signif(frame$dev, digits = digits)),
               yval, term)

    omit <- x$na.action
    if (length(omit))
    cat("n=", n[1], " (", naprint(omit), ")\n\n", sep="")
    else cat("n=", n[1], "\n\n")

    #This is stolen, unabashedly, from print.tree
    if (x$method=="class")
         cat("node), split, n, loss, yval, (yprob)\n")
    else cat("node), split, n, deviance, yval\n")
    cat("      * denotes terminal node\n\n")

    cat(z, sep = "\n")
    return(invisible(x))
    #end of the theft
    }
#SCCS  @(#)printcp.s	1.6 01/20/00
# print out the cptable, along with some summary of the tree
printcp <- function(x, digits=getOption("digits")-2)
{
    if (!inherits(x, 'rpart')) stop ("Must be an rpart x")
    cat(switch(x$method,anova = "\nRegression tree:\n" ,
			class = "\nClassification tree:\n" ,
			poisson="\nRates regression tree:\n",
			exp = "\nSurvival regression tree:\n")
        )

    if(!is.null(cl <- x$call)) {
	dput(cl)
	cat("\n")
    }
    frame <- x$frame
    leaves <- frame$var == "<leaf>"
    used <- unique(frame$var[!leaves])

    if(!is.null(used)) {
        cat("Variables actually used in tree construction:\n")
        print(sort(as.character(used)), quote=FALSE)
        cat("\n")
    }


    cat("Root node error: ", format(frame$dev[1], digits=digits), '/',
        frame$n[1], ' = ',
        format(frame$dev[1]/frame$n[1], digits=digits),
        '\n\n', sep='')


    n <- x$frame$n
    omit <- x$na.action
    if (length(omit))
    cat("n=", n[1], " (", naprint(omit), ")\n\n", sep="")
    else cat("n=", n[1], "\n\n")

    print (x$cptable, digits=digits)
    invisible(x$cptable)
}


#SCCS @(#)prune.rpart.s	1.9 10/30/01
prune.rpart <- function(tree, cp, ...)
{
    ff <- tree$frame
    id <- as.integer(row.names(ff))
    toss <- id[ff$complexity <= cp &  ff$var!='<leaf>']#not a leaf
    if (length(toss)==0) return(tree)   #all the tree is retained

    newx <- snip.rpart(tree, toss)

    ## Now cut down the CP table
    temp <- pmax(tree$cptable[,1], cp)
    keep <- match(unique(temp), temp)
    newx$cptable <- tree$cptable[keep,,drop=FALSE]
    newx$cptable[max(keep),1] <- cp

    newx
}
# SCCS @(#)prune.s	1.2 02/12/98
# This should be part of Splus proper -- make prune a method
prune <- function(tree, ...)  UseMethod("prune")
#SCCS  %W% %G%
residuals.rpart <- function(object, type = c("usual", "pearson", "deviance"), ...)
    {
    if(!inherits(object, "rpart"))
	    stop("Not legitimate rpart object")

    y <- object$y
    if (is.null(y)) y <- model.extract(model.frame(object), "response")
    frame <- object$frame
    type <- match.arg(type)
    if(is.na(match(type, c("usual", "pearson", "deviance"))))
                stop("Don't know about this type of residual")

    if (object$method=='class') {
	ylevels <- attr(object, "ylevels")
	nclass <- length(ylevels)

        if(type == "usual") {
                yhat <- frame$yval[object$where]
		loss <- object$parms$loss
		}
        else {
	    yprob <- frame$yval2[object$where, 1 + nclass + 1:nclass]
	    yhat <- yprob[cbind(seq(y), unclass(y))]
	    }
        resid  <- switch(type,
                usual = loss[cbind(y, yhat)],
                pearson = (1 - yhat)/yhat,
                deviance = -2 * log(yhat))
       }

    else if (object$method=='poisson' || object$method=='exp') {
	lambda <- (object$frame$yval)[object$where]
	time   <- y[,1]  # observation time in new data
	events <- y[,2]  # number of events, in new data
	expect <- lambda * time #expected number of events
	temp <- ifelse(expect==0, .0001, 0)  #failsafe for log(0)

	resid <- switch(type,
			usual = events - expect,
			pearson = (events - expect)/sqrt(temp),
			deviance= sign(events- expect) *
			   sqrt(2*(events*log(events/temp) - (events-expect)))
			)
	}

    else  resid <- y - frame$yval[object$where]


    names(resid) <- names(y)
    #Expand out the missing values in the result
    if (!is.null(object$na.action))
	resid <- naresid(object$na.action, resid)

    resid
    }
#SCCS @(#)rpart.anova.s	1.4 05/02/01
rpart.anova <- function(y, offset, parms, wt) {
    if (!is.null(offset)) y <- y-offset
    list(y=y, parms=0, numresp=1, numy=1,
	 summary= function(yval, dev, wt, ylevel, digits ) {
	     paste("  mean=", formatg(yval, digits),
		   ", MSE=" , formatg(dev/wt, digits),
		   sep='')
	     },
	 text= function(yval, dev, wt, ylevel, digits, n, use.n ) {
	     if(use.n) {paste(formatg(yval,digits),"\nn=", n,sep="")} else
	               {paste(formatg(yval,digits))}}

	 )
    }
#SCCS @(#)rpart.branch.s	1.2 01/25/97
#
# Compute the "branches" to be drawn for an rpart object
#
rpart.branch <- function(x, y, node, branch) {
    if (missing(branch)) {
	if (exists(parms <-paste(".rpart.parms", dev.cur(), sep="." ),
                   envir=.GlobalEnv)) {
#	    parms <- get(parms, frame=0)
            parms <- get(parms, envir=.GlobalEnv)
            branch <- parms$branch
	    }
	else branch <- 0
        }

    # Draw a series of horseshoes, left son, up, over, down to right son
    #   NA's in the vector cause lines() to "lift the pen"
    is.left <- (node%%2 ==0)        #left hand sons
    node.left <- node[is.left]
    parent <- match(node.left/2, node)
    sibling <- match(node.left+1, node)
    temp <- (x[sibling] - x[is.left])*(1-branch)/2
    xx <- rbind(x[is.left], x[is.left]+ temp,
                x[sibling]- temp, x[sibling], NA)
    yy <- rbind(y[is.left], y[parent], y[parent], y[sibling], NA)
    list(x=xx, y=yy)
    }
#SCCS @(#)rpart.class.s	1.7 07/05/01

rpart.class <- function (y, offset, parms, wt) 
{
    if (!is.null(offset)) 
        stop("No offset allowed in classification models")
    fy <- as.factor(y)
    y <- as.integer(fy)
    numclass <- max(y[!is.na(y)])
    counts <- tapply(wt, factor(y, levels = 1:numclass), sum)
    counts <- ifelse(is.na(counts), 0, counts)
    numresp <- 1 + numclass
    if (missing(parms) || is.null(parms)) 
        parms <- list(prior = counts/sum(counts), loss = matrix(rep(1, 
            numclass^2) - diag(numclass), numclass), split = 1)
    else if (is.list(parms)) {
        if (is.null(names(parms))) 
            stop("The parms list must have names")
        temp <- pmatch(names(parms), c("prior", "loss", "split"), 
            nomatch = 0)
        if (any(temp == 0)) 
            stop(paste("parms component not matched:", (names(parms))[temp == 
                0]))
        names(parms) <- c("prior", "loss", "split")[temp]
        if (is.null(parms$prior)) 
            temp <- c(counts/sum(counts))
        else {
            temp <- parms$prior
            if (sum(temp) != 1) 
                stop("Priors must sum to 1")
            if (any(temp < 0)) 
                stop("Priors must be >= 0")
            if (length(temp) != numclass) 
                stop("Wrong length for priors")
        }
        if (is.null(parms$loss)) 
            temp2 <- 1 - diag(numclass)
        else {
            temp2 <- parms$loss
            if (length(temp2) != numclass^2) 
                stop("Wrong length for loss matrix")
            temp2 <- matrix(temp2, ncol = numclass)
            if (any(diag(temp2) != 0)) 
                stop("Loss matrix must have zero on diagonals")
            if (any(temp2 < 0)) 
                stop("Loss matrix cannot have negative elements")
            if (any(rowSums(temp2) == 0)) 
                stop("Loss matrix has a row of zeros")
        }
        if (is.null(parms$split)) 
            temp3 <- 1
        else {
            temp3 <- pmatch(parms$split, c("gini", "information"))
            if (is.null(temp3)) 
                stop("Invalid splitting rule")
        }
        parms <- list(prior = temp, loss = matrix(temp2, numclass), 
            split = temp3)
    }
    else stop("Parameter argument must be a list")
    list(y = y, parms = parms, numresp = numclass + 1, counts = counts, 
        ylevels = levels(fy), numy = 1, print = function(yval, 
            ylevel, digits) {
            if (is.null(ylevel)) temp <- as.character(yval[, 
                1]) else temp <- ylevel[yval[, 1]]
            nclass <- (ncol(yval) - 1)/2
            if (nclass < 5) {
                yprob <- format(yval[, 1 + nclass + 1:nclass], 
                  digits = digits, nsmall = digits)
            } else yprob <- formatg(yval[, 1 + nclass + 1:nclass], 
                digits = 2)
            if (is.null(dim(yprob))) yprob <- matrix(yprob, ncol = length(yprob))
            temp <- paste(temp, " (", yprob[, 1], sep = "")
            for (i in 2:ncol(yprob)) temp <- paste(temp, yprob[, 
                i], sep = " ")
            temp <- paste(temp, ")", sep = "")
            temp
        }, summary = function(yval, dev, wt, ylevel, digits) {
            nclass <- (ncol(yval) - 1)/2
            group <- yval[, 1]
            counts <- yval[, 1 + (1:nclass)]
            yprob <- yval[, 1 + nclass + 1:nclass]
            if (!is.null(ylevel)) group <- ylevel[group]
            temp1 <- formatg(counts, format = "%5g")
            temp2 <- formatg(yprob, format = "%5.3f")
            if (nclass > 1) {
                temp1 <- apply(matrix(temp1, ncol = nclass), 
                  1, paste, collapse = " ")
                temp2 <- apply(matrix(temp2, ncol = nclass), 
                  1, paste, collapse = " ")
            }
            paste("  predicted class=", format(group, justify = "left"), 
                "  expected loss=", formatg(dev/wt, digits), 
                "\n", "    class counts: ", temp1, "\n", "   probabilities: ", 
                temp2, sep = "")
        }, text = function(yval, dev, wt, ylevel, digits, n, 
            use.n) {
            nclass <- (ncol(yval) - 1)/2
            group <- yval[, 1]
            counts <- yval[, 1 + (1:nclass)]
            if (!is.null(ylevel)) group <- ylevel[group]
            temp1 <- formatg(counts, digits)
            if (nclass > 1) {
                temp1 <- apply(matrix(temp1, ncol = nclass), 
                  1, paste, collapse = "/")
            }
            if (use.n) {
                out <- paste(format(group, justify = "left"), 
                  "\n", temp1, sep = "")
            } else {
                out <- format(group, justify = "left")
            }
            return(out)
        }, bar = function(yval2){
            yval2[,2:((ncol(yval2) + 1)/2)]
        }
        )
}


#SCCS @(#)rpart.control.s	1.10 07/05/01

rpart.control <- function (minsplit = 5, minbucket = round(minsplit/3), cp = 0.01, 
    maxcompete = 4, maxsurrogate = 5, usesurrogate = 2, xval = 10, 
    surrogatestyle = 0, maxdepth = 30, ...) 
{
    if (maxcompete < 0) {
        warning("The value of maxcompete supplied is <0; the value 0 was used instead")
        maxcompete <- 0
    }
    if (any(xval < 0)) {
        warning("The value of xval supplied is <0; the value 0 was used instead")
        xval <- 0
    }
    if (maxdepth > 30) 
        stop("Maximum depth is 30")
    if (maxdepth < 1) 
        stop("Maximum depth must be at least 1")
    if (missing(minsplit) && !missing(minbucket)) 
        minsplit <- minbucket * 3
    if ((usesurrogate < 0) || (usesurrogate > 2)) {
        warning("The value of usesurrogate supplied was out of range,", 
            "the default value of 2 is used instead.")
        usesurrogate <- 2
    }
    if ((surrogatestyle < 0) || (surrogatestyle > 1)) {
        warning("The value of surrogatestyle supplied was out of range,", 
            "the default value of 0 is used instead.")
        surrogatestyle <- 0
    }
    list(minsplit = minsplit, minbucket = minbucket, cp = cp, 
        maxcompete = maxcompete, maxsurrogate = maxsurrogate, 
        usesurrogate = usesurrogate, surrogatestyle = surrogatestyle, 
        maxdepth = maxdepth, xval = xval)
}


#SCCS @(#)rpart.exp.s	1.8 07/05/01
# rescaled exponential splitting
#  The survival object 'y' is rescaled so that
#    a. overall death rate is 1.0
#    b. death rate within small intervals of time is 1
#  The first makes printouts easier to read, since the rates in subnodes are
#  numbers like "1.23" (23% higher event rate than the root node) instead
#  of "0.0014" (which requires looking back at the root node rate and
#  dividing).  The second makes the data appear to be Poisson, and causes
#  the early splits at least to be equivalent to the local full likelihood
#  method of LeBlanc and Crowley
#
rpart.exp <- function(y, offset, parms, wt) {

    if (!inherits(y, "Surv"))
	   stop("Response must be a survival object - use the Surv() function")

    ny <- ncol(y)
    n  <- nrow(y)

    status <- y[,ny]
    if (any(y[,1]<=0)) stop("Observation time must be >0")
    if (all(status==0)) stop("No deaths in data set")
    time <- y[ ,ny-1]

    # Make my table of time intervals.  The first goes from 0 to the first
    #   death, the next from death 2 to death 3, ..., and the last from
    #   "next to last death" to "max time in dataset".
    # We also need to avoid a pathological condition in some data sets, where
    #   two death times differ by a trivial amount, e.g., 10^-16, perhaps due
    #   to roundoff error in creating the input data.  Ammalgamate such
    #   intervals.  This turns out to be hard to do in S, but easy in C
    dtimes <- sort(unique(time[status==1]))        # unique death times
    temp <- .C('rpartexp2',
	       as.integer(length(dtimes)),
	       as.double(dtimes),
	       as.double(.Machine$double.eps),
	       keep=integer(length(dtimes)), PACKAGE="mvpart")$keep
    dtimes <- dtimes[temp==1]

    # For the sake of speed, restrict the number of intervals to be <1000.
    #   (Actually, anything >100 is probably overkill for the
    #   actual task at hand, which is to approximately scale to exponential).
    if (length(dtimes) > 1000) dtimes <- quantile(dtimes, 0:1000/1000)

    # The last interval goes to the max time in the data set
    itable <- c(0, dtimes[-length(dtimes)], max(time)) # set of intervals

    drate1 <- function(n, ny, y, wt, itable) {
	# Compute the death rate within each of the intervals
	#  The pyears2 routine is part of the survival library
	ngrp <- length(itable) -1
	temp <- .C('pyears2',
		   as.integer(n),
		   as.integer(ny),
		   as.integer(1),
		   as.double (y),
		   as.double(wt),
		   as.integer(1),
		   as.integer(0),
		   as.integer(ngrp),
		   as.double(itable),
		   as.double(rep(0., n)),
		   pyears = double(ngrp),
		   pn     = double(ngrp),
		   pcount = double(ngrp),
		   offtable= double(1), PACKAGE="survival")[11:14]
	rates <- temp$pcount / temp$pyears
	rates
	}

    drate2 <- function(n, ny, y, wt, itable) {
	# An alternative to the drate1 function
	# Why?  The pyears2 routine changed in 6/2001, with the inclusion
	#  of case weights.  We need the newer version.  If you have the
	#  older version of the survival library, the above will crash S.
	# How to tell -- list the pyears function, and see whether it's
	#  call to pyears2 has weights in the argument list.
	#
	time <- y[, ny-1]
	status <- y[,ny]
	ilength <- diff(itable)                   #lengths of intervals
	ngrp <- length(ilength)                   #number of intervals

	# The code below is as opaque as any I've written, all in the
	#  service of "no for loops".
	# First, 'index' gives the time interval (as defined by itable)
	#  in which the end of each observation's follow-up (time) lies.
	#  Then 'itime' will be the amount of time spent in that last
	#  interval, which is of course somewhat less than ilength.
	index <- unclass(cut(time, itable, include.lowest=TRUE))
	itime <- time - itable[index]
	if (ny ==3) {
	    # there is both a start time and a stop time
	    #  compute the amount of time NOT spent in the interval that
	    #  the start time lies in.
	    stime <- y[,1]   #start time for each interval
	    index2<- unclass(cut(stime, itable, include.lowest=TRUE))
	    itime2<- stime - itable[index2]
	    }

	# Compute the amount of person-years in each of the intervals
	#   This is:  (width of interval) * (number of "time" elements that
	#                                     end in an interval farther right)
	#            + (ending times in this interval)
	# By construction, I know that there is at least 1 obs in each of the
	#  intervals, so tab1 is of a determined length
	tab1 <- table(index)
	temp <- rev(cumsum(rev(tab1)))  #cumsum, counting from the right
	pyears <- ilength * c(temp[-1], 0) +	 tapply(itime, index, sum)
	if (ny==3) {
	    #subtract off the time before "start"
	    tab2 <- table(index2, levels=1:ngrp) #force the length of tab2
	    temp <- rev(cumsum(rev(tab2)))
	    py2  <-  ilength * c(0, temp[-ngrp]) +  tapply(itime2, index2, sum)
	    pyears <- pyears - py2
	    }

	deaths <- tapply(status, index, sum)
	rate <- deaths/pyears   #hazard rate in each interval
	rate
	}

    #
    # Now, compute the "new y" for each observation.
    #  This is a stretching of the time axis
    # The cumulative hazard over each interval is rate*length(interval),
    #  and is the basis of the rescaling.
    rate <- drate2(n, ny, y, wt, itable)
    cumhaz <- cumsum(c(0, rate*diff(itable)))
    newy <- approx(itable, cumhaz, time)$y
    if (ny==3) {
	newy <- newy - approx(itable, cumhaz, stime)$y
	}

    if (length(offset)==n)  newy <- newy * exp(offset)

    if (missing(parms)) parms <- c(shrink=1, method=1)
    else {
	parms <- as.list(parms)
        if(is.null(names(parms))) stop("You must input a named list for parms")
        parmsNames <- c("method", "shrink")
        indx <- pmatch(names(parms), parmsNames, nomatch= 0)
        if (any(indx==0))
            stop(paste("parms component not matched: ",
		       names(parms)[indx==0]))
	else names(parms) <- parmsNames[indx]

	if (is.null(parms$method)) method <- 1
	else method <- pmatch(parms$method, c("deviance", "sqrt"))
	if (is.na(method)) stop("Invalid error method for Poisson")

	if (is.null(parms$shrink)) shrink <- 2-method
	else shrink <- parms$shrink
	if (!is.numeric(shrink) || shrink < 0)
		stop("Invalid shrinkage value")
	parms <- c(shrink=shrink, method=method)
	}
    list(y=cbind(newy, y[,2]), parms=parms, numresp=2, numy=2,
	 summary= function(yval, dev, wt, ylevel, digits) {
	     paste("  events=", formatg(yval[,2]),
		",  estimated rate=" , formatg(yval[,1], digits),
		" , mean deviance=",formatg(dev/wt, digits),
		sep = "")
	     },
	 text= function(yval, dev, wt, ylevel, digits, n, use.n) {
	     if(use.n) {paste(formatg(yval[,1],digits),"\n",
				formatg(yval[,2]),"/",n,sep="")} else
		    {paste(formatg(yval[,1],digits))}
	     })
    }
#SCCS  @(#)rpart.matrix.s	1.6 04/02/01
#
# This differs from tree.matrix in xlevels -- we don't keep NULLS in
#   the list for all of the non-categoricals
#
rpart.matrix <- function(frame)
    {
    if(!inherits(frame, "data.frame"))
	    return(as.matrix(frame))
    frame$"(weights)" <- NULL
    terms <- attr(frame, "terms")
    if(is.null(terms)) predictors <- names(frame)
    else {
	a <- attributes(terms)
	predictors <- as.character(a$variables)[-1] # R change
	removals <- NULL
	if((TT <- a$response) > 0) {
	    removals <- TT
	    frame[[predictors[TT]]] <- NULL
	    }
	if(!is.null(TT <- a$offset)) {
	    removals <- c(removals, TT)
	    frame[[predictors[TT]]] <- NULL
	    }
	if(!is.null(removals)) predictors <- predictors[ - removals]
        labels <- a$term.labels
	if(abs(length(labels)-length(predictors))>0)
	  predictors <- predictors[match(labels,predictors)]
	}

    factors <- sapply(frame, function(x) !is.null(levels(x)))
    characters <- sapply(frame, is.character)
    if(any(factors | characters)) {
	# change characters to factors
	for (preds in predictors[characters])
		frame[[preds]] <- as.factor(frame[[preds]])
        factors <- factors | characters
        column.levels <- lapply(frame[factors], levels)

	# Now make them numeric
	for (preds in predictors[factors])
	     frame[[preds]] <- as.numeric(frame[[preds]])
	x <- as.matrix(frame)
	attr(x, "column.levels") <- column.levels
	}
    else x <- as.matrix(frame[predictors])
    class(x) <- "rpart.matrix"
    x
    }


#SCCS @(#)rpart.mrt.s	1.4 05/07/03

rpart.mrt <- function (y, offset, parms, wt) 
{
    if (!is.null(offset)) 
        y <- y - offset
	    ny <- ifelse(is.null(dim(y)), 1, dim(y)[2])
    	list(y = y, parms = 0, numresp = ny, numy = ny, summary = function(yval, 
        dev, wt, ylevel, digits) {
        paste("  Means=", apply(formatg(yval, digits - 3), 1, 
            paste, collapse = ",", sep = ""), ", Summed MSE=", 
            formatg(dev/wt, digits), sep = "")
    }, text = function(yval, dev, wt, ylevel, digits, n, use.n) {
        if (use.n) {
            paste(formatg(dev, digits), " : n=", 
                n, sep = "")
        } 
          else {
            paste(formatg(dev, digits))
        }
    }, bar = function(yval2) {
            yval2
    })
}


#SCCS @(#)rpart.dist.s	1.4 05/07/03

rpart.dist <- function (y, offset, parms, wt) 
{
    if (!is.null(offset)) 
        y <- y - offset
    	ny <-  1
    	list(y = y, parms = 0, numresp = ny, numy = ny, summary = function(yval, 
        dev, wt, ylevel, digits) {
        paste("  Means=", apply(formatg(yval, digits - 3), 1, 
            paste, collapse = ",", sep = ""), ", Summed MSE=", 
            formatg(dev/wt, digits), sep = "")
    }, text = function(yval, dev, wt, ylevel, digits, n, use.n) {
        if (use.n) {
            paste(formatg(dev, digits), " : n=", 
                n, sep = "")
        } 
          else {
            paste(formatg(dev, digits))
        }
    }, bar = function(yval2) {
            yval2
    })
}



#SCCS @(#)rpart.poisson.s	1.6 07/05/01
rpart.poisson <- function(y, offset, parms, wt) {
    if (is.matrix(y)) {
	if (ncol(y)!=2) stop("response must be a 2 column matrix or a vector")
	if (!is.null(offset)) y[,1] <- y[,1] * exp(offset)
	}
    else {
	if (is.null(offset)) y <- cbind(1,y)
	else  y <- cbind( exp(offset), y)
	}
    if (any(y[,1] <=0)) stop("Observation time must be >0")
    if (any(y[,2] <0))  stop("Number of events must be >=0")

    if (missing(parms)) parms <- c(shrink=1, method=1)
    else {
	parms <- as.list(parms)
	if(is.null(names(parms))) stop("You must input a named list for parms")
	parmsNames <- c("method", "shrink")
	indx <- pmatch(names(parms), parmsNames, nomatch= 0)
	if (any(indx==0))
               stop(paste("parms component not matched: ",
			  names(parms)[indx==0]))
	else names(parms) <- parmsNames[indx]

	if (is.null(parms$method)) method <- 1
	else method <- pmatch(parms$method, c("deviance", "sqrt"))
	if (is.null(method)) stop("Invalid error method for Poisson")

	if (is.null(parms$shrink)) shrink <- 2- method
	else shrink <- parms$shrink

	if (!is.numeric(shrink) || shrink <0)
		stop("Invalid shrinkage value")
	parms <- c(shrink=shrink, method=method)
	}

    list(y=y, parms=parms, numresp=2, numy=2,
	 summary= function(yval, dev, wt, ylevel, digits) {
	     paste("  events=", formatg(yval[,2]),
		",  estimated rate=" , formatg(yval[,1], digits),
		" , mean deviance=",formatg(dev/wt, digits),
		sep = "")
	     },
	 text= function(yval, dev, wt, ylevel, digits, n, use.n) {
	     if(use.n) {paste(formatg(yval[,1],digits),"\n",
				formatg(yval[,2]),"/",n,sep="")} else
		    {paste(formatg(yval[,1],digits))}}
	 )

    }
# SCCS  @(#) rpart.s	1.35 07/05/01
#
#  The recursive partitioning function, for S
#

rpart <- function (formula, data = NULL, weights, subset, na.action = na.rpart, 
    method, model = FALSE, x = FALSE, y = TRUE, parms, control, 
    cost, ...) 
{
    call <- match.call()
    if (is.data.frame(model)) {
        m <- model
        model <- FALSE
    }
    else {
        m <- match.call(expand = FALSE)
        m$model <- m$method <- m$control <- NULL
        m$x <- m$y <- m$parms <- m$... <- NULL
        m$cost <- NULL
        m$na.action <- na.action
        m[[1]] <- as.name("model.frame.default")
        m <- eval(m, parent.frame())
    }
    Terms <- attr(m, "terms")
    if (any(attr(Terms, "order") > 1)) 
        stop("Trees cannot handle interaction terms")
    Y <- model.extract(m, "response")
    wt <- model.extract(m, "weights")
    if (length(wt) == 0) 
        wt <- rep(1, nrow(m))
    offset <- attr(Terms, "offset")
    X <- rpart.matrix(m)
    nobs <- nrow(X)
    nvar <- ncol(X)
    if (missing(method)) {
        if (is.factor(Y) || is.character(Y)) 
            method <- "class"
        else if (inherits(Y, "Surv")) 
            method <- "exp"
        else if (is.matrix(Y)) 
            method <- "mrt"
        else method <- "anova"
    }
    if (is.list(method)) {
        mlist <- method
        method <- "user"
        if (missing(parms)) 
            init <- mlist$init(Y, offset, wt = wt)
        else init <- mlist$init(Y, offset, parms, wt)
        method.int <- 6
        keep <- rpartcallback(mlist, nobs, init)
    }
    else {
        method.int <- pmatch(method, c("anova", "mrt", "poisson", 
            "class", "dist", "exp"))
        if (is.na(method.int)) 
            stop("Invalid method")
        method <- c("anova", "mrt", "poisson", "class", "dist", 
            "exp")[method.int]
        if (method.int == 6) 
            method.int <- 3
        if (missing(parms)) 
            init <- (get(paste("rpart", method, sep = ".")))(Y, 
                offset, , wt)
        else init <- (get(paste("rpart", method, sep = ".")))(Y, 
            offset, parms, wt)
    }
    Y <- init$y
    if (method == "dist") {
        Y <- Y[row(Y) > col(Y)]
        init$y <- init$y[row(init$y) > col(init$y)]
    }
    xlevels <- attr(X, "column.levels")
    cats <- rep(0, ncol(X))
    if (!is.null(xlevels)) {
        cats[match(names(xlevels), dimnames(X)[[2]])] <- unlist(lapply(xlevels, 
            length))
    }
    extraArgs <- list(...)
    if (length(extraArgs)) {
        controlargs <- names(formals(rpart.control))
        indx <- match(names(extraArgs), controlargs, nomatch = 0)
        if (any(indx == 0)) 
            stop(paste("Argument", names(extraArgs)[indx == 0], 
                "not matched"))
    }
    controls <- rpart.control(...)
    if (!missing(control)) 
        controls[names(control)] <- control
    xval <- controls$xval
    if (is.null(xval) || (length(xval) == 1 && xval == 0) || 
        method == "user" || method == "dist") {
        xgroups <- 0
#   Set xval to 0 for dist splitting and reset controls$xval -- GD 12/03
        xval <- 0
		controls$xval <- xval 
           }
    else if (length(xval) == 1) {
        xgroups <- sample(rep(1:xval, length = nobs), nobs, replace = FALSE)
    }
    else if (length(xval) == nobs) {
        xgroups <- xval
        xval <- length(unique(xgroups))
    }
    else {
        if (!is.null(attr(m, "na.action"))) {
            temp <- as.integer(attr(m, "na.action"))
            xval <- xval[-temp]
            if (length(xval) == nobs) {
                xgroups <- xval
                xval <- length(unique(xgroups))
            }
            else stop("Wrong length for xval")
        }
        else stop("Wrong length for xval")
    }
    if (missing(cost)) 
        cost <- rep(1, nvar)
    else {
        if (length(cost) != nvar) 
            stop("Cost vector is the wrong length")
        if (any(cost <= 0)) 
            stop("Cost vector must be positive")
    }
    tfun <- function(x) {
        if (is.matrix(x)) 
            rep(is.ordered(x), ncol(x))
        else is.ordered(x)
    }
    isord <- unlist(lapply(m[attr(Terms, "term.labels")], tfun))
    rpfit <- .C("s_to_rp", n = as.integer(nobs), nvarx = as.integer(nvar), 
        ncat = as.integer(cats * (!isord)), method = as.integer(method.int), 
        as.double(unlist(controls)), parms = as.double(unlist(init$parms)), 
        as.integer(xval), as.integer(xgroups), as.double(t(init$y)), 
        as.double(X), as.integer(!is.finite(X)), error = character(1), 
        wt = as.double(wt), as.integer(init$numy), as.double(cost), 
        NAOK = TRUE, PACKAGE = "mvpart")
    if (rpfit$n == -1) 
        stop(rpfit$error)
    nodes <- rpfit$n
    nsplit <- rpfit$nvarx
    numcp <- rpfit$method
    ncat <- rpfit$ncat[1]
    numresp <- init$numresp
    if (nsplit == 0) 
        xval <- 0
    cpcol <- if (xval > 0 && nsplit > 0) 
        5
    else 3
    if (ncat == 0) 
        catmat <- 0
    else catmat <- matrix(integer(1), ncat, max(cats))
    rp <- .C("s_to_rp2", as.integer(nobs), as.integer(nsplit), 
        as.integer(nodes), as.integer(ncat), as.integer(cats * 
            (!isord)), as.integer(max(cats)), as.integer(xval), 
        which = integer(nobs), cptable = matrix(double(numcp * 
            cpcol), nrow = cpcol), dsplit = matrix(double(1), 
            nsplit, 3), isplit = matrix(integer(1), nsplit, 3), 
        csplit = catmat, dnode = matrix(double(1), nodes, 3 + 
            numresp), inode = matrix(integer(1), nodes, 6), PACKAGE = "mvpart")
    tname <- c("<leaf>", dimnames(X)[[2]])
    if (cpcol == 3) 
        temp <- c("CP", "nsplit", "rel error")
    else temp <- c("CP", "nsplit", "rel error", "xerror", "xstd")
    dimnames(rp$cptable) <- list(temp, 1:numcp)
    dn1 <- if (nsplit == 0) 
        character(0)
    else tname[rp$isplit[, 1] + 1]
    splits <- matrix(c(rp$isplit[, 2:3], rp$dsplit), ncol = 5, 
        dimnames = list(dn1, c("count", "ncat", "improve", "index", 
            "adj")))
    index <- rp$inode[, 2]
    nadd <- sum(isord[rp$isplit[, 1]])
    if (nadd > 0) {
        newc <- matrix(integer(1), nadd, max(cats))
        cvar <- rp$isplit[, 1]
        indx <- isord[cvar]
        cdir <- splits[indx, 2]
        ccut <- floor(splits[indx, 4])
        splits[indx, 2] <- cats[cvar[indx]]
        splits[indx, 4] <- ncat + 1:nadd
        for (i in 1:nadd) {
            newc[i, 1:(cats[(cvar[indx])[i]])] <- -1 * as.integer(cdir[i])
            newc[i, 1:ccut[i]] <- as.integer(cdir[i])
        }
        if (ncat == 0) 
            catmat <- newc
        else catmat <- rbind(rp$csplit, newc)
        ncat <- ncat + nadd
    }
    else catmat <- rp$csplit
    if (nsplit == 0) {
        frame <- data.frame(row.names = 1, var = "<leaf>", n = rp$inode[, 
            5], wt = rp$dnode[, 3], dev = rp$dnode[, 1], yval = rp$dnode[, 
            4], complexity = rp$dnode[, 2], ncompete = pmax(0, 
            rp$inode[, 3] - 1), nsurrogate = rp$inode[, 4])
    }
    else {
        temp <- ifelse(index == 0, 1, index)
        svar <- ifelse(index == 0, 0, rp$isplit[temp, 1])
        frame <- data.frame(row.names = rp$inode[, 1], var = factor(svar, 
            0:ncol(X), tname), n = rp$inode[, 5], wt = rp$dnode[, 
            3], dev = rp$dnode[, 1], yval = rp$dnode[, 4], complexity = rp$dnode[, 
            2], ncompete = pmax(0, rp$inode[, 3] - 1), nsurrogate = rp$inode[, 
            4])
        if (method == "mrt") 
            frame$yval <- apply(rp$dnode[, -c(1:3)], 1, mean)
    }
    if (method.int == 4) {
        numclass <- init$numresp - 1
        temp <- rp$dnode[, -(1:4)] %*% diag(init$parms$prior * 
            sum(init$counts)/pmax(1, init$counts))
        yprob <- temp/rowSums(temp)
        yval2 <- matrix(rp$dnode[, -(1:3)], ncol = numclass + 
            1)
        frame$yval2 <- cbind(yval2, yprob)
    }
    else if (method.int == 2) 
        frame$yval2 <- rp$dnode[, -(1:3)]
    if (is.null(init$summary)) 
        stop("Initialization routine is missing the summary function")
    if (is.null(init$print)) 
        functions <- list(summary = init$summary)
    else functions <- list(summary = init$summary, print = init$print)
    if (!is.null(init$text)) 
        functions <- c(functions, list(text = init$text))
    if (!is.null(init$bar)) 
        functions <- c(functions, list(bar = init$bar))
    if (method == "user") 
        functions <- c(functions, mlist)
    where <- rp$which
    names(where) <- row.names(m)
    if (nsplit == 0) {
        ans <- list(frame = frame, where = where, call = call, 
            terms = Terms, cptable = t(rp$cptable), method = method, 
            parms = init$parms, control = controls, functions = functions)
    }
    else {
        ans <- list(frame = frame, where = where, call = call, 
            terms = Terms, cptable = t(rp$cptable), splits = splits, 
            method = method, parms = init$parms, control = controls, 
            functions = functions)
    }
    if (ncat > 0) 
        ans$csplit <- catmat + 2
    if (model) {
        ans$model <- m
        if (missing(y)) 
            y <- FALSE
    }
    if (y) 
        ans$y <- Y
    if (x) {
        ans$x <- X
        ans$wt <- wt
    }
    ans$ordered <- isord
    if (!is.null(attr(m, "na.action"))) 
        ans$na.action <- attr(m, "na.action")
    if (!is.null(xlevels)) 
        attr(ans, "xlevels") <- xlevels
    if (method == "class") 
        attr(ans, "ylevels") <- init$ylevels
    class(ans) <- "rpart"
    ans
}


#  SCCS %W% %G%
#This routine sets up the callback code for user-written split
#  routines in rpart
#
rpartcallback <- function(mlist, nobs, init)
{
    if (length(mlist) < 3)
        stop("User written methods must have 3 functions")
    if (is.null(mlist$init) || typeof(mlist$init) != 'closure')
        stop("User written method does not contain an init function")
    if (is.null(mlist$split) || typeof(mlist$split) != 'closure')
        stop("User written method does not contain a split function")
    if (is.null(mlist$eval) || typeof(mlist$eval) != 'closure')
        stop("User written method does not contain an eval function")

    user.eval <- mlist$eval
    user.split <- mlist$split

    numresp <- init$numresp
    numy <-  init$numy
    parms <- init$parms

    #
    # expr2 is an expression that will call the user "evaluation"
    #   function, and check that what comes back is valid
    # expr1 does the same for the user "split" function
    #
    # For speed in the C interface, yback, xback, and wback are
    #  fixed S vectors of a fixed size, and nback tells us how
    #  much of the vector is actually being used on this particular
    #  callback.
    #
    if (numy==1) {
        expr2 <- quote({
            temp <- user.eval(yback[1:nback], wback[1:nback], parms)
            if (length(temp$label) != numresp)
                stop("User eval function returned invalid label")
            if (length(temp$deviance) !=1)
                stop("User eval function returned invalid deviance")
            as.numeric(as.vector(c(temp$deviance, temp$label)))
        })
        expr1 <- quote({
            if (nback <0) { #categorical variable
                n2 <- -1*nback
                temp  <- user.split(yback[1:n2], wback[1:n2],
                                    xback[1:n2], parms, FALSE)
                ncat <- length(unique(xback[1:n2]))
                if (length(temp$goodness) != ncat-1 ||
                    length(temp$direction) != ncat)
                    stop("Invalid return from categorical split fcn")
            }

            else {
                temp <- user.split(yback[1:nback], wback[1:nback],
                                   xback[1:nback], parms, TRUE)
                if (length(temp$goodness) != (nback-1))
                    stop("User split function returned invalid goodness")
                if (length(temp$direction) != (nback-1))
                    stop("User split function returned invalid direction")
            }
            as.numeric(as.vector(c(temp$goodness, temp$direction)))
        })
    }
    else {
        expr2 <- quote({
            tempy <- matrix(yback[1:(nback*numy)], ncol=numy)
            temp <- user.eval(tempy, wback[1:nback], parms)
            if (length(temp$label) != numresp)
                stop("User eval function returned invalid label")
            if (length(temp$deviance) != 1)
                stop("User eval function returned invalid deviance")
            as.numeric(as.vector(c(temp$deviance, temp$label)))
        })
        expr1 <- quote({
            if (nback <0) { #categorical variable
                n2 <- -1*nback
                tempy <- matrix(yback[1:(n2*numy)], ncol=numy)
                temp  <- user.split(tempy, wback[1:n2], xback[1:n2],
                                    parms, FALSE)
                ncat <- length(unique(xback[1:n2]))
                if (length(temp$goodness) != ncat-1 ||
                    length(temp$direction) != ncat)
                    stop("Invalid return from categorical split fcn")
            }
            else {
                tempy <- matrix(yback[1:(nback*numy)], ncol=numy)
                temp <- user.split(tempy, wback[1:nback],xback[1:nback],
                                   parms, TRUE)
                if (length(temp$goodness) != (nback-1))
                    stop("User split function returned invalid goodness")
                if (length(temp$direction) != (nback-1))
                    stop("User split function returned invalid direction")
            }
            as.numeric(as.vector(c(temp$goodness, temp$direction)))
        })
    }
    #
    #  The vectors nback, wback, xback and yback will have their
    #  contents constantly re-inserted by C code.  It's one way to make
    #  things very fast.  It is dangerous to do this, so they
    #  are tossed into a separate frame to isolate them.  Evaluations of
    #  the above expressions occur in that frame.
    #
    rho <- new.env()
    assign("nback", integer(1), envir = rho)
    assign("wback", double(nobs), envir = rho)
    assign("xback", double(nobs), envir = rho)
    assign("yback", double(nobs*numy), envir = rho)
    assign("user.eval", user.eval, envir = rho)
    assign("user.split", user.split, envir = rho)
    assign("numy", numy, envir = rho)
    assign("numresp", numresp, envir = rho)
    assign("parms", parms, envir = rho)
    .Call("init_rpcallback", rho, as.integer(numy), as.integer(numresp),
          expr1, expr2, PACKAGE = "mvpart")
    list(expr1 = expr1, expr2 = expr2, rho = rho)
}
#SCCS @(#)rpartco.s	1.7 02/07/00
# Compute the x-y coordinates for a tree

rpartco <- function(tree, parms =  paste(".rpart.parms", dev.cur(), sep = "."))
    {

    frame <- tree$frame
    method <- tree$method
    node <- as.numeric(row.names(frame))
    depth <- tree.depth(node)
    is.leaf <- (frame$var == '<leaf>')
    if (exists(parms, envir=.GlobalEnv)) {
	parms <- get(parms, envir=.GlobalEnv)
	uniform <- parms$uniform
	nspace <-parms$nspace
	minbranch <- parms$minbranch
	}
    else {
	uniform <- FALSE
	nspace <- -1
	minbranch <- .3
        }

    if(uniform) y <- (1 + max(depth) -depth) / max(depth,4)
    else {                    #make y- (parent y) = change in deviance
	y <- dev <- frame$dev
        temp <- split(seq(node), depth)     #depth 0 nodes, then 1, then ...
        parent <- match(floor(node/2), node)
        sibling <- match(ifelse(node %% 2, node - 1, node + 1), node)

	# assign the depths
        for(i in temp[-1]) {
	    temp2 <- dev[parent[i]] - (dev[i] + dev[sibling[i]])
            y[i] <- y[parent[i]] - temp2
	    }
	#
	# For some problems, classification & loss matrices in particular
	#   the gain from a split may be 0.  This is ugly on the plot.
	# Hence the "fudge" factor of  .3* the average step
	#
	fudge <-  minbranch * diff(range(y)) / max(depth)
        for(i in temp[-1]) {
	    temp2 <- dev[parent[i]] - (dev[i] + dev[sibling[i]])
	    haskids <- !(is.leaf[i] & is.leaf[sibling[i]])
	    y[i] <- y[parent[i]] - ifelse(temp2<=fudge & haskids, fudge, temp2)
	    }
	y <- y / (max(y))
        }

    # Now compute the x coordinates, by spacing out the leaves and then
    #   filling in
    x   <-  double(length(node))         #allocate, then fill it in below
    x[is.leaf] <- seq(sum(is.leaf))      # leaves at 1, 2, 3, ....
    left.child <- match(node * 2, node)
    right.child <- match(node * 2 + 1, node)

    # temp is a list of non-is.leaf, by depth
    temp <- split(seq(node)[!is.leaf], depth[!is.leaf])
    for(i in rev(temp))
            x[i] <- 0.5 * (x[left.child[i]] + x[right.child[i]])

    if (nspace < 0) return(list(x=x, y=y))

    #
    # Now we get fancy, and try to do overlapping
    #
    #  The basic algorithm is, at each node:
    #      1: get the left & right edges, by depth, for the left and
    #           right sons, of the x-coordinate spacing.
    #      2: find the minimal free spacing.  If this is >0, slide the
    #           right hand son over to the left
    #      3: report the left & right extents of the new tree up to the
    #           parent
    #   A way to visualize steps 1 and 2 is to imagine, for a given node,
    #      that the left son, with all its descendants, is drawn on a
    #      slab of wood.  The left & right edges, per level, give the
    #      width of this board.  (The board is not a rectangle, it has
    #      'stair step' edges). Do the same for the right son.  Now
    #      insert some spacers, one per level, and slide right hand
    #      board over until they touch.  Glue the boards and spacer
    #      together at that point.
    #
    #  If a node has children, its 'space' is considered to extend left
    #    and right by the amount "nspace", which accounts for space
    #    used by the arcs from this node to its children.  For
    #    horseshoe connections nspace usually is 1.
    #
    #  To make it global for a recursive function, the x coordinate list
    #    is written into frame 0.
    #
    compress <- function(me, depth) {
        lson <- me +1
	x <- x
	if (is.leaf[lson]) left <- list(left=x[lson], right=x[lson],
						depth=depth+1, sons=lson)
        else               left <- compress(me+1, depth+1)

        rson <- me + 1 + length(left$sons)        #index of right son
	if (is.leaf[rson]) right<- list(left=x[rson], right=x[rson],
						depth=depth+1, sons=rson)
	else               right<- compress(rson, depth+1)

	maxd <- max(left$depth, right$depth) - depth
        mind <- min(left$depth, right$depth) - depth

	# Find the smallest distance between the two subtrees
	#   But only over depths that they have in common
	# 1 is a minimum distance allowed
	slide <- min(right$left[1:mind] - left$right[1:mind]) -1
	if (slide >0) { # slide the right hand node to the left
	    x[right$sons] <- x[right$sons] - slide;
	    x[me] <- (x[right$sons[1]] + x[left$sons[1]])/2
#	    assign("x", x)
            x <<- x
	    }
	else slide <- 0

	# report back
        if (left$depth > right$depth) {
	    templ <- left$left
            tempr <- left$right
            tempr[1:mind] <- pmax(tempr[1:mind], right$right -slide)
	    }
        else {
	    templ <- right$left  - slide
	    tempr <- right$right - slide
	    templ[1:mind] <- pmin(templ[1:mind], left$left)
	    }

	list(left = c(x[me]- nspace*(x[me] -x[lson]), templ),
	     right= c(x[me]- nspace*(x[me] -x[rson]), tempr),
	     depth= maxd+ depth, sons=c(me, left$sons, right$sons))
	}
#    assign('compress', compress)
#    assign('x', x)
#    assign('is.leaf', is.leaf)
#    assign('nspace', nspace)

#    temp <-
    compress(1, 1)
#    x <- get('x')
#    remove(c('compress', 'x', 'is.leaf', 'nspace'))
    list(x = x, y = y)
}

# SCCS @(#)rpconvert.s	1.3 06/08/01
# Convert from the orginial style rpart object to the newer
#  style object (the changes made when user-written splits were added)
#

rpconvert <- function(x)
{
    if (!inherits(x, "rpart"))
        stop("x does not appear to be an rpart object")
    ff <- x$frame
    if (is.null(ff$splits)) {
        # this appears to be a new style one already
	warning("x not converted")
	return(x)
    }
    ff$splits <- NULL
    ff$wt <- ff$n

    xlev <- attr(x, "xlevels")
    if (length(xlev) >0) {
	zz <- as.numeric(names(xlev))
	names(xlev) <- attr(x$terms, "term.labels")[zz]
	attr(x, "xlevels") <- xlev
    }

    if (x$method=="class") {
	temp <- cbind(ff$yval, ff$yval2, ff$yprob)
	dimnames(temp) <- NULL
	ff$yval2 <- temp
	ff$yprob <- NULL
	x$frame <- ff

	temp <- rpart.class(c(1,1,2,2), NULL, wt=c(1,1,1,1))#dummy call
	x$functions <- list(summary=temp$summary, print=temp$print,
			    text = temp$text)
    }

    else if (x$method=="anova") {
	x$frame <- ff

	temp <- rpart.anova(1:5, NULL, wt=rep(1,5))#dummy call
	x$functions <- list(summary=temp$summary, text = temp$text)
    }

    else {  #either exp or poisson (they have the same summary/text pair)
	ff$yval2 <- cbind(ff$yval, ff$yval2)
	x$frame <- ff

	temp <- rpart.poisson(1:5, NULL, wt=rep(1,5))#dummy call
	x$functions <- list(summary=temp$summary, text = temp$text)
    }

    class(x) <- "rpart"
    x
}
## This function plots the approximate r-square for the different
## splits (assumes using anova method).

## SCCS @(#)rsq.rpart.s	1.6 08/28/97

rsq.rpart <- function(x) {

  if(!inherits(x,'rpart')) stop("Not legitimate rpart")

  p.rpart <- printcp(x)
  xstd <- p.rpart[,5]
  xerror <- p.rpart[,4]
  rel.error <- p.rpart[,3]
  nsplit <- p.rpart[,2]
  method <- x$method

  if(!method=='anova') cat("May not be applicable for this method\n")

  plot(nsplit, 1-rel.error, xlab='Number of Splits', ylab='R-square',
       ylim=c(0,1), type='o')
  par(new=TRUE)
  plot(nsplit, 1-xerror, type='o', ylim=c(0,1),lty=2, xlab=' ', ylab=' ')
  legend(0,1, c('Apparent','X Relative'), lty=1:2)


  ylim <- c(min(xerror-xstd) -.1, max(xerror + xstd) + .1)
  plot(nsplit, xerror, xlab='Number of Splits', ylab='X Relative Error',
       ylim=ylim, type='o')
  segments(nsplit, xerror - xstd, nsplit, xerror + xstd)
  invisible()

  }

# SCCS %W% %G%
#
#  Interactively snip off part of a tree
#

snip.rpart.mouse <- function(tree,
		      parms=paste(".rpart.parms", dev.cur(), sep = ".")) {
    xy <- rpartco(tree)
    toss <- NULL
    ff <- tree$frame
    if (exists(parms, envir=.GlobalEnv)) {
        parms <- get(parms, envir=.GlobalEnv)
	branch <- parms$branch
	}
    else branch <- 1

    node <- as.numeric(row.names(tree$frame))
    draw <- rpart.branch(xy$x,xy$y, node, branch)

    lastchoice <- 0
    while (length(choose <- identify(xy, n=1, plot=FALSE)) >0 ) {
	if (ff$var[choose] == '<leaf>') {
		cat("Terminal node -- try again\n")
		next
		}

	if (choose != lastchoice) {
	    # print out some info on the click
	    cat("node number:", node[choose], " n=", ff$n[choose], "\n")
	    cat("    response=", format(ff$yval[choose]))
	    if (is.null(ff$yval2)) cat ("\n")
	    else if (is.matrix(ff$yval2))
		  cat(" (", format(ff$yval2[choose,]), ")\n")
	    else  cat(" (", format(ff$yval2[choose]), ")\n")
	    cat("    Error (dev) = ", format(ff$dev[choose]), "\n")
	    lastchoice <- choose
	    }
	else {
	    # second click-- erase all of the descendants
	    #   (stolen from snip.tree)
	    id  <- node[choose]
	    id2 <- node
	    while (any(id2>1)) {
		id2 <- floor(id2/2)
		temp  <- (match(id2, id, nomatch=0) >0)
  	        id <- c(id, node[temp])
		id2[temp] <- 0
		}
	    temp <- match(id, node[ff$var != '<leaf>'], nomatch=0)
	    lines(c(draw$x[,temp]), c(draw$y[,temp]), col=0)
	    toss <- c(toss, node[choose])
	    }
	}
    toss
    }
#SCCS  @(#)snip.rpart.s	1.10 10/30/01
#
#  This routine "throws away" branches
#

snip.rpart <- function(x, toss) {
    if (!inherits(x, 'rpart')) stop("Not an rpart object")

    if (missing(toss) || length(toss)==0) {
        toss <- snip.rpart.mouse(x)
	if (length(toss)==0) return(x)
	}

    where <- x$where
    ff   <- x$frame
    id    <- as.numeric(row.names(ff))
    index <- ff$index
    ff.n  <- length(id)

    toss <- unique(toss)
    toss.idx <- match(toss, id, nomatch=0) #the rows of the named nodes
    if (any(toss.idx ==0)) {
	warning(paste("Nodes", toss[toss.idx==0], "are not in this tree"))
	toss <- toss[toss.idx>0]
        toss.idx <- toss.idx[toss.idx>0]
        }

#    if (any(toss==1))  {
#	# a special case that causes grief later
#	warning("Can't prune away the root node and still have a tree!")
#        return(NULL)
#	}

    # Now add all of the descendants of the selected nodes
    #   We do this be finding all node's parents.
    #        (Division by 2 gives the parent of any node.)
    #   At each step we make id2 <- parent(id2), and augment 'toss' with
    #     found children.  The loop should take <  log_2(maxdepth)/2 steps
    id2 <- id
    while (any(id2>1)) {
	id2 <- floor(id2/2)
	xx <- (match(id2, toss, nomatch=0) >0)
	toss <- c(toss, id[xx])
        id2[xx] <- 0
	}

    # Now "toss" contains all of the nodes that should not be splits
    temp <- match(floor(toss/2) , toss, nomatch=0)  #which are leaves?
    newleaf <- match(toss[temp==0], id)             # row numbers, leaves
    keepit <- (1:ff.n)[is.na(match(id,toss))]  # row numbers to be let be

    # Compute the parent row for each row in the splits structure
    #  Then "thin out" the splits and csplit components
    n.split <- rep((1:ff.n), ff$ncompete + ff$nsurrogate+ 1*(ff$var!='<leaf>'))
    split <- x$splits[match(n.split, keepit, nomatch=0) >0, ,drop=FALSE]
    temp <- split[,2] >1      #which rows point to categoricals?
    if (any(temp)) {
        x$csplit <- x$csplit[split[temp,4], , drop=FALSE]
	split[temp,4] <- 1
        if(is.matrix(x$csplit)) split[temp,4] <- 1:nrow(x$csplit)
	}
    else x$csplit <- NULL
    x$splits <- split

    # Thin out unneeded rows in the frame component
    ff$ncompete[newleaf] <- ff$nsurrogate[newleaf] <- 0
    ff$var[newleaf]     <- "<leaf>"
    x$frame <- ff[sort(c(keepit, newleaf)),]

    # Now do the 'parents' loop one more time, to fix up the "where"
    #   vector
    # This pass requires log_2(depth) iterations
    #
    id2 <- id[x$where]         #the list of old leaf nodes
    id3 <- id[sort(c(keepit, newleaf))]
    temp <- match(id2, id3, nomatch=0)
    while (any(temp==0)) {
	id2[temp==0] <- floor(id2[temp==0]/2)
	temp <- match(id2, id3, nomatch=0)
	}
    x$where <- match(id2, id3)

    x
    }
#SCCS  @(#)summary.rpart.s	1.18 07/05/01

summary.rpart <- function(object, cp=0, digits=getOption("digits"), file,  ...)
{
    if(!inherits(object, "rpart")) stop("Not legitimate rpart object")

    # If this is an older-style rpart object, convert it
    #  either way, rename it to "x" to save typing
    if (!is.null(object$frame$splits)) x <- rpconvert(object)
    else  x <- object

    if (!missing(file)) {
	  sink(file)
	  on.exit(sink())
	  }

    if(!is.null(x$call)) {
        cat("Call:\n")
        dput(x$call)
        }

    omit <- x$na.action
    n <- x$frame$n
    if (length(omit))
          cat("  n=", n[1], " (", naprint(omit), ")\n\n", sep="")
    else cat("  n=", n[1], "\n\n")

    print(x$cptable, digits=digits)
    ff <- x$frame
    ylevel <- attr(x,'ylevels')
    id <- as.integer(row.names(ff))
    parent.id <- ifelse(id==1,1, floor(id/2))
    parent.cp <- ff$complexity[match(parent.id, id)]
    rows <- (1:length(id))[parent.cp > cp]
    if (length(rows)>0) rows <- rows[order(id[rows])]
    else rows <- 1
    is.leaf <- (ff$var=='<leaf>')
    index <- cumsum(c(1, ff$ncompete + ff$nsurrogate + 1*(!is.leaf)))

    if(!all(is.leaf)) {  #skip these lines for a "no splits" tree
        sname <- dimnames(x$splits)[[1]]
        cuts <- vector(mode='character', length=nrow(x$splits))
        temp <- x$splits[ ,2]
        for (i in 1:length(cuts)) {
            if (temp[i] == -1)
                cuts[i] <-paste("<", format(signif(x$splits[i,4], digits=digits)))
            else if (temp[i] ==1)
                cuts[i] <-paste("<", format(signif(x$splits[i,4], digits=digits)))
            else cuts[i]<- paste("splits as ",
                                 paste(c("L", "-", "R")[x$csplit[x$splits[i,4], 1:temp[i]]],
                                       collapse='', sep=''), collapse='')
        }
    # S-PLUS 4.0 can't handle null vectors here
        if(any(temp<2)) cuts[temp<2 ] <- format(cuts[temp<2],justify="left")
        cuts <- paste(cuts, ifelse(temp >=2, ",",
                                   ifelse(temp==1, " to the right,", " to the left, ")),
                      sep = '')
    }

    if (is.null(ff$yval2))
        tprint <- x$functions$summary(ff$yval[rows], ff$dev[rows],
                                      ff$wt[rows], ylevel, digits)
    else
        tprint <- x$functions$summary(ff$yval2[rows,], ff$dev[rows],
                                      ff$wt[rows], ylevel, digits)

    for (ii in 1:length(rows)) {
	i <- rows[ii]
	nn <- ff$n[i]
	twt <- ff$wt[i]
	cat("\nNode number ", id[i], ": ", nn, " observations", sep='')
	if (ff$complexity[i] < cp || is.leaf[i]) cat("\n")
	else cat(",    complexity param=",
                 format(signif(ff$complexity[i], digits)), "\n", sep="")

	cat(tprint[ii], "\n")
	if (ff$complexity[i] > cp && !is.leaf[i] ){
	    sons <- 2*id[i] + c(0,1)
	    sons.n <- ff$n[match(sons, id)]
	    cat("  left son=", sons[1], " (", sons.n[1], " obs)",
		" right son=", sons[2], " (", sons.n[2], " obs)", sep='')
	    j <- nn - (sons.n[1] + sons.n[2])
	    if (j>1) cat(", ", j, " observations remain\n", sep='')
	    else if (j==1) cat(", 1 observation remains\n")
	    else     cat("\n")
	    cat("  Primary splits:\n")
	    j <- seq(index[i], length=1+ff$ncompete[i])
	    if (all(nchar(cuts[j]) < 25))
                temp <- format(cuts[j], justify="left")
	    else  temp <- cuts[j]
	    cat(paste("      ", format(sname[j], justify="left"), " ", temp,
		      " improve=", format(signif(x$splits[j,3], digits)),
		      ", (", nn - x$splits[j,1], " missing)", sep=''),
                sep="\n")
	    if (ff$nsurrogate[i] >0) {
		cat("  Surrogate splits:\n")
		j <- seq(1 +index[i] + ff$ncompete[i], length=ff$nsurrogate[i])
		agree <- x$splits[j,3]
		if (all(nchar(cuts[j]) < 25))
                    temp <- format(cuts[j], justify="left")
		else  temp <- cuts[j]
		if (ncol(x$splits)==5) {
		    adj   <- x$splits[j,5]
		    cat(paste("      ", format(sname[j], justify="left"), " ",
			      temp,
			      " agree=", format(round(agree, 3)),
			      ", adj=" , format(round(adj, 3)),
			      ", (", x$splits[j,1], " split)", sep=''),
			sep="\n")
                }
		else {                  #an older style rpart object -- no adj value present
		    cat(paste("      ", format(sname[j], justify="left"), " ",
			      temp,
			      " agree=", format(round(agree, 3)),
			      ", (", x$splits[j,1], " split)", sep=''),
			sep="\n")
                }
            }
        }
    }
    cat("\n")
    invisible(x)
}
# SCCS @(#)text.rpart.s	1.12 06/06/01
# This is a modification of text.tree.
# Fancy option has been added in (to mimic post.tree)
#

text.rpart <- function (x, splits = TRUE, which = 1, label = "yval", FUN = text, 
    all = FALSE, pretty = NULL, digits = getOption("digits") - 2, tadj = 0.65,
    stats = TRUE, use.n = FALSE, bars = TRUE, xadj = 1, yadj = 1, bord = FALSE,
    big.pts = FALSE, ...) 
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
    xy <- rpartco(x)
    node <- as.numeric(row.names(x$frame))
    is.left <- (node%%2 == 0)
    node.left <- node[is.left]
    parent <- match(node.left/2, node)
    add.bars <- bars & is.matrix(frame$yval2)
    text.adj <- ifelse(add.bars, yadj * diff(range(xy$y))/12, 
        0)
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
    leaves <- if (all) 
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
    if (add.bars) {
        bar.vals <- x$functions$bar(yval2 = frame$yval2)
        sub.barplot(xy$x, xy$y, bar.vals, leaves, xadj = xadj, 
            yadj = yadj, bord = bord, line = TRUE , col = c("lightblue", 
                "blue", "darkblue"))
        rx <- range(xy$x)
        ry <- range(xy$y)
        if (!is.null(ylevels)) 
            bar.labs <- ylevels
        else bar.labs <- dimnames(x$y)[[2]]
        legend(min(xy$x) - 0.1 * rx, max(xy$y) + 0.05 * ry, bar.labs, 
            col = c("lightblue", "blue", "darkblue"), pch = 15, 
            bty = "n", ...)
    }
    if (big.pts) 
        points(xy$x[leaves], xy$y[leaves], pch = 16, cex=3*par()$cex, col = 2:(sum(leaves) + 
            1))
    invisible()
}


# SCCS @(#)xpred.rpart.s	1.18 07/05/01
#
#  Get a set of cross-validated predictions

xpred.rpart <- function(fit, xval=10, cp)
{
    if (!inherits(fit, 'rpart')) stop("Invalid fit object")

    method <- fit$method
    method.int <- pmatch(method, c("anova", "poisson", "class", "user", "exp"))
    if (method.int==5) method.int <- 2
    Terms <- fit$terms

    Y <- fit$y
    X <- fit$x
    wt<- fit$wt
    if (is.null(Y) || is.null(X)) {
	m <- fit$model
	if (is.null(m)) {
	    m <-fit$call[match(c("", 'formula', 'data', 'weights', 'subset',
                                 'na.action'),
                               names(fit$call), nomatch=0)]
	    if (is.null(m$na.action)) m$na.action<- na.rpart
	    m[[1]] <- as.name("model.frame.default")
	    m <- eval(m, parent.frame())
        }
	if (is.null(X)) X <- rpart.matrix(m)
	if (is.null(wt)) wt <- model.extract(m, "weights")
	if (is.null(Y)) {
	    yflag <- TRUE
	    Y <- model.extract(m, "response")
            offset <- attr(Terms, "offset")
	    if (method != user) {
		init <- (get(paste("rpart", method, sep='.')))(Y,offset, NULL)
		Y <- as.matrix(init$y)
		numy <- ncol(Y)
            }
        }
	else {
	    yflag <- FALSE
	    Y <- as.matrix(Y)
	    numy <- ncol(Y)
	    offset <- 0
        }
    }
    else {
	yflag <- FALSE
	Y <- as.matrix(Y)
	numy <- ncol(Y)
	offset <- 0
    }

    nobs <- nrow(X)
    nvar <- ncol(X)
    if (length(wt)==0) wt <- rep(1.0, nobs)

    cats <- rep(0, nvar)
    xlevels <- attr(fit, "xlevels")
    if (!is.null(xlevels)){
        cats[match(names(xlevels), dimnames(X)[[2]])] <-
            unlist(lapply(xlevels, length))
    }

    controls <- fit$control
    if (missing(cp)) {
	cp<- fit$cptable[,1]
	cp <- sqrt(cp * c(10, cp[-length(cp)]))
	cp[1] <- (1+fit$cptable[1,1])/2
    }
    ncp <- length(cp)

    if (length(xval)==1) {
                                        # make random groups
	xgroups <- sample(rep(1:xval, length=nobs), nobs, replace=FALSE)
    }
    else if (length(xval) == nrow(Y)) {
	xgroups <- xval
	xval <- length(unique(xgroups))
    }
    else {
        ## Check to see if observations were removed due to missing
	if (!is.null(fit$na.action)) {
            ## if na.rpart was used, then na.action will be a vector
	    temp <- as.integer(fit$na.action)
	    xval <- xval[-temp]
	    if (length(xval) == nobs) {
		xgroups <- xval
		xval <- length(unique(xgroups))
            }
	    else stop("Wrong length for xval")
        }
	else stop("Wrong length for xval")
    }

    costs <- fit$call$costs
    if (is.null(costs)) costs <- rep(1.0, nvar)

    parms <- fit$parms
    if (method=='user') {
	mlist <- fit$functions
	if (length(parms)==0) init <- mlist$init(Y, offset, wt=wt)
	else                  init <- mlist$init(Y, offset, parms, wt)

        ## assign this to avoid garbage collection
        keep <- rpartcallback(mlist, nobs, init)
    }

    rpfit <- .C("s_xpred",
                n = as.integer(nobs),
                nvarx = as.integer(nvar),
                ncat = as.integer(cats * !fit$ordered),
                method= as.integer(method.int),
                as.double(unlist(controls)),
                parms = as.double(unlist(parms)),
                as.integer(xval),
                as.integer(xgroups),
                as.double(t(Y)),
                as.double(X),
                as.integer(is.na(X)),
                pred = double(ncp* nobs),
                as.integer(ncp),
                as.double(cp * fit$frame[1,"dev"]),
                error = character(1),
                wt = as.double(wt),
                as.integer(numy),
                as.double(costs),
                NAOK=TRUE, PACKAGE = "mvpart")
    if (rpfit$n == -1)  stop(rpfit$error)

    matrix(rpfit$pred, ncol=ncp, byrow=TRUE,
           dimnames=list(dimnames(X)[[1]], format(cp)) )
}

.First.lib <- function(lib, pkg) {
	library.dynam("mvpart", pkg, lib)
	cat(" mvpart package loaded: extends rpart to include\n") 
	cat(" multivariate and distance-based partitioning\n\n")
}

.noGenerics <- TRUE

tree.depth <- function (nodes)
{
    depth <- floor(log(nodes, base = 2) + 1e-7)
    as.vector(depth - min(depth))
}

string.bounding.box <- function(s)
{
    s2 <- strsplit(s, "\n")
    rows <- sapply(s2, length)
    columns <- sapply(s2, function(x) max(nchar(x)))
    list(columns=columns, rows=rows)
}

node.match <- function(nodes, nodelist, leaves, print.it = TRUE)
{
    node.index <- match(nodes, nodelist, nomatch = 0)
    bad <- nodes[node.index == 0]
    if(length(bad) > 0 & print.it)
        warning(paste("supplied nodes", paste(bad, collapse = ","),
                      "are not in this tree"))
    good <- nodes[node.index > 0]
    if(!missing(leaves) && any(leaves <- leaves[node.index])) {
        warning(paste("supplied nodes",
                      paste(good[leaves], collapse = ","), "are leaves"))
        node.index[node.index > 0][!leaves]
    }
    else node.index[node.index > 0]
}

descendants <- function(nodes, include = TRUE)
{
    n <- length(nodes)
    if(n == 1) return(matrix(TRUE, 1, 1))
    ind <- 1:n
    desc <- matrix(FALSE, n, n)
    if(include) diag(desc) <- TRUE
    parents <- match((nodes %/% 2), nodes)
    lev <- floor(log(nodes, base = 2))
    desc[1, 2:n] <- TRUE
    for(i in max(lev):2) {
        desc[cbind(ind[parents[lev == i]], ind[lev == i])] <- TRUE
        parents[lev == i] <- parents[parents[lev == i]]
        lev[lev == i] <- i - 1
    }
    return(desc)
}

rpart.pca <- function (tree, pts = TRUE, plt.allx = TRUE, speclab = TRUE,
    specvecs = TRUE, wgt.ave = FALSE, add.tree = TRUE,  
    cv1 = 1, cv2 = 2, interact = FALSE, ...) 
{
    if (tree$method != "mrt") 
        stop("Only for multivariate trees !! \n")
    if (nrow(tree$frame) < 4) 
        stop("Only 2 terminal nodes -- PCA not done !! \n")
    old.par <- par(mar = rep(2, 4))
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
    xall <- z$u[, 1:maxd, drop = FALSE ] %*% d
    x <- amat %*% (z$v)[, 1:maxd, drop = FALSE ]
    xlv <- x[leaves, ]
    if (!wgt.ave) 
        y <- z$v[, 1:maxd, drop = FALSE ]
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
    ppts <- rep("Ö", length(treegrps))
    xx <- (scale(as.matrix(data), center = TRUE , scale = FALSE ) %*% 
        z$v)[, 1:maxd, drop = FALSE ]
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
    barpts <- rep("Ö", ln)
    myy <- sqrt(apply(y[, c(cv1, cv2)]^2, 1, sum))
    sc <- ifelse(wgt.ave, 1, max(mxx)/max(myy))
    repeat {
        eqscplt(c(sc * y[, cv1], x[, cv1]), c(sc * y[, cv2], 
            x[, cv2]), axes = FALSE , xlab = "", ylab = "", type = "n", 
            tol = 0.15, sym = FALSE )
        points(0, 0, pch = 4, col = 7, cex = 2)
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
        if (plt.allx) {
            unitg <- sort(unique(treegrps))
            for (i in 1:length(unitg)) points(xx[unitg[i] == 
                treegrps, cv1], xx[unitg[i] == treegrps, cv2], 
                pch = 16, col = i + 1, cex = par()$cex)
        }
        if (pts) {
            lvnode <- sort(node[leaves])
            for (i in 1:length(lvnode)) points(xlv[, cv1][lvnode[i] == 
                lvnode], xlv[, cv2][lvnode[i] == lvnode], pch = 16, 
                cex = 2 * par()$cex, col = i + 1)
        }
        if (add.tree) {
            pp <- match(c(even.node, even.node + 1), node)
            nn <- length(even.node)
            from <- pp[1:nn]
            to <- pp[(nn + 1):(2 * nn)]
            segments(x[from, cv1], x[from, cv2], x[to, cv1], 
                x[to, cv2])
        }
        if (speclab) 
            text(sc * y[, cv1], sc * y[, cv2] + 0.5 * adj * specvecs * 
                (y[, cv2] > 0), specs, col = "black", cex = par()$cex)
        points(0, 0, pch = 3, cex = par()$cex * 3, col = 1)
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
    invisible()
}


gdist <- function(x, method = "bray", keepdiag = FALSE , full = FALSE, sq = FALSE)
{
	METHODS <- c("manhattan", "euclidean", "canberra", "bray", "kulczynski", "gower", 
	"maximum", "binary", "chisq", "chord")
	method <- pmatch(method, METHODS)
	if(is.na(method))
		stop("invalid distance method")
	N <- nrow(x <- as.matrix(x))
	if(method == 6) {
		mx <- apply(x, 2, max)
		mn <- apply(x, 2, min)
		x <- x/rep(mx - mn, N)
		method <- 2
	}
	else if(method == 9) {
		rr <- apply(x, 1, sum)
		cc <- apply(x, 2, sum)
		x <- diag(1/sqrt(rr)) %*% x %*% diag(1/sqrt(cc))
		method <- 2
	}
	else if(method == 10) {
		mns <- sqrt(apply(x^2, 1, sum))
		x <- x/(mns * sqrt(2))
		method <- 2
	}
	d <- .C("gdistance",
		x = as.double(x),
		nr = N,
		nc = ncol(x),
		d = double((N * (N - 1))/2),
		keepdiag = as.integer(FALSE),
		method = as.integer(method))$d
	attr(d, "Size") <- N
	class(d) <- "dist"
	if (full) d <- distfull(d)
	if (sq) d <- d^2
	d
}


xdiss <- function(data, dcrit = 1, dauto = TRUE , dinf = 0.5, method = "man", use.min = TRUE, 
		eps = 0.0001, replace.neg = TRUE, big = 10000, sumry = TRUE, full = FALSE, sq = FALSE)
{ 
		scale.row <- function(data, p = 1) 		{
		tmp <- apply(data, 1, sum, na.rm = TRUE )
		if(any(t0 <- (tmp == 0)))
		cat(sum(t0), " rows with sum = 0 !!  -- these rows untransformed\n")
		if(p == 1)
			data[!t0,  ] <- data[!t0,  ]/apply(data[!t0,  ], 1, sum, na.rm = TRUE )
		else if(p == 2)
			data[!t0,  ] <- data[!t0,  ]/(apply(data[!t0,  ]^2, 1, sum, na.rm = TRUE ))^0.5
		data
		}	

		METHODS <- c("manhattan", "euclidean", "canberra", "bray", "kulczynski", "gower", 
		"maximum", "binary", "chisq", "chord")
		method <- METHODS[pmatch(method, METHODS)]
		if(is.na(method))
		stop("invalid distance method")
        if(any(data < 0))
                data <- apply(data, 2, function(x)
                x - min(x))
        n <- dim(data)[1]
        if(method == "chisq" | method == "gower" | method == "maximum") {
                method <- "manhattan"
                cat("This dissimilarity is not suitable -- switching to Manhattan\n")
        }
        if(method == "manhattan") {
                cat("Using Extended Dissimilarity : Manhattan (Site Standardised by Mean)\n")
                data <- scale.row(as.matrix(data), p = 1)/2
                d <- gdist(data, method = "man")
        }
        else if(method == "chord") {
                cat("Using Extended Dissimilarity : Chord \n")
                data <- scale.row(as.matrix(data), p = 2)/sqrt(2)
                d <- gdist(data, method = "euc")
        }
        else if(method == "euclidean") {
                cat("Using Extended Dissimilarity : Euc (Site Standardised by SS) \n")
                data <- scale.row(as.matrix(data), p = 2)/sqrt(2)
                d <- gdist(data, method = "euc")
        }
        else if(method == "bray") {
                cat("Using Extended Dissimilarity : Bray \n")
                d <- gdist(data, method = "bra")
        }
        else if(method == "canberra") {
                cat("Using Extended Dissimilarity : Canberra \n")
                d <- gdist(data, method = "can")
        }
        else if(method == "binary") {
                cat("Using Extended Dissimilarity : Binary \n")
                d <- gdist(data, method = "bin")
        }
        else if(method == "kulczynski") {
                cat("Using Extended Dissimilarity : Kulczynski \n")
                d <- gdist(data, method = "kul")
        }
        cat("Maximum distance = ", round(max(d), 4), "\n")
        if(dauto) {
                dcrit <- max(apply(distfull(d) + diag(1, n), 1, min))
                dcrit <- dcrit * (1 - dinf) + dinf
                cat("Critical distance = ", signif(dcrit, 4), "\n")
        }
        cat("% Distances > Crit Dist = ", round(100 * mean(d > dcrit), 2), "\n")
        use.min <- ifelse(use.min, 1, 0)
        storage.mode(d) <- "double"
        storage.mode(n) <- "integer"
        storage.mode(dcrit) <- "double"
        storage.mode(use.min) <- "integer"
        storage.mode(eps) <- "double"
        storage.mode(big) <- "double"
        dnew <- .C("xdists",
                d = d,
                n,
                dcrit,
                use.min,
                eps,
                big)$d
        if(any(dnew == -1, na.rm = TRUE ))
                attr(dnew, "ok") <- F
        else attr(dnew, "ok") <- TRUE 
        if(any(dnew == -1))
                cat("WARNING : Data disconnected\n")
        if(replace.neg)
                dnew[dnew == -1] <- max(dnew)
        if(sumry) {
                cat("Summary of Extended Dissimilarities\n")
                print(summary(dnew))
        }
	    attr(dnew, "Size") <- n
        class(dnew) <- "dist"
        if (full) dnew <- distfull(dnew)
		if (sq) dnew <- dnew^2
        dnew
}

scaler <- function(x, col = NULL, row = NULL)
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


cmds.diss <- function (data, k = ncol(data), x.use = FALSE, method = "man", 
    zero.chk = TRUE, plt = FALSE, plot.subset = FALSE , plot.subn = 5) 
{
    if (x.use) {
        xdists <- xdiss(data, method = method, dauto = TRUE , dcrit = 0.6)
        xds <- cmdscale(xdists, k = k)
		colnames(xds) <- paste("s",1:ncol(xds),sep="")
            }
    else {
        xdists <- gdist(data, method = method)
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

distfull <- function(dis)
{
	n <- attr(dis, "Size")
	full <- matrix(0, n, n)
	full[lower.tri(full)] <- dis
	full + t(full)
}

trclcomp <- function(x, method = "com")
{
    if(class(x) != "rpart")
        stop("Rpart object needed!")
    if(x$method != "mrt")
        stop("Multivariate tree needed!")
    cpt <- x$cptable
    size <- cpt[, 2] + 1
    nr <- nrow(cpt)
    tree.err <- cpt[, 3]
    clust.err <- rep(1, nr)
    sst <- sum(scale(x$y, scale = FALSE)^2)
    d <- dist(x$y)
    if(any(is.na(d))) {
        cat("Warning -- NA distances in cluster -- replacing by 0\n")
        d[is.na(d)] <- 0
    }
        hclout <- hclust(d, method = method)
    for(i in 2:nr) {
        grp <- factor(cutree(hclout, k = size[i]))
        clust.err[i] <- sum(resid(lm(x$y ~ factor(grp), sin = TRUE))^2)/sst
        NULL
    }
    minerr <- min(c(tree.err, clust.err))
    plot(size, tree.err, type = "n", ylim = c(minerr, 1), xlab = "Size", ylab = "Resubstition Error")
    points(size, tree.err, type = "o", col = 2, pch = 16)
    points(size, clust.err, type = "o", col = 4, pch = 16)
    legend(mean(size), 1, c("Tree", "Cluster"), col = c(2, 4), lty = 1, bty = "n")
    title("Comparison of tree and cluster errors across size")
    cat("Tree error     : ", signif(tree.err, 3), "\n")
    cat("Cluster error  : ", signif(clust.err, 3), "\n")
	invisible(list(tree.err = tree.err, clust.err = clust.err))
    }

