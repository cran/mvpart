"rpart" <-
function (formula, data = NULL, weights, subset, na.action = na.rpart,
    method, dissim, model = FALSE, x = FALSE, y = TRUE, parms, control,
    cost, ...)
{
    call <- match.call()
    if (is.data.frame(model)) {
        m <- model
        model <- FALSE
    }
    else {
        m <- match.call(expand.dots = FALSE)
        m$model <- m$method <- m$control <- NULL
        m$x <- m$y <- m$parms <- m$... <- NULL
        m$cost <- m$dissim <- NULL
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
### added 15/04/13
    	parms <- init$parms
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

    if(missing(dissim) || is.null(dissim))  dissim <- "euclidean"
    dissim.int <- pmatch(dissim, c("euclidean", "manhattan"))
    if(is.na(dissim.int))
        stop("Invalid dissimilarity")
    dissim <- c("euclidean", "manhattan")[dissim.int]
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
        as.double(X), as.integer(dissim.int), as.integer(!is.finite(X)),
        error = character(1), wt = as.double(wt), as.integer(init$numy), as.double(cost),
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
            2], ncompete = pmax(0, rp$inode[, 3] - 1), nsurrogate = rp$inode[, 4])
    }
    if (method == "class") {
        numclass <- init$numresp - 1
        temp <- rp$dnode[, -(1:4), drop = FALSE] %*% diag(init$parms$prior *
            sum(init$counts)/pmax(1, init$counts))
        yprob <- temp/rowSums(temp)
        yval2 <- matrix(rp$dnode[, -(1:3), drop = FALSE], ncol = numclass +
            1)
        frame$yval2 <- cbind(yval2, yprob)
    }
    else if (method == "mrt") {
        frame$yval <- apply(rp$dnode[, -c(1:3), drop = FALSE], 1, mean)
        frame$yval2 <- rp$dnode[, -(1:3), drop = FALSE]
    }
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
            terms = Terms, cptable = t(rp$cptable), method = method, dissim = dissim,
            parms = init$parms, control = controls, functions = functions)
    }
    else {
        ans <- list(frame = frame, where = where, call = call,
            terms = Terms, cptable = t(rp$cptable), splits = splits,
            method = method, dissim = dissim, parms = init$parms, control = controls,
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
