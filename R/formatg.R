"formatg" <-
function(x, digits= unlist(options('digits')),
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
    else matrix(temp,nrow=1)
    }

