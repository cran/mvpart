"xdiss" <-
function(data, dcrit = 1, dauto = TRUE , dinf = 0.5, method = "man", use.min = TRUE, 
        eps = 0.0001, replace.neg = TRUE, big = 10000, sumry = TRUE, full = FALSE, sq = FALSE)
{ 
        scale.row <- function(data, p = 1)      {
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
                big,
                PACKAGE="mvpart")$d
        if(any(dnew == -1, na.rm = TRUE ))
                attr(dnew, "ok") <- FALSE
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

