\name{cmds.diss}
\alias{cmds.diss}
\title{Classical Scaling of Dissimilarity Measures}
\description{
  The function first computes the dissimilarity matrix according to the specified method --
  see \code{\link[mvpart]{gdist}} or \code{\link[mvpart]{xdiss}}. The dissimilarities are then
  scaled using classical scaling -- see \code{\link[stats]{cmdscale}}. The returned matrix can be input into
  \code{\link[mvpart]{rpart}} or \code{\link[mvpart]{mvpart}} for multivariate regression tree splitting.

}
\usage{
cmds.diss(data, k = ncol(data), x.use = FALSE, zero.chk = TRUE,
            plt = FALSE, plot.subset = FALSE, plot.subn = 5, ...)
}

\arguments{
  \item{data}{ Data matrix}
  \item{k}{ Number of vectors to be returned}
  \item{x.use}{ Use extended dissimilarity?}
  \item{zero.chk}{ Check for zero row sums -- if zero ignore these rows according to method}
  \item{plt}{ Plot the relationship between the dissimilarities and the distances calculated
  from the scaled output vectors. }
  \item{plot.subset}{ Plot a subset of the points -- useful for large data sets. }
  \item{plot.subn}{ Controls how many points are plotted when \code{plot.subset}=TRUE.
  The number of points plotted is 750 + N * plot.subn where N = number of rows in \code{data}.}
  \item{\dots}{ arguments passed to either \code{xdiss} or \code{gdist}}
}

\details{
  The function knows the same dissimilarity indices as \code{\link[mvpart]{gdist}}.
  Plotting the relationship between the dissimilarities and the distances calculated
  from the scaled output vectors is useful in assessing potential loss of information.
  If the loss is high then the results from partitioning directly from the dissimilarity
  matrix using distance-base partitioning (see \code{dist} in \code{\link[mvpart]{rpart}}),and those
  obtained from partitioning the output of cmds.diss using multivariate regression trees
  (see \code{mrt} in \code{\link[mvpart]{rpart}}) can be substantial.

}

\author{ Glenn De'ath }

\examples{
data(spider)
dist.vecs <- cmds.diss(spider)

# comparing splitting using "dist" and "mrt" methods
# for euclidean distance the answers are indentical :
# first using "mrt" on the data directly
mvpart(data.matrix(spider[,1:12])~water+twigs+reft+herbs+moss+sand,
spider,method="mrt",size=5)

# now using "dist" -- note we need the full distance matrix squared
mvpart(gdist(spider[,1:12],meth="euc",full=TRUE,sq=TRUE)~water+twigs
+reft+herbs+moss+sand,spider,method="dist",size=5)

# finally using "mrt" from the scaled dissimilarities.
mvpart(cmds.diss(spider[,1:12],meth="euc")~water+twigs+reft
+herbs+moss+sand,spider,method="mrt",size=5)

# try with some other measure of dissimilarity
# eg extended bray-curtis -- the result will differ
# between methods

}

\keyword{ multivariate }%-- one or more ...
