\name{xdiss}
\alias{xdiss}
\title{Extendend Dissimilarity Measures}
\description{
  The function computes extended dissimilarity indices which are for long gradients
  have better good rank-order relation with gradient separation and are
  thus efficient in community ordination with multidimensional scaling.
}
\usage{
xdiss(data, dcrit = 1, dauto = TRUE, dinf = 0.5, method = "man", 
    use.min = TRUE, eps = 1e-04, replace.neg = TRUE, big = 10000,
    sumry = TRUE, full = FALSE, sq = FALSE) 
}

\arguments{
  \item{data}{ Data matrix}
  \item{dcrit}{ Dissimilarities < \code{dcrit} are considered to have no species in common
   and are recalculated.}
  \item{dauto}{ Automatically select tuning parameters -- recommended.}
  \item{method}{ Dissimilarity index }
  \item{use.min}{ Minimum dissimilarity of pairs of distances used -- recommended.}
  \item{dinf, eps, replace.neg, big}{ Internal parameters -- leave as is usually.}
  \item{sumry}{ Print summary of extended dissimilarities? }
  \item{full}{ Return the square dissimilarity matrix. }
  \item{sq}{ Square the dissimilarities -- useful for distance-based partitionong. }
}
\details{
  The function knows the same dissimilarity indices as \code{\link[mvpart]{gdist}}.
    
}

\value{
  Returns an object of class distance with attributes "Size" and "ok".
  "ok" is TRUE if rows are not disconnected (De'ath 1999).
}
\references{ 
  De'ath, G. (1999) 
  Extended dissimilarity: a method of robust estimation of ecological distances
  from high beta diversity data. \emph{Plant Ecology} 144(2):191-199.

  Faith, D.P, Minchin, P.R. and Belbin, L. (1987)
  Compositional dissimilarity as a robust measure of ecological
  distance. \emph{Vegetatio} 69, 57-68. }
\author{ Glenn De'ath }

\examples{
data(spider)
spider.dist <- xdiss(spider)
}
\keyword{ multivariate }%-- one or more ...
