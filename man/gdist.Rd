\name{gdist}
\alias{gdist}
\title{Dissimilarity Measures}
\description{
  The function computes useful dissimilarity indices which are known
  to have a good rank-order relation with gradient separation and are
  thus efficient in community ordination with multidimensional scaling.
}
\usage{
gdist(x, method="bray", keepdiag=FALSE, full=FALSE, sq=FALSE)
}

\arguments{
  \item{x}{ Data matrix}
  \item{method}{ Dissimilarity index }
  \item{keepdiag}{ Compute amd keep diagonals}
  \item{full}{ Return the square dissimilarity matrix}
  \item{sq}{ Square the dissimilarities -- useful for distance-based partitioning}
}

\details{
  The function knows the following dissimilarity indices:
  \tabular{ll}{
    \code{euclidean}
    \tab \eqn{d_{jk} = \sqrt{\sum_i (x_{ij}-x_{ik})^2}}{d[jk] = sqrt(sum (x[ij]-x[ik])^2)}
    \cr
    \code{manhattan}
    \tab \eqn{d_{jk} = \sum_i |x_{ij} - x_{ik}|}{d[jk] = sum(abs(x[ij] -
      x[ik]))}
    \cr
    \code{gower}
    \tab \eqn{d_{jk} = \sum_i \frac{|x_{ij}-x_{ik}|}{\max_i-\min_i}}{d[jk] = sum (abs(x[ij]-x[ik])/(max(i)-min(i))}
    \cr
    \code{canberra}
    \tab \eqn{d_{jk}=\frac{1}{N-Z} \sum_i
      \frac{|x_{ij}-x_{ik}|}{x_{ij}+x_{ik}}}{d[jk] = (1/NZ) sum ((x[ij]-x[ik])/(x[ij]+x[ik]))}
    \cr
    \code{bray}
    \tab \eqn{d_{jk} = \frac{\sum_i |x_{ij}-x_{ik}|}{\sum_i (x_{ij}+x_{ik})}}{d[jk] = (sum abs(x[ij]-x[ik])/(sum (x[ij]+x[ik]))}
    \cr
    \code{kulczynski}
    \tab \eqn{d_{jk} = 1-0.5(\frac{\sum_i \min(x_{ij},x_{ik})}{\sum_i x_{ij}} +
      \frac{\sum_i \min(x_{ij},x_{ik})}{\sum_i x_{ik}} )}{d[jk] 1 - 0.5*((sum min(x[ij],x[ik])/(sum x[ij]) + (sum
      min(x[ij],x[ik])/(sum x[ik]))}
	\cr
  	\code{maximum}
    \tab \eqn{d_{jk} = \max_i |x_{ij} - x_{ik}|}{d[jk] = max(abs(x[ij] -
      x[ik]))}
    \cr
	\code{binary}
    \tab \eqn{d_{jk} = \sum_i |x_{ij}>0 - x_{ik}>0|}{d[jk] = max(abs(x[ij]>0 -
      x[ik]>0))}
    \cr
    \code{chord}
    \tab \eqn{d_{jk} = \sqrt{\sum_i (x_{ij}-x_{ik})^2} / {\sum_i (x_{ij}+x_{ik})^2}}{d[jk] = sqrt((sum (x[ij]-x[ik])^2)/(sum (x[ij]+x[ik])^2))}
    \cr
    }
    
  where \eqn{N-Z}{NZ} is the number of non-zero entries.

  Infamous ''double zeros'' are removed in Canberra dissimilarity.

  Euclidean and Manhattan dissimilarities are not good in gradient
  separation without proper standardization but are still included for
  comparison and special needs.

  Some of indices become identical or rank-order similar after some
  standardizations. 
}
\value{
  Should be interchangeable with \code{\link[stats]{dist}} and returns a
  distance object of the same type.
}
\references{ Faith, D.P, Minchin, P.R. and Belbin, L. (1987)
  Compositional dissimilarity as a robust measure of ecological
  distance. \emph{Vegetatio} 69, 57-68. }

\author{ Jari Oksanen  -- modified Glenn De'ath (Dec 03)}

\note{The  function is an alternative to \code{\link[stats]{dist}} adding
  some ecologically meaningful indices.  Both methods should produce
  similar types of objects which can be interchanged in any method
  accepting either.  Manhattan and Euclidean dissimilarities should be
  identical in both methods, and Canberra dissimilary may be similar.
}

\examples{
data(spider)
spider.dist <- gdist(spider[1:12,])
}
\keyword{ multivariate }%-- one or more ...
