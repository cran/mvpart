\name{trclcomp}
\alias{trclcomp}
\title{Tree-Clustering Comparison}
\description{
  This function compares the within-group variation for groups formed by
  tree partitioning and unconstrained clustering.
  The results are plotted and returned invisibly.
}
\usage{
trclcomp(x, method = "com")
}

\arguments{
  \item{x}{ Rpart object with method "mrt" -- see \code{\link{rpart}}}
  \item{method}{ The clustering method for the unconstrained clustering
  -- see \code{\link[mva]{hclust}}}.
}

\details{
  The within-group variation for groups formed by tree partitioning and 
  unconstrained clustering are compared for all sizes of the hierarchy of tree
  partitions,   
}

\value{
  Returns a list (invisibly) of the within-tree and within-cluster variation for all tree sizes.
}

\references{ 

De'ath G. (2002)  
Multivariate Regression Trees : A New Technique for Constrained Classification Analysis.
Ecology 83(4):1103-1117.
  
  \author{ Glenn De'ath }
}
  
\examples{
data(spider)
fit <- mvpart(data.matrix(spider[,1:12])~herbs+reft+moss+sand+twigs+water,spider)  
trclcomp(fit)
}

\keyword{ multivariate }%-- one or more ...
