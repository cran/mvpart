\name{text.rpart}
\alias{text.rpart}
\title{
  Place Text on a Dendrogram 
}
\description{
  Labels the current plot of the tree dendrogram with text.
}
\usage{
\method{text}{rpart}(x, splits = TRUE, which = 4, label = "yval", FUN = text, 
    all.leaves = FALSE, pretty = NULL, digits = getOption("digits") - 2,
    tadj = 0.65, stats = TRUE, use.n = FALSE, bars = TRUE, 
    legend = FALSE, xadj = 1, yadj = 1, bord = FALSE, big.pts = FALSE,
    uniform = FALSE, branch = 1, nspace = -1, minbranch = 0.3, \dots)
}

\arguments{
\item{x}{
  fitted model object of class \code{rpart}.  This is assumed to be the result
  of some function that produces an object with the same named
  components as that returned by the \code{rpart} function.
}
\item{splits}{
  logical flag.  If \code{TRUE} (default), then the splits in the tree are
  labeled with the criterion for the split.
}
\item{which}{
  labels splits 1 = center, 2 = left, 3 = right , 4 = both.  
}
\item{label}{
a column name of \code{x$frame};  values of this
will label the nodes.  For the \code{"class"} method,
\code{label="yval"} results in the factor levels
being used, \code{"yprob"} results in the probability of the winning factor
level being used, and 'specific yval level' results in the probability
of that factor level.   
}
\item{FUN}{
  the name of a labeling function, e.g. \code{text}.  
}
\item{all.leaves}{
  Logical.  If \code{TRUE}, all nodes are labeled, otherwise just
  terminal nodes.
}
\item{tadj}{
Adjustment of text above (or below) splits.
}

\item{pretty}{
  an integer denoting the extent to which factor levels in split labels
  will be abbreviated.  A value of (0) signifies no abbreviation.  A
  \code{NULL}, the default, signifies using elements of letters to represent
  the different factor levels.
}
\item{digits}{
  number of significant digits to include in numerical labels.
}
\item{stats}{
  If \code{TRUE} adds statistics to nodes.
}
\item{use.n}{
  If \code{TRUE} adds N to labels. 
  (\#events level1/ \#events level2/etc. for \code{class},
  \code{n} for \code{anova}, and \#events/n for \code{poisson} and \code{exp}).
}
\item{bars}{
  If \code{TRUE} adds barplots for multivariate regression trees.
}
\item{legend}{
  If \code{TRUE} adds legends for multivariate regression trees.
}
\item{xadj, yadj}{
 varies the size of barplots for multivariate regression trees.
}
\item{bord}{
 Adds borders (boxes) to barplots for multivariate regression trees.
}
\item{big.pts}{
 Adds color coded points to nodes. Useful to track groups to PCA plot
(see \code{rpart.pca}).
}
\item{uniform}{
 Uniform spacing of tree branches; defualt is FALSE.
}
\item{branch}{
 branch parameter default = 1.
}
\item{nspace}{
 branch parameter default = 1.
}
\item{minbranch}{
 minimum branch parameter default = 0.3.
}
\item{\dots}{
  Graphical parameters may also be supplied as arguments to this
  function (see \code{par}).
}}
\section{Side Effects}{
  the current plot of a tree dendrogram is labeled.
}
\seealso{
  \code{\link{text}}, \code{\link{plot.rpart}}, \code{\link{rpart}},
  \code{\link{post.rpart}}, \code{\link{abbreviate}} 
}
\examples{
data(car.test.frame)
z.auto <- rpart(Mileage ~ Weight, car.test.frame)
plot(z.auto)
text(z.auto, use.n=TRUE, all=TRUE)
}
\keyword{tree}
% docclass is function
% Converted by Sd2Rd version 1.20.
