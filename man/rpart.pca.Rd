\name{rpart.pca}
\alias{rpart.pca}
\title{
Principle Components Plot of a Multivariate Rpart Object
}
\description{
Plots a PCA of the rpart object on the current graphics device.
}
\usage{
\method{rpart}{pca}(tree, pts = TRUE, plt.allx = TRUE, speclabs = TRUE,
    specvecs = TRUE, wgt.ave = FALSE, add.tree = TRUE,  
    cv1 = 1, cv2 = 2, chulls = TRUE, interact = FALSE, \dots) 
}
\arguments{
\item{tree}{
A fitted object of class \code{rpart} containing a multivariate 
regression tree.
}
\item{pts}{
If \code{TRUE}, large points representing the leaf means are plotted.
}
\item{plt.allx}{
If \code{TRUE}, small points representing individual cases are plotted.
}
\item{speclabs}{
If \code{TRUE} the labels of the response variables are plotted.
}
\item{specvecs}{
If \code{TRUE} the vectors of the response variables are plotted provided
\code{wgt.ave} is \code{FALSE} .
}
\item{wgt.ave}{
If \code{TRUE} use weighted averages of responses not vectors.
}
\item{add.tree}{
If \code{TRUE} add the tree structure to the plot.
}
\item{cv1}{
Defines the principal component to plot horizontally -- but see \code{interact}.
}
\item{cv2}{
Defines the principal component to plot vertically -- but see \code{interact}.
}
\item{chulls}{
If \code{TRUE} adds convex hulls to thr tree groups.
}
\item{interact}{
If \code{TRUE} the plot can be viewed in dimensions by left-clicking to 
top-left, bottom-right or bottom-left (reset).
}
\item{\dots}{
  arguments to be passed to or from other methods.
}}
\value{
NULL
}
\section{Side Effects}{
A PCA biplot plot is produced on the current graphics device.
}
\details{
This function plots a PCA biplot of the group means (leaves) of
 multivariate regresssion objects of class \code{rpart}.
The responses and group means and indivdual cases can be shown on the plot.
If responses are positive (eg species-environment data) weighted averages of 
responses can be plotted.

}
\seealso{
\code{\link{rpart}}, \code{\link{text.rpart}}
}
\examples{
data(spider)
fit<-mvpart(data.matrix(spider[,1:12])~herbs+reft+moss+sand+twigs+water,spider)
rpart.pca(fit)
rpart.pca(fit,wgt.ave=TRUE,interact=TRUE)
}
\keyword{tree}
% docclass is function
% Converted by Sd2Rd version 1.20.
