\name{mvpart}
\alias{mvpart}
\title{
Recursive Partitioning and Regression Trees
}
\description{
Wrapper function for fitting and plotting \code{rpart} models
}

\usage{
mvpart(form, data, minauto = TRUE, size, xv = c("1se", "min",
    "pick", "none"), xval = 10, xvmult = 0, xvse = 1, snip = FALSE,
    plot.add = TRUE, text.add = TRUE, digits = 3, margin = 0,
    uniform = FALSE, which = 4, pretty = TRUE, use.n = TRUE,
    all.leaves = FALSE, bars = TRUE, legend, bord = FALSE,
    xadj = 1, yadj = 1, prn = FALSE, branch = 1, rsq = FALSE,
    big.pts = FALSE, pca = FALSE, interact.pca = FALSE,
    wgt.ave.pca = FALSE, keep.y = TRUE, \dots)
}

\arguments{
\item{form}{
As for \code{rpart} function.
Arguments to rpart can be passed by \dots.
}
\item{data}{
Optional data frame in which to interpret the variables named in the
formula
}
\item{minauto}{
If \code{TRUE} uses smart minsplit and minbucket based on N cases.
}
\item{size}{
The size of tree to be generated.
}
\item{xv}{
Selection of tree by cross-validation: \code{"1se"} - gives best tree within one SE of the overall best,
\code{"min"} - the best tree, \code{"pick"} - pick the tree size interactively,
\code{"none"} - no cross-validation.
}
\item{xval}{
Number of cross-validations or vector defining cross-validation groups.
}

\item{xvmult}{
Number of multiple cross-validations.
}
\item{xvse}{
Multiplier for the number of SEs used for \code{xv = "1se"}.
}
\item{plot.add}{
Plot the tree and (optionally) add text.
}
\item{text.add}{
Add output of \code{text.rpart} to tree.
}
\item{snip}{
Interactively prune the tree.
}
\item{digits}{
Number of digits on labels.
}
\item{margin}{
Margin around plot, 0.1 gives an extra 10 percent space around the plot.
}
\item{uniform}{
Uniform lengths to the branches of the tree.
}
\item{which}{
Which split labels and where to plot them, 1=centered, 2 = left, 3 = right and 4 = both.
}
\item{pretty}{
Pretty labels or full labels.
}
\item{use.n}{
Add number of cases at each node.
}
\item{all.leaves}{
Annotate all nodes.
}
\item{bars}{
If \code{TRUE} adds barplots to nodes.
}
\item{legend}{
If \code{TRUE} adds legend for mrt and classification trees.
}
\item{bord}{
Border (box) around the barplots.
}
\item{xadj, yadj}{
Adjust the size of the individual barplots (default = 1).
}
\item{prn}{
If \code{TRUE} prints tree details.
}
\item{branch}{
Controls spread of branches: 1=vertical lines, 0=maximum slope.
}
\item{rsq}{
If \code{TRUE} gives "rsq" plot.
}
\item{big.pts}{
Plot colored points at leaves -- useful to link to PCA plot.
}
\item{pca}{
If \code{TRUE} plots PCA of group means and add species and site information.
}
\item{interact.pca}{
If \code{TRUE} runs interactive PCA. See \code{rpart.pca}.
}
\item{wgt.ave.pca}{
If \code{TRUE} plot weighted averages acorss sites for species.
}
\item{keep.y}{
If \code{TRUE} y values are returned.
}
\item{\dots}{
\dots other arguments passed to \code{rpart}.
}
}

\value{
an object of class \code{rpart}, a superset of class tree.
}

\seealso{
  \code{\link{rpart}}, \code{\link{rpart.pca}},
}
\examples{
data(spider)
mvpart(data.matrix(spider[,1:12])~herbs+reft+moss+sand+twigs+
water,spider)       # defaults
mvpart(data.matrix(spider[,1:12])~herbs+reft+moss+sand+twigs+
water,spider,xv="p")  # pick the tree size
# pick cv size and do PCA
fit <- mvpart(data.matrix(spider[,1:12])~herbs+reft+moss+sand+
twigs+water,spider,xv="1se",pca=TRUE)
rpart.pca(fit,interact=TRUE,wgt.ave=TRUE)
# interactive PCA plot of saved multivariate tree
}
\keyword{tree}
