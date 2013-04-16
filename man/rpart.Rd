\name{rpart}
\alias{rpart}
\alias{rpartcallback}
\title{
Recursive Partitioning and Regression Trees
}
\description{
Fit a \code{rpart} model
}
\usage{
rpart(formula, data=NULL, weights, subset, na.action=na.rpart, method,
      dissim, model=FALSE, x=FALSE, y=TRUE, parms, control, cost, \dots)
}
\arguments{
\item{formula}{
a formula, as in the \code{lm} function.
}
\item{data}{
an optional data frame in which to interpret the variables named in the
formula
}
\item{weights}{
optional case weights.
}
\item{subset}{
optional expression saying that only a subset of the rows of the data
should be used in the fit.
}
\item{na.action}{
The default action deletes all observations for which \code{y} is missing,
but keeps those in which one or more predictors are missing.
}
\item{method}{
one of \code{"anova"}, \code{"poisson"}, \code{"class"}, \code{"mrt"},
\code{"dist"}, or \code{"exp"}.
If \code{method} is missing then the routine tries to make an intellegent guess.
If \code{y} is a survival object, then \code{method="exp"} is assumed,
if \code{y} is a matrix then \code{method="mrt"} is assumed,
if \code{y} is a factor then \code{method="class"} is assumed, otherwise \code{method="anova"}
is assumed.  It is wisest to specifiy the method directly, especially as
more criteria are added to the function.

For method=\code{"dist"} the response must be a square symmetric distance matrix; e.g.
returned by \code{\link[mvpart]{gdist}} or \code{\link[mvpart]{xdiss}}. Weights and
cross-validation are currently not implemented for method=\code{"dist"}.

Alternatively, \code{method} can be a list of functions named
\code{init}, \code{split} and \code{eval}.
}
\item{dissim}{
used when method=\code{"anova"} or method=\code{"mrt"}. Dissimilarity types are either
\code{"euc"} for Euclidean (sums of squares about the mean) or \code{"man"} for Manhattan
(sums of absolute deviations about the mean. The latter is experimental and has proved useful
for ecological data.
}
\item{model}{
keep a copy of the model frame in the result.
If the input value for \code{model} is a model frame (likely from an earlier
call to the \code{rpart} function), then this frame is used rather than constructing
new data.
}
\item{x}{
keep a copy of the \code{x} matrix in the result.
}
\item{y}{
keep a copy of the dependent variable in the result.
}
\item{parms}{
optional parameters for the splitting function.
Anova splitting has no parameters.
Poisson splitting has a single parameter, the coefficient of variation of
the prior distribution on the rates.  The default value is 1.
Exponential splitting has the same parameter as Poisson.
For classification splitting, the list can contain any of:
the vector of prior probabilities (component \code{prior}), the loss matrix
(component \code{loss}) or the splitting index (component \code{split}).  The
priors must be positive and sum to 1.  The loss matrix must have zeros
on the diagnoal and positive off-diagonal elements.  The splitting
index can be \code{gini} or \code{information}.  The default priors are
proportional to the data counts, the losses default to 1,
and the split defaults to \code{gini}.
}
\item{control}{
options that control details of the \code{rpart} algorithm.
}
\item{cost}{
a vector of non-negative costs, one for each variable in the
model. Defaults to one for all variables.  These are scalings to be
applied when considering splits, so the improvement on splitting on a
variable is divided by its cost in deciding which split to choose.
}
\item{\dots}{
arguments to \code{rpart.control} may also be specified in the call to
\code{rpart}.  They are checked against the list of valid arguments.
}
}
\value{
an object of class \code{rpart}, a superset of class tree.
}
\details{
This differs from the tree function mainly in its handling of surrogate
variables.  In most details it follows Breiman et. al. quite closely.
}
\references{
Breiman, Friedman, Olshen, and Stone. (1984)
\emph{Classification and Regression Trees.}
Wadsworth.

De'ath G. (2002)
Multivariate Regression Trees : A New Technique for Constrained Classification Analysis.
Ecology 83(4):1103-1117.
}
\seealso{
  \code{\link{rpart.control}}, \code{\link{rpart.object}},
  \code{\link{summary.rpart}}, \code{\link{print.rpart}}
}
\examples{

data(car.test.frame)
z.auto <- rpart(Mileage ~ Weight, car.test.frame)
summary(z.auto)
plot(z.auto); text(z.auto)

data(spider)
fit1 <- rpart(data.matrix(spider[,1:12])~water+twigs+reft+herbs+
moss+sand,spider,method="mrt")
plot(fit1); text(fit1)
fit2 <- rpart(data.matrix(spider[,1:12])~water+twigs+reft+herbs+
moss+sand,spider,method="mrt",dissim="man")
plot(fit2); text(fit2)
fit3 <- rpart(gdist(spider[,1:12],meth="bray",full=TRUE,sq=TRUE)
~water+twigs+reft+herbs+moss+sand,spider,method="dist")
plot(fit3); text(fit3)
}
\keyword{tree}
