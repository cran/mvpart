\name{rpart.dist}
\alias{rpart.dist}
\title{
  Set up for `dist' Method for an Rpart Model
}
\description{
  Set up for \code{dist} method for an \code{rpart} model.
}
\usage{
rpart.dist(y, offset, parms, wt)
}
\arguments{
  \item{y}{the responses.}
  \item{offset}{an offset, or \code{NULL}.}
  \item{parms}{a list of parameters, usually empty.}
  \item{wt}{case weights -- ignored.}
}
\value{
  Returns a list consisting of
  \item{y}{(subtracting \code{offset} if necessary),}
  \item{parms}{as input,}
  \item{numresp}{the number of responses,}
  \item{summary}{a function to be invoked by \code{\link{summary.rpart}},}
  \item{text}{a function to be invoked by \code{\link{text.rpart}}.}
}
\details{
  \code{rpart.dist} is to be used only by the function \code{rpart}.
}
\keyword{internal}

