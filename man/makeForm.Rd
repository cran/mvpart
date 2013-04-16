\name{makeform}
\alias{makeform}
\title{Formula Constructor}
\description{
  This function constructs a formula from a data frame given the locations of the
  response and other variables. It is useful for long or repetitive formulas.
}

\usage{
makeform(data, ycol, xcol, zcol, FUN, maxy = 20, extra)
}

\arguments{
  \item{data}{ The data frame or matrix for which the formula is to be contructed, }
  \item{ycol}{ The locations (or names) of the response variables in \code{data}. }
  \item{xcol}{ The locations (or names) of the explanatory variables in \code{data}. }
  \item{zcol}{ The locations (or names) of the conditioning variables in \code{data};
   for use in the \code{vegan} package. }
  \item{FUN}{ A function to apply to the response.}
  \item{maxy}{ The number of multivariate responses before abbreviated notation is used.
   see examples below.}
  \item{extra}{ Extra term(s) for the RHS of the formula.}
}

\details{
makeform constructs a formula and is useful for long and/or repetitive formulae.
See examples below.
}

\value{
  Returns a formula.
}

  
\examples{
data(spider)
makeform(spider,1,13:18)
# arct.lute ~ water + sand + moss + reft + twigs + herbs

makeform(spider,1:12,13:18)
# cbind(arct.lute, pard.lugu, zora.spin, pard.nigr, pard.pull, 
# aulo.albi, troc.terr, alop.cune, pard.mont, alop.acce, alop.fabr, 
# arct.peri) ~ water + sand + moss + reft + twigs + herbs

makeform(spider,1:12,13:15,16:18)
# cbind(arct.lute, pard.lugu, zora.spin, pard.nigr, pard.pull, 
# aulo.albi, troc.terr, alop.cune, pard.mont, alop.acce, alop.fabr, 
# arct.peri) ~ water + sand + moss + Condition(reft + twigs + herbs)

makeform(spider,1:12,13:15,maxy=6)
# as.matrix(spider[, 1:12]) ~ water + sand + moss

makeform(spider,1:3,13:15,FUN=sqrt)
# sqrt(cbind(arct.lute, pard.lugu, zora.spin)) ~ water + sand + moss

makeform(spider,1:3,13:15,FUN=sqrt,extra="I(water^2)+")
# sqrt(cbind(arct.lute, pard.lugu, zora.spin)) ~ I(water^2) + water + 
# sand + moss

}
\keyword{models}

