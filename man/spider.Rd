\name{spider}
\alias{spider}
\title{Spider Data}

\description{
Data set on abundances of spiders and environmental predictors.
All variables are rated on a 0-9 scale.
}
\usage{data(spider)}
\format{
  A data frame with 28 observations with 12 species and
  6 environmental predictors.
}

\source{
   Van der Aart, P. J. and N. Smeeck-Enserink. 1975.
   Correlations between distributions of hunting spiders
   (Lycosidae, Ctenidae) and environmental characteristics
   in a dune area.
   Netherlands Journal of Zoology. 25:1-45.

   These data were analysed using multivariate trees in
   De'ath, G. 2002. Multivariate Regression Trees:
   A New Technique for Modelling Species-Environment Relationships.
   Ecology. 83(4):1103-1117
}

\examples{
data(spider)
fit<-mvpart(as.matrix(spider[,1:12])~water+twigs+reft+herbs+moss+sand,spider)
summary(fit)
}
\keyword{datasets}
