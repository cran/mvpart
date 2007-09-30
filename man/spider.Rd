\name{spider}
\alias{spider}
\title{Spider Data}
\description{
  The \code{spider} data frame has 28 rows and 18 columns. 
  The first 12 columns are abundances of different species of spiders
  and the next 6 are environmental data.
}
\usage{
data(spider)
}
\format{
This data frame contains the following columns:
\describe{
\item{\code{arct.lute, pard.lugu, zora.spin, pard.nigr, pard.pull,
 aulo.albi, troc.terr, alop.cune, pard.mont, alop.acce,
 alop.fabr, arct.peri }}{
numeric vectors giving the abundances of 12 species of spider
}
\item{\code{water, sand, moss, reft, twigs, herbs}}{
numeric vectors giving the values of 6 environmental characteristics
}}
}

\source{
   Van der Aart, P. J. and N. Smeeck-Enserink. 1975. 
   Correlations between distributions of hunting spiders (Lycosidae, Ctenidae) 
   and environmental characteristics in a dune area.  
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
