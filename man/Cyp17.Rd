\name{Cyp17}
\alias{Cyp17}
\docType{data}
\title{Cyp17 expression data}
\description{Observed expression measured at several dose levels for multiple runs.}
\usage{data(Cyp17)}
\format{
  A data frame with 63 observations on the following 3 variables.
  \describe{
    \item{\code{run}}{ID of 3 different runs}
    \item{\code{dose}}{5 dose levels (0, 0.1, 10, 100, 500)}
    \item{\code{expression}}{observed expression}
  }
}


\examples{
data(Cyp17)
# transforming the vector of observed expression 
Cyp17$logexpression <- log(Cyp17$expression) + 5

# assuming a 4-parameter log-logistic model
# upper asymptote fixed at 2.59
# random run effects on the lower asymptote and the ED50
Cyp17.LL4.mixed <- medrm(logexpression ~ dose, 
                         random = c + e ~ 1 | run, 
                         fct = LL.4(fixed = c(NA,NA,-2.41 + 5,NA)), 
                         data = Cyp17, 
                         start=c(-0.5, -4 + 5, 700)) 
  
plot(Cyp17.LL4.mixed, logx=TRUE) + 
    theme_classic()
plot(Cyp17.LL4.mixed, logx=TRUE, ranef=TRUE) + 
    theme_classic()
    
summary(Cyp17.LL4.mixed)
}



\keyword{datasets}
