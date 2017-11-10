\name{meW}

\alias{meW1.4}
\alias{meW2.4}
\alias{meW1.3}
\alias{meW2.3}
\alias{meW1.2}
\alias{meW2.2}

\title{The four-parameter Weibull functions}

\description{
  'W1.X' and 'W2.X' provide the X-parameter Weibull functions.   
}

\usage{
  meW1.4(x, b, c, d, e)
  
  meW1.3(x, b, d, e)
  
  meW1.2(x, b, e)
  
  meW2.4(x, b, c, d, e)
  
  meW2.3(x, b, d, e)
  
  meW2.2(x, b, e)
}

\arguments{
  \item{x}{numeric dose vector.}
  \item{b}{steepness}
  \item{c}{lower limit}
  \item{d}{upper limit}
  \item{e}{ED50}  
}

\details{
  As pointed out in Seber and Wild (1989), there exist two different parameterisations of the Weibull model. They
  do not yield the same fitted curve for a given dataset (see under Examples).

  The four-parameter Weibull type 1 model ('weibull1') is
  \deqn{ f(x) = c + (d-c) \exp(-\exp(b(\log(x)-\log(e)))).}
 
  Thw four-parameter Weibull type 2 model ('weibull2') is
  \deqn{ f(x) = c + (d-c) (1 - \exp(-\exp(b(\log(x)-\log(e))))).}  
  
  Both four-parameter model functions are asymmetric with inflection point at the dose equal \eqn{e}. 
}



\references{
  Seber, G. A. F. and Wild, C. J (1989) \emph{Nonlinear Regression}, New York: Wiley \& Sons (pp. 330--331).
  
  Ritz, C (2009) Towards a unified approach to dose-response modeling in ecotoxicology \emph{To appear in Environ Toxicol Chem}.
}

\author{Daniel Gerhard}


\keyword{models}
\keyword{nonlinear}
