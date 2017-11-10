\name{meL}

\alias{meL.2}
\alias{meL.3}
\alias{meL.4}
\alias{meL.5}


\title{The logistic model}

\description{
  The general asymmetric five-parameter logistic model for describing dose-response relationships.
}

\usage{
  meL.5(x, b, c, d, e, f)

  meL.4(x, b, c, d, e)
  
  meL.3(x, b, d, e)
  
  meL.2(x, b, e)
}


\arguments{
  \item{x}{numeric dose vector.}
  \item{b}{steepness}
  \item{c}{lower limit}
  \item{d}{upper limit}
  \item{e}{ED50}
  \item{f}{asymmetry} 
}

\details{
  The default arguments yields the five-parameter logistic mean function given by the expression
  
  \deqn{ f(x) = c + \frac{d-c}{(1+\exp(b(x - e)))^f}}
 
  The model is sometimes referred to as the Boltzmann model.
}


\author{Daniel Gerhard}

\keyword{models}
\keyword{nonlinear}
