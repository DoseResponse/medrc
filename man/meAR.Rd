\name{meAR}

\alias{meAR.2}
\alias{meAR.3}

\title{Asymptotic regression model}

\description{
  Providing the mean function and the corresponding self starter function for the asymptotic regression model.
}

\usage{
  meAR.2(x, d, e)

  meAR.3(x, c, d, e)
}

\arguments{
  \item{x}{numeric vector specifying the dose.}
  \item{c}{lower limit at at \eqn{x=0}}
  \item{d}{upper limit}
  \item{e}{\eqn{e>0} is determining the steepness of the increase as \eqn{x}}
}

\details{
  The asymptotic regression model is a three-parameter model with mean function:
  
  \deqn{ f(x) = c + (d-c)(1-\exp(-x/e))}
  
  The parameter \eqn{c} is the lower limit (at \eqn{x=0}), the parameter \eqn{d} is the upper limit 
  and the parameter \eqn{e>0} is determining the steepness of the increase as \eqn{x}.
}

\value{
  A list of class \code{drcMean}, containing the mean function, the self starter function,
  the parameter names and other components such as derivatives and a function for calculating ED values.
}

%\references{ ~put references to the literature/web site here ~ }

\author{Daniel Gerhard}

\keyword{models}
\keyword{nonlinear}

