\name{meLN}

\alias{meLN.2}
\alias{meLN.3}
\alias{meLN.4}
\alias{meLLN.4}

\title{Log-normal dose-response model}

\description{
  \code{LN.X} and the accompanying convenience functions provide a general framework for specifying 
  the mean function of the decreasing or incresing log-normal dose-response model.
}

\usage{
  meLN.4(x, b, c, d, e)

  meLLN.4(x, b, c, d, e)
  
  meLN.3(x, b, d, e)
  
  meLN.2(x, b, e)
}

\arguments{
  \item{x}{numeric dose vector.}
  \item{b}{steepness}
  \item{c}{lower limit}
  \item{d}{upper limit}
  \item{e}{ED50}
}

\details{
  For the case where log(ED50), denoted \eqn{e} in the equation below, is a parameter in the model, 
  the mean function is:
  
  \deqn{f(x) = c + (d-c)(\Phi(b(\log(x)-e)))}

  and the mean function is:
  
  \deqn{f(x) = c + (d-c)(\Phi(b(\log(x)-\log(e))))}  

  in case ED50, which is also denoted \eqn{e}, is a parameter in the model. If the former model is fitted
  any estimated ED values will need to be back-transformed subsequently in order to obtain effective doses
  on the original scale.

  The mean functions above yield the same models as those described by Bruce and Versteeg (1992),
  but using a different parameterisation (among other things the natural logarithm is used).
  
  For the case \eqn{c=0} and \eqn{d=1}, the log-normal model reduces the classic probit model (Finney, 1971) 
  with log dose as explanatory variable (mostly used for quantal data). This special case is available through
  the convenience function \code{meLN.2}.
  
  The case \eqn{c=0} is available as the function \code{meLN.3}. The full four-parameter model is available
  through \code{meLN.4}.
}


\references{
  Finney, D. J. (1971) \emph{Probit analysis}, London: Cambridge University Press.

  Bruce, R. D. and Versteeg, D. J. (1992) A statistical procedure for modeling continuous toxicity data,
  \emph{Environ. Toxicol. Chem.}, \bold{11}, 1485--1494.
}

\author{Daniel Gerhard}



\keyword{models}
\keyword{nonlinear}
