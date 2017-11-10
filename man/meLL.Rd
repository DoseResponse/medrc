\name{meLL}

\alias{meLL.5}
\alias{meLL.4}
\alias{meLL.3}
\alias{meLL.2}

\title{The log-logistic function}

\description{
  'meLL.X' provides the X-parameter log-logistic function.
}

\usage{
  meLL.5(x, b, c, d, e, f)

  meLL.4(x, b, c, d, e)
  
  meLL.3(x, b, d, e)
  
  meLL.2(x, b, e)
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
  The five-parameter logistic function is given by the expression
  \deqn{ f(x) = c + \frac{d-c}{(1+\exp(b(\log(x)-\log(e))))^f}}
  
  or in another parameterisation
  \deqn{ f(x) = c + \frac{d-c}{(1+\exp(b(\log(x)-e)))^f}}
  
  The function is asymmetric for \eqn{f} different from 1.
  
}

\references{
  Finney, D. J. (1979) Bioassay and the Practise of Statistical Inference, 
  \emph{Int. Statist. Rev.}, \bold{47}, 1--12.
}

\author{Daniel Gerhard}


\keyword{models}
\keyword{nonlinear}
