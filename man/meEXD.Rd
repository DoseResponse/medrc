\name{meEXD}

\alias{meEXD.2}
\alias{meEXD.3}

\title{Exponential decay model}

\description{
  Exponential decay model with or without a nonzero lower limit.
}

\usage{
  meEXD.2(x, d, e)

  meEXD.3(x, c, d, e)
}

\arguments{
  \item{x}{numeric dose vector.}
  \item{c}{lower limit}
  \item{d}{upper limit}
  \item{e}{steepness}
}

\details{
  The exponential decay model is a three-parameter model with mean function:
  
  \deqn{f(x) = c + (d-c)(\exp(-x/e))}
  
  The parameter init is the upper limit (attained at \eqn{x=0}), the parameter plateau is the lower limit 
  reached for x going to infinity and the parameter \eqn{e>0} is determining the steepness of the 
  decay. The curve is monotonously decreasing in \eqn{x}.
}


\references{
  Organisation for Economic Co-operation and Development (OECD) (2006)
  \emph{Current approaches in the statistical analysis of ecotoxicity data: A guidance to application - annexes},
  Paris: OECD (p. 80).
}

\author{Daniel Gerhard}

%\note{}

\keyword{models}
\keyword{nonlinear}
