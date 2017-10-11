\name{BMDadjresp}
\alias{BMDadjresp}
\title{Adjusted ED response levels for BMD estimation based on medrc or glsdrc models}
\description{Calculates adjusted response levels for estimation of the BMD}
\usage{
BMDadjresp(object, respLev, bmd = c("additional", "extra"), background = 0.05)
}
\arguments{
  \item{object}{an medrc object}
  \item{respLev}{a numeric vector containing the benchmark response levels}
  \item{bmd}{benchmark dose estimation (smallest dose resulting in a
    probability of an abnormal response)}
  \item{background}{probability of an abnormal response}
}

\author{Christian Ritz, Daniel Gerhard}

\keyword{ models }
\keyword{ htest }
