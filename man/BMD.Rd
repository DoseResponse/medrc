\name{BMD}
\alias{BMD}
\title{BMD estimation for medrc or glsdrc models}
\description{Estimates the benchmark dose for fixed effects of medrc or glsdrc objects}
\usage{
BMD(object, respLev, interval = c("none", "delta", "fls", "tfls"), 
    level = 0.95, bmd = c("additional", "extra"), 
    background = 0.05, display=TRUE)
}
\arguments{
  \item{object}{an medrc object}
  \item{respLev}{a numeric vector containing the benchmark response levels}
  \item{interval}{character string specifying the type of confidence
    intervals to be supplied.}
  \item{level}{confidence level, default at 0.95}
  \item{bmd}{benchmark dose estimation (smallest dose resulting in a
    probability of an abnormal response)}
  \item{background}{probability of an abnormal response}
  \item{display}{logical. If TRUE results are displayed; otherwise they are not}
}

\author{Christian Ritz, Daniel Gerhard}

\keyword{ models }
\keyword{ htest }
