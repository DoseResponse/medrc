\name{BMDmarg}
\alias{BMDmarg}
\title{Marginal BMD estimation for medrc models}
\description{Estimates the marginal benchmark dose, conditional on the estimated variance components.}
\usage{
BMDmarg(object, respLev, 
        interval = c("none", "delta", "fls", "tfls"), 
        clevel=NULL, 
        level = ifelse(!(interval == "none"), 0.95, NULL), 
        bmd = c("additional", "extra"), 
        background = 0.05, nGQ=5, rfinterval=c(0, 1000), 
        display=TRUE, ...)
}
\arguments{
  \item{object}{an medrc object}
  \item{respLev}{a numeric vector containing the benchmark response levels}
  \item{interval}{character string specifying the type of confidence
    intervals to be supplied.}
  \item{clevel}{character string specifying the curve id in case on estimates for a specific curve or compound is requested. By default estimates are shown for all curves.}
  \item{level}{confidence level, default at 0.95}
  \item{bmd}{benchmark dose estimation (smallest dose resulting in a
    probability of an abnormal response)}
  \item{background}{probability of an abnormal response}
  \item{nGQ}{integer. Specifies the number nof nodes for Gauss-Hermite quadrature.}
  \item{rfinterval}{numeric vector. Interval for root finding (uniroot) to search for ED values. Default is [0, 1000].}
  \item{display}{logical. If TRUE results are displayed; otherwise they are not}
  \item{...}{additional arguments}
}

\author{Christian Ritz, Daniel Gerhard}

\keyword{ models }
\keyword{ htest }
