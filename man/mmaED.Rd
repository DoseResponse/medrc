\name{mmaED}
\alias{mmaED}
\title{ED estimation for averaged medrc and glsdrc models}
\description{Estimates ED for averaged fixed effects of several medrc or glsdrc objects}
\usage{
mmaED(object, ..., respLev, ic = c("AIC", "BIC"), interval = c("none", "buckland", "kang"), 
     level = 0.95, marginal=FALSE, nGQ=5, rfinterval=c(0, 1000))
}
\arguments{
  \item{object}{an medrc object}
  \item{\dots}{further medrc objects}
  \item{respLev}{a numeric vector containing the response levels}
  \item{ic}{the information criterion used for model averaging}
  \item{interval}{character string specifying the type of confidence
    intervals to be supplied.}
  \item{level}{confidence level, default at 0.95}
  \item{marginal}{logical. If TRUE, marginalized ED estimates are obtained. See EDmarg for further information.}
  \item{nGQ}{integer. Specifies the number nof nodes for Gauss-Hermite quadrature.}
  \item{rfinterval}{numeric vector. Interval for root finding (uniroot) to search for ED values. Default is [0, 1000].}
}

\author{Christian Ritz, Daniel Gerhard}

\keyword{ models }
\keyword{ htest }
