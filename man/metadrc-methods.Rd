\name{metadrc-methods}
\alias{df.residual.metadrc}
\alias{vcov.metadrc}
\alias{coef.metadrc}
\alias{coefficients.metadrc}
\alias{print.metadrc}
\alias{summary.metadrc}
\alias{confint.metadrc}

\title{metadrc methods}
\description{Methods for objects of class metadrc}
\usage{
\method{df.residual}{metadrc}(object, ...)
\method{vcov}{metadrc}(object, ...)
\method{coef}{metadrc}(object, ...)
\method{coefficients}{metadrc}(object, ...)
\method{print}{metadrc}(x, digits = max(3, getOption("digits") - 3), ...)
\method{summary}{metadrc}(object, digits = max(3, getOption("digits") - 3), ...)
\method{confint}{metadrc}(object, parm, level=0.95, fixed=TRUE, random=FALSE, ...)
}

\arguments{
  \item{object}{An object of class metadrc}
  \item{x}{An object of class metadrc}
  \item{parm}{parm argument in confint function}
  \item{digits}{An integer giving the number of digits of the parameter coefficients. Default is 3.}
  \item{level}{confidence level, default 0.95}
  \item{fixed}{if TRUE, confidence limits for population parameters are returned}
  \item{random}{if TRUE, confidence limits for variance parameters are returned, see metafor package}  
  \item{...}{further arguments}
}

\author{Daniel Gerhard}



\keyword{ methods }
