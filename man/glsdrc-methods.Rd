\name{glsdrc-methods}
\alias{df.residual.glsdrc}
\alias{predict.glsdrc}
\alias{vcov.glsdrc}
\alias{summary.glsdrc}
\alias{residuals.glsdrc}
\alias{AIC.glsdrc}
\alias{BIC.glsdrc}
\alias{logLik.glsdrc}
\alias{print.glsdrc}

\title{gdrc methods}
\description{Methods for objects of class gdrc}
\usage{
\method{df.residual}{glsdrc}(object, ...)
\method{predict}{glsdrc}(object, ..., newdata = NULL)
\method{residuals}{glsdrc}(object, ...)
\method{vcov}{glsdrc}(object, ...)
\method{summary}{glsdrc}(object, ...)
\method{AIC}{glsdrc}(object, ..., k = 2)
\method{BIC}{glsdrc}(object, ...)
\method{logLik}{glsdrc}(object, REML = FALSE, ...)
\method{print}{glsdrc}(x, ..., digits = max(3, getOption("digits") - 3))
}

\arguments{
  \item{object}{An object of class glsdrc}
  \item{x}{An object of class glsdrc}
  \item{newdata}{an optional data frame to be used for obtaining the predictions. All variables used in the fixed and random effects models, as well as the grouping factors, must be present in the data frame. If missing, the fitted values are returned.}
  \item{k}{numeric, the penalty per parameter to be used; the default k = 2 is the classical AIC.}
  \item{REML}{an optional logical value. If TRUE the restricted log-likelihood is returned, else, if FALSE, the log-likelihood is returned. Defaults to FALSE.}  
  \item{digits}{minimal number of significant digits}
  \item{...}{further arguments}
}

\author{Daniel Gerhard}



\keyword{ methods }
