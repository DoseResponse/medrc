\name{medrc-methods}
\alias{df.residual.medrc}
\alias{predict.medrc}
\alias{vcov.medrc}
\alias{summary.medrc}
\alias{residuals.medrc}
\alias{ranef.medrc}
\alias{VarCorr.medrc}
\alias{AIC.medrc}
\alias{BIC.medrc}
\alias{logLik.medrc}
\alias{print.medrc}

\title{medrc methods}
\description{Methods for objects of class medrc}
\usage{
\method{df.residual}{medrc}(object, ...)
\method{predict}{medrc}(object, ..., newdata = NULL, level = NULL, type=c("conditional", "marginal"), nGQ=5)
\method{residuals}{medrc}(object, ...)
\method{vcov}{medrc}(object, ...)
\method{summary}{medrc}(object, ...)
\method{ranef}{medrc}(object, ...)
\method{VarCorr}{medrc}(x, sigma = 1, ...)
\method{AIC}{medrc}(object, ..., k = 2)
\method{BIC}{medrc}(object, ...)
\method{logLik}{medrc}(object, REML = FALSE, ...)
\method{print}{medrc}(x, ..., digits = max(3, getOption("digits") - 3))
}

\arguments{
  \item{object}{An object of class medrc}
  \item{x}{An object of class medrc}
  \item{newdata}{an optional data frame to be used for obtaining the predictions. All variables used in the fixed and random effects models, as well as the grouping factors, must be present in the data frame. If missing, the fitted values are returned.}
  \item{level}{an optional integer vector giving the level(s) of grouping to be used in obtaining the predictions. Level values increase from outermost to innermost grouping, with level zero corresponding to the population predictions. Defaults to the highest or innermost level of grouping.}
  \item{type}{character, either "conditional" or "marginal" for predictions conditional on the random effects or marginal predictions (conditional on the estimated variance components).}
  \item{nGQ}{integer. Number of nodes for Gauss-Hermite quadrature, when type = "marginal".}
  \item{k}{numeric, the penalty per parameter to be used; the default k = 2 is the classical AIC.}
  \item{REML}{an optional logical value. If TRUE the restricted log-likelihood is returned, else, if FALSE, the log-likelihood is returned. Defaults to FALSE.}
  \item{sigma}{an optional numeric value used as a multiplier for the standard deviations. Default is 1.}
  \item{digits}{minimal number of significant digits}
  \item{...}{further arguments}
}

\author{Daniel Gerhard}



\keyword{ methods }
