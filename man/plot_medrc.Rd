\name{plot.medrc}
\alias{plot.medrc}

\title{Plot medrc objects}
\description{Plot predicted dose-response curves for medrc objects}
\usage{
\method{plot}{medrc}(x, ..., ndose=25, ranef=FALSE, level=NULL, logx = FALSE)
}

\arguments{
  \item{x}{An object of class medrc}
  \item{ndose}{Number of points to interpolate the dose response curve}
  \item{ranef}{If FALSE the predictions are based only on fixed effects, if TRUE the random effects are presented as additional lines around the fixed effects curve. Only implemented for curveid=NULL.}
  \item{level}{Hierarchical level of random effects, which are presented when ranef=TRUE.}
  \item{logx}{If TRUE, plot x-axis on a logarithmic scale}
  \item{...}{further arguments}
}

\author{Daniel Gerhard}



\keyword{ methods }
