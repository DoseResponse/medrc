\name{mmplot}
\alias{mmplot}

\title{Plot multiple medrc objects}
\description{Plot multiple predicted dose-response curves based on fixed effect estimates from multiple medrc objects}
\usage{
mmplot(x, ..., ndose=25, logx = FALSE)
}

\arguments{
  \item{x}{An object of class medrc}
    \item{...}{further objects of class medrc}
  \item{ndose}{Number of points to interpolate the dose response curve}
  \item{logx}{If TRUE, plot x-axis on a logarithmic scale}
}

\author{Daniel Gerhard}

\keyword{ graphics }
