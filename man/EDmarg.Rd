\name{EDmarg}

\alias{EDmarg}


\title{Estimating marginal effective doses for medrc objects}

\description{

  \code{ED} estimates marginal effective doses (ECp/EDp/ICp) for given reponse levels, conditional on the estimated variance components.
}

\usage{
  EDmarg(object, respLev, 
         interval = c("none", "delta", "fls", "tfls"), 
         clevel=NULL,
         level = ifelse(!(interval == "none"), 0.95, NULL), 
         reference = c("control", "upper"), 
         type = c("relative", "absolute"), nGQ=5, 
         rfinterval=c(0, 1000), lref, uref, bound = TRUE, 
         display = TRUE, logBase = NULL, ...)
}

\arguments{
  \item{object}{an object of class 'medrc'.}
  \item{respLev}{a numeric vector containing the response levels.}
  \item{interval}{character string specifying the type of confidence intervals to be supplied. The default is "none". 
  Use "delta" for asymptotics-based confidence intervals (using the delta method and the t-distribution).
  Use "fls" for from logarithm scale based confidence intervals (in case the parameter in the model is log(ED50) as for
  the \code{\link{llogistic2}}) models. The only alternative for model-robust fits is using inverse regression.}
  \item{clevel}{character string specifying the curve id in case on estimates for a specific curve or compound is requested. By default estimates are shown for all curves.}
  \item{level}{numeric. The level for the confidence intervals. The default is 0.95.}
  \item{reference}{character string. Is the upper limit or the control level the reference?}
  \item{type}{character string. Whether the specified response levels are absolute or relative (default).}
  \item{nGQ}{integer. Specifies the number nof nodes for Gauss-Hermite quadrature.}
  \item{rfinterval}{numeric vector. Interval for root finding (uniroot) to search for ED values. Default is [0, 1000].}
  \item{lref}{numeric value specifying the lower limit to serve reference.}
  \item{uref}{numeric value specifying the upper limit to serve reference (eg. 100\%).}  
  \item{bound}{logical. If TRUE only ED values between 0 and 100\% are allowed. FALSE is useful for hormesis models.}
  \item{display}{logical. If TRUE results are displayed. Otherwise they are not (useful in simulations).}
  \item{logBase}{numeric. The base of the logarithm in case logarithm transformed dose values are used.}  
  \item{...}{see the details section below.}
}


\value{
  A matrix with two or more columns, containing the estimates 
  and the corresponding estimated standard errors and possibly lower and upper confidence limits.
}

%\references{ ~put references to the literature/web site here ~ }
\author{Christian Ritz, Daniel Gerhard}

%\note{This function is only implemented for the built-in functions of class 'braincousens', 'gompertz', 'logistic' and 'mlogistic'.}

\keyword{models}
\keyword{nonlinear}
