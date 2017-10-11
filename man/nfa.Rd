\name{nfa}
\alias{nfa}
\docType{data}
\title{Network Formation Assay Data}
\description{Neurotoxicity test using a network formation assay studying
  the inhibition of network formation at acrylamide exposure.}
\usage{data(nfa)}
\format{
  A data frame with 45 observations on the following 4 variables.
  \describe{
    \item{\code{chip}}{chip ID}
    \item{\code{conc}}{7 concentrations of acrylamide, ranging from
      0-5mM}
    \item{\code{experiment}}{factor with levels 1 or 2 denoting two
      consecutive experiments}
    \item{\code{response}}{Number of connections [\%]}
  }
}
\references{
  Frimat, JP, Sisnaiske, J, Subbiah, S, Menne, H, Godoy, P, Lampen, P,
  Leist, M, Franzke, J, Hengstler, JG, van Thriel, C, West, J. The
  network formation assay: a spatially standardized neurite outgrowth
  analytical display for neurotoxicity screening. Lab Chip 2010; 10:701-709.
}

\examples{
\dontrun{
data(nfa)

ggplot(nfa, aes(x=log(conc), y=response, colour=experiment)) + 
  geom_point()

# starting values for fixed effects
fix <- coefficients(drm(response ~ conc, fct=LL.4(), data=nfa))

# starting values for random experiment effects
rexp <- drm(response ~ conc, curveid=experiment, fct=LL.4(), data=nfa)
cmatexp <- matrix(coefficients(rexp), ncol=4)
mexp <- apply(cmatexp, 2, mean)
names(mexp) <- letters[2:5]
rmatexp <- t(apply(cmatexp, 1, function(x) x-mexp))
rownames(rmatexp) <- levels(nfa$experiment)
colnames(rmatexp) <- letters[2:5]

start <- list(fixed=fix, random=list(experiment=rmatexp[,3, drop=FALSE]))


## set of nonlinear mixed models
nfa.LL4.mixed <- medrm(response ~ conc, fct=LL.4(),
                       data=nfa, random=d ~ 1 | experiment, start=start)
nfa.LN4.mixed <- medrm(response ~ conc, fct=LN.4(),
                       data=nfa, random=d ~ 1 | experiment, start=start)
nfa.W14.mixed <- medrm(response ~ conc, fct=W1.4(),
                       data=nfa, random=d ~ 1 | experiment, start=start)
nfa.W24.mixed <- medrm(response ~ conc, fct=W2.4(),
                       data=nfa, random=d ~ 1 | experiment, start=start)

## information criteria
AIC(nfa.LL4.mixed, nfa.LN4.mixed,
    nfa.W14.mixed, nfa.W24.mixed)

## BMDL estimation for BMR 1, 5, 10, 20
mmaBMD(nfa.LL4.mixed, nfa.LN4.mixed,
      nfa.W14.mixed, nfa.W24.mixed,
      respLev=c(1,5,10,20),
      interval = "buckland", bmd = "additive")
}
}

\keyword{datasets}
