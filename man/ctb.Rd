\name{ctb}
\alias{ctb}
\docType{data}
\title{CellTiter-Blue Cell Viability Assay Data}
\description{Neurotoxicity test using the CellTiter-Blue Cell Viability
  Assay on SH-SY5Y cells for increasing concentrations of acrylamide.}
\usage{data(ctb)}
\format{
  A data frame with 647 observations on the following 5 variables.
  \describe{
    \item{\code{well}}{well ID of a 96 well plate}
    \item{\code{conc}}{12 concentrations of acrylamide, ranging from
      0-500mM}
    \item{\code{fluorescence}}{measured fluorescence after adding the
      resazurin reagent into the wells}
    \item{\code{day}}{integer denoting 3 different days}
    \item{\code{plate}}{factor with 7 levels representing the plate ID}
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
data(ctb)
ctb$day <- as.factor(ctb$day)
ctb$dayplate <- as.factor(with(ctb, paste(day, plate, sep="/")))

ggplot(ctb, aes(x=log(conc), y=fluorescence, colour=day, group=day:plate)) + 
  geom_point()

# starting values for fixed effects
fix <- coefficients(drm(fluorescence ~ conc, fct=LL.4(), data=ctb))

# starting values for random day effects
rday <- drm(fluorescence ~ conc, curveid=day, fct=LL.4(), data=ctb)
cmatday <- matrix(coefficients(rday), ncol=4)
mday <- apply(cmatday, 2, mean)
names(mday) <- letters[2:5]
rmatday <- t(apply(cmatday, 1, function(x) x-mday))
rownames(rmatday) <- levels(ctb$day)
colnames(rmatday) <- letters[2:5]

start <- list(fixed=fix, random=list(day=rmatday[,3, drop=FALSE]))


## set of nonlinear mixed models
ctb.LL4.mixed <- medrm(fluorescence  ~ conc, fct=LL.4(), 
                       data=ctb, random=d  ~ 1 | day/plate, start=start)
ctb.LN4.mixed <- medrm(fluorescence  ~ conc, fct=LN.4(), 
                       data=ctb, random=d  ~ 1 | day/plate, start=start)
ctb.W14.mixed <- medrm(fluorescence  ~ conc, fct=W1.4(), 
                       data=ctb, random=d  ~ 1 | day/plate, start=start)
ctb.W24.mixed <- medrm(fluorescence  ~ conc, fct=W2.4(), 
                       data=ctb, random=d  ~ 1 | day/plate, start=start)
ctb.FPL4b.mixed <- medrm(fluorescence ~ conc, fct=FPL.4(-1, 1), 
                         data=ctb, random=d ~ 1 | day/plate, start=start)
ctb.FPL4c.mixed <- medrm(fluorescence ~ conc, fct=FPL.4(-1, 2), 
                         data=ctb, random=d  ~ 1 | day/plate, start=start)
ctb.FPL4d.mixed <- medrm(fluorescence ~ conc, fct=FPL.4(-0.5, 3), 
                         data=ctb, random=d  ~ 1 | day/plate, start=start)

## information criteria
AIC(ctb.LL4.mixed, ctb.LN4.mixed, ctb.W14.mixed,
    ctb.W24.mixed, ctb.FPL4b.mixed,
    ctb.FPL4c.mixed, ctb.FPL4d.mixed)


}
}



\keyword{datasets}
