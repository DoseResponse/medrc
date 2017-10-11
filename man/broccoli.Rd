\name{broccoli}
\alias{broccoli}
\docType{data}
\title{The Effects of Drought Stress on Leaf Development in a \emph{Brassica oleracea} population}
\description{The effect of drought stress on \emph{Brassica oleracea} should be investigated, selecting drought stress resistant out of a population of different DH genotypes. The study was carried out on 48 DH lines developed from F1 plants of a cross between the rapid cycling chinese kale (\emph{Brassica oleracea} var. \emph{alboglabra} (L.H. Bailey) Musil) and broccoli (\emph{Brassica oleracea} var. \emph{italica} Plenck). 2 stress treatments (not watered and a watered control) are randomly assigned to 4 plants per genotype (2 per treatment) resulting in 192 plants in total. For the genotypes 5, 17, 31, 48, additional 12 plants (6 per treatment) are included into the completely randomized design, which results in a total of 240 plants. For each plant the length of the youngest leaf at the beginning of the experiment is measured daily for a period of 16 days. For the additional 12 plants of the 4 genotypes the leaf water potential was measured as a secondary endpoint (omitted here); due to these destructive measurements some dropouts occur.}
\usage{data(broccoli)}
\format{
  A data frame with 3689 observations on the following 5 variables.
  \describe{
    \item{\code{LeafLength}}{Length of the youngest leaf [cm]}
    \item{\code{ID}}{Plant identifier for 240 plants}
    \item{\code{Stress}}{Drought stress treatment with 2 levels (control/drought)}
    \item{\code{Genotype}}{Genotype ID with 48 levels}
    \item{\code{Day}}{Day of repeated measurement (1,2,...,16)}
  }
}
\references{
Uptmoor, R.; Osei-Kwarteng, M.; Guertler, S. & Stuetzel, H. Modeling the Effects of Drought Stress on Leaf Development in a Brassica oleracea Doubled Haploid Population Using Two-phase Linear Functions. Journal of the American Society for Horticultural Science, 2009, 134, 543-552.
}

\examples{
data(broccoli)

\dontrun{
# 5-parameter log-logistic model
# different steepness, asymptotes, and inflection point per stress treatment level
# random asymptotes and inflection point per genotype level and plant within genotypes
# no additional structure on random effects, estimating the pairwise random effect correlation 
#   at each hierarchical level

bm <- medrm(LeafLength ~ Day, curveid=b + c + e + d ~ Stress,
            random=c + e + d ~ 1|Genotype/ID,
            data=broccoli, fc=LL.5()) 

# estimated coefficients
print(bm)
summary(bm)

# show estimated variance components
VarCorr(bm)

# plot curves (based on fixed effects conditional on random effects of 0)
plot(bm)
plot(bm, logx=TRUE) # log tranformed scale on x-axis

# estimation of ED25, ED50, and ED75 conditional on random effects of 0
ED(bm, c(25, 50, 75), interval="tfls")

# testing difference of ED50 values
SI(bm, c(50, 50))

#################################################
# marginal model by generalized least squares
bmgls <- glsdrm(LeafLength ~ Day, curveid=c + d + e ~ Stress,
            correlation=corAR1(0.95, form= ~ Day|ID),
            data=broccoli, fc=LL.5())
  
# estimated coefficients
print(bmgls)
summary(bmgls)  
  
# plot marginal curves
plot(bmgls)
plot(bmgls, logx=TRUE) # log tranformed scale on x-axis            
}

}


\keyword{datasets}
