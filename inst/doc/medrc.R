## ----packages, echo=FALSE, warning=FALSE, message=FALSE------------------
library(medrc)
library(multcomp)

## ----vinclozolindata-----------------------------------------------------
data(vinclozolin)

## ----vinclozolinmixedmodel-----------------------------------------------
m <- medrm(effect ~ conc, data=vinclozolin,
           random=b + d + e ~ 1|exper,
           fct=LL.3(), start=c(0.5, 2000, 0.05))
print(m)

## ----vinclognls----------------------------------------------------------
mg <- glsdrm(effect ~ conc, data=vinclozolin,
             fct=LL.3(), start=c(0.5, 2000, 0.05), 
             correlation=corCompSymm(form=~1|exper), 
             control=gnlsControl(tolerance=0.01))
print(mg)

## ----vincloplot, warning=FALSE, message=FALSE, fig.width=6---------------
nd <- expand.grid(conc=exp(seq(log(0.01), log(3.2), length=25)))
nd1 <- nd2 <- nd3 <- nd
nd1$p <- predict(m, type="marginal", newdata=nd)
nd1$model <- "Marginalized"
nd2$p <- predict(mg, newdata=nd)
nd2$model <- "GNLS"
nd3$p <- m$fct$fct(nd$conc, rbind(coefficients(m)))
nd3$model <- "Assay-specific"
nd <- rbind(nd1, nd2, nd3)
nd4 <- expand.grid(conc=exp(seq(log(0.01), log(3.2), length=25)), 
                   exper=unique(vinclozolin$exper))
nd4$p <- predict(m, newdata=nd4)
vinclozolin2 <- vinclozolin
vinclozolin2$conc[vinclozolin2$conc == 0] <- 0.01

ggplot(vinclozolin2, aes(x=conc, y=effect, group=exper, colour=exper)) + 
  geom_point() + 
  geom_line(data=nd4, aes(y=p, colour=exper), linetype=3) +
  geom_line(data=nd, aes(y=p, linetype=model, group=NULL, colour=NULL), size=1.1) +
  coord_trans(x="log") + 
  theme_classic() +
  scale_linetype_manual("Model", values=c(2,1,4)) +
  ylab("Chemiluminescence") +
  xlab("Concentration [muM]") +
  scale_x_continuous(breaks=c(0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 1, 2, 3), 
                    labels=c(0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 1, 2, 3))

## ----vincloED------------------------------------------------------------
ED(m, respLev=c(10, 25, 50, 75, 90))

## ----vinclomarg----------------------------------------------------------
EDmarg(m, respLev=c(10, 25, 50, 75, 90))
ED(mg, respLev=c(10, 25, 50, 75, 90))

## ----spinachdata---------------------------------------------------------
data(spinach)
spinach$CURVE <- as.factor(spinach$CURVE)

## ----spinachmodel--------------------------------------------------------
### NLME
spm <- medrm(SLOPE ~ DOSE, 
           curveid=b + c + d + e ~ HERBICIDE, 
           data=spinach, 
           random=b + c + d + e ~ 1|CURVE,
           fct=LL.4())

print(spm)

### GNLS
spmg <- glsdrm(SLOPE ~ DOSE, 
           curveid=b + c + d + e ~ HERBICIDE, 
           data=spinach, 
           correlation=corCompSymm(form=~1|CURVE),
           fct=LL.4())

print(spmg)

## ----spinachplot, warning=FALSE, message=FALSE, fig.width=6--------------
snd <- expand.grid(DOSE=exp(seq(log(0.01), log(150), length=100)), 
                   HERBICIDE=levels(spinach$HERBICIDE))
snd1 <- snd2 <- snd
snd1$p <- predict(spm, type="marginal", newdata=snd)
snd1$model <- "Marginalized"
snd2$p <- predict(spmg, newdata=snd)
snd2$model <- "GNLS"
snd <- rbind(snd1, snd2)
snd3 <- expand.grid(DOSE=exp(seq(log(0.01), log(150), length=100)), 
                    HERBICIDE=levels(spinach$HERBICIDE), 
                    CURVE=as.factor(1:5))
snd3$p <- predict(spm, newdata=snd3)
spinach2 <- spinach
spinach2$DOSE[spinach2$DOSE == 0] <- 0.01

ggplot(spinach2, aes(x=DOSE, y=SLOPE, group=HERBICIDE, shape=HERBICIDE, colour=CURVE)) + 
  geom_point(size=2.5) + 
  geom_line(data=snd3, aes(y=p, group=CURVE:HERBICIDE), linetype=3) +
  geom_line(data=snd, aes(y=p, linetype=model, group=HERBICIDE:as.factor(snd$model), 
                          shape=NULL, colour=NULL), size=1.1) +
  coord_trans(x="log") + 
  theme_classic() +
  scale_linetype_manual("Model", values=c(2,1,4)) +
  scale_shape_discrete("Herbicide") +
  ylab("Oxygen consumption of thylakoid membranes") +
  xlab("Concentration [muM]") +
  scale_x_continuous(breaks=c(0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 1,2.5, 5,10,25,100), 
                     labels=c(0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 1,2.5, 5,10,25,100))

## ----spinachSI-----------------------------------------------------------
edm <- EDmarg(spm, respLev=c(10,25,50,75,90))[[2]]
library(mratios)
Cnum <- cbind(diag(5),matrix(0,5,5))
Cden <- cbind(matrix(0,5,5),diag(5))
rownames(Cnum) <- rownames(Cden) <- paste("bentazon/diuron ED:", c(10,25,50,75,90))
gsci.ratio(edm$coef, edm$vcov, Cnum, Cden, degfree=0)

## ----ctbdata-------------------------------------------------------------
data(ctb)
ctb$day <- as.factor(ctb$day)
ctb$dayplate <- as.factor(with(ctb, paste(day, plate, sep="/")))

## ----ctbmodel, warning=FALSE---------------------------------------------
## medrm fit
mv <- medrm(fluorescence  ~ conc, fct=LL.5(fixed=c(NA, 0, NA, NA, NA)), 
            data=ctb, random=b+d+e+f~1|day/plate, 
            weights=varExp(form = ~conc|plate), 
            control=nlmeControl(maxIter=200))

summary(mv)

## ----ctbplot, warning=FALSE, message=FALSE, fig.width=6------------------
cnd <- expand.grid(conc=exp(seq(log(0.001), log(500), length=100)))
cnd1 <- cnd
cnd1$p <- predict(mv$fit, level=0, newdata=cnd)
cnd2 <- expand.grid(conc=exp(seq(log(0.001), log(500), length=100)), 
                    day=as.factor(levels(ctb$day)))
cnd2$p <- as.vector(apply(coef(mv$fit, level=1), 1, 
                          function(x) mv$fct$fct(cnd$conc, rbind(x))))
cnd3 <- expand.grid(conc=exp(seq(log(0.001), log(500), length=100)), 
                    dayplate=levels(ctb$dayplate))
cnd3$p <- as.vector(apply(coef(mv$fit, level=2), 1, 
                          function(x) mv$fct$fct(cnd$conc, rbind(x))))
ctb2 <- ctb
ctb2$conc[ctb2$conc == 0] <- 0.001

ggplot(ctb2, aes(x=conc, y=fluorescence, group=dayplate, shape=day, colour=dayplate)) + 
  geom_point(size=1.5) + 
  geom_line(data=cnd3, aes(y=p, group=dayplate, shape=NULL), linetype=2) +
  geom_line(data=cnd2, aes(y=p, group=day, colour=NULL), size=1.1, linetype=2, colour="grey2") +
  geom_line(data=cnd1, aes(y=p, group=NULL, shape=NULL, colour=NULL), size=1.2) +
  coord_trans(x="log") + 
  theme_classic() +
  guides(shape="none") +
  ylab("Fluorescence") +
  xlab("Concentration [mM]") +
  scale_x_continuous(breaks=c(0.001, 0.01, 0.05, 0.25, 1, 5,25,100,500), 
                     labels=c(0.001, 0.01, 0.05, 0.25, 1, 5,25,100,500))

## ----ctbBMDL-------------------------------------------------------------
BMDmarg(mv, respLev=c(5, 10, 25, 50), nGQ=2, interval="tfls", level=0.9)

## ----mdradata------------------------------------------------------------
data(mdra)

## ----mdramodel, fig.width=6----------------------------------------------
mdramod <- medrm(Response ~ Concentration, data=mdra, fct=LL.3(), 
           random=d + e ~ 1|LabID/ExperimentID, 
           weights=varExp(form=~Concentration)) 
plot(mdramod, logx=TRUE, ndose=250, ranef=TRUE) + theme_classic()

## ----mdrabmd-------------------------------------------------------------
bmdra <- BMDmarg(mdramod, respLev=c(1, 5, 10), 
                 interval="tfls", rfinterval=c(0,10), nGQ=5)

## ----broccoli, fig.width=6-----------------------------------------------
data(broccoli)
str(broccoli)

ggplot(broccoli, aes(x=Day, y=LeafLength, group=ID, colour=Stress)) +
  geom_line() +
  facet_wrap(~ Genotype, ncol=8)

## ----removeobs-----------------------------------------------------------
bro <- droplevels(subset(broccoli, ID != "110" & ID != "125"))

## ----5pl-----------------------------------------------------------------
m5pl <- medrm(LeafLength ~ Day, data=bro, 
              fct=L.5(), 
              curveid=b + d + e + f ~ Stress, 
              random=d + e ~ 1|Genotype/ID)

print(m5pl)

## ----vcorr---------------------------------------------------------------
VarCorr(m5pl)
# same as
VarCorr(m5pl$fit)

## ----ranef, fig.width=6--------------------------------------------------
re <- ranef(m5pl)[[1]]
head(re)

panellab <- function(x, y, ...){
  abline(h=0, lty=2)
  abline(v=0, lty=2)
  text(x, y, rownames(re))
}
pairs(re, panel = panellab)

## ----plotfixed, fig.width=6----------------------------------------------
plot(m5pl) +
  geom_line(data=bro, aes(group=ID), linetype=2, alpha=0.2) +
  theme_bw()  

## ----resvsfitted---------------------------------------------------------
plot(residuals(m5pl) ~ fitted(m5pl))
abline(h=0, lty=2, col="red3")

## ----multimodel, eval=FALSE----------------------------------------------
#  # 3 parameter logistic with lower asymptote fixed at 0
#  mod1 <- medrm(LeafLength ~ Day, data=bro,
#                fct=L.3(),
#                curveid=b + d + e ~ Stress,
#                random=d + e ~ 1|Genotype/ID)
#  # 4 parameter logistic
#  # with the same lower asymptote for both stress treatments
#  mod2 <- medrm(LeafLength ~ Day, data=bro,
#                fct=L.4(fixed=c(NA, 5, NA, NA)),
#                curveid=b + d + e ~ Stress,
#                random=d + e ~ 1|Genotype/ID)
#  # 4 parameter Weibull model
#  mod3 <- medrm(LeafLength ~ Day, data=bro,
#                fct=W1.4(),
#                curveid=b + d + e ~ Stress,
#                random=d + e ~ 1|Genotype/ID)
#  # 2nd parameterization of 4 parameter Weibull model
#  mod4 <- medrm(LeafLength ~ Day, data=bro,
#                fct=W2.4(),
#                curveid=b + d + e ~ Stress,
#                random=d + e ~ 1|Genotype/ID)
#  # even a 4 parameter logistic with same parameters
#  # for both stress treatments is available with a onesided curveid formula
#  mod5 <- medrm(LeafLength ~ Day, data=bro,
#                fct=L.4(),
#                curveid= ~ Stress,
#                random=d + e ~ 1|Genotype/ID)
#  # or a 4p-log-logistic model with a lower asymptote fixed at 5
#  mod6 <- medrm(LeafLength ~ Day, data=bro,
#                fct=L.4(fixed=c(NA, 5, NA, NA)),
#                curveid=d + e ~ Stress,
#                random=d + e ~ 1|Genotype/ID)

## ----multimodelplot, eval=FALSE------------------------------------------
#  mmplot(mod1, mod2, mod3, mod4, mod5, mod6, ndose=50)

## ----muco, message=FALSE-------------------------------------------------
library(multcomp)
K <- rbind("drought-control | b"=c(-1, 1,  0, 0,  0, 0,  0, 0,  0, 0),
           "drought-control | d"=c( 0, 0,  0, 0, -1, 1,  0, 0,  0, 0),
           "drought-control | e"=c( 0, 0,  0, 0,  0, 0, -1, 1,  0, 0),
           "drought-control | f"=c( 0, 0,  0, 0,  0, 0,  0, 0, -1, 1))
gg <- glht(m5pl, linfct=K)
summary(gg)

## ----EDdrc---------------------------------------------------------------
ED(m5pl, respLev=c(25, 50, 75), interval="delta")

## ----SIdrc---------------------------------------------------------------
SI(m5pl, percVec=c(25, 25), interval="delta")

## ----mmaED, eval=FALSE---------------------------------------------------
#  mmaED(mod1, mod2, mod3, mod4, mod5, mod6, respLev=c(25, 50, 75), interval="kang")
#  

## ----EDSImarg------------------------------------------------------------
EDmarg(m5pl, respLev=c(25, 50, 75), interval="delta", nGQ=3)

