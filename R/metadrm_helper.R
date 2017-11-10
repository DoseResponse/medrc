### fixed/random/mixed-effects multivariate/multilevel model with:
###    - possibly one or multiple random intercepts (sigma2) with potentially known correlation matrices
###    - possibly correlated random effects for arms/groups/levels within studies (tau2 and rho for 1st term, gamma2 and phi for 2nd term)
### model also allows for correlated sampling errors via non-diagonal V matrix

# V      = variance-covariance matrix of the sampling errors
# sigma2 = (preset) value(s) for the variance of the random intercept(s)
# tau2   = (preset) value(s) for the variance of the random effects
# rho    = (preset) value(s) for the correlation(s) between random effects
# gamma2 = (preset) value(s) for the variance of the random effects
# phi    = (preset) value(s) for the correlation(s) between random effects

### structures when there is a (~ inner | outer) term in the random argument:
### - CS   (compound symmetry)
### - HCS  (heteroscedastic compound symmetry)
### - UN   (general positive-definite matrix with no structure)
### - UNHO (homoscedastic general positive-definite matrix with single tau2/gamma2 value and unstructured correlation matrix)
### - AR   (AR1 structure with a single tau2/gamma2 value and autocorrelation rho/phi)
### - HAR  (heteroscedastic AR1 structure with multiple tau2/gamma2 values and autocorrelation rho/phi)
### - ID   (same as CS but with rho/phi=0)
### - DIAG (same as HCS but with rho/phi=0)

rmadrc <- function(yi, V, W, mods, random, struct="CS", intercept=TRUE, data, slab, subset, ### add ni as argument in the future
method="REML", test="z", level=95, digits=4, btt, R, Rscale="cor", sigma2, tau2, rho, gamma2, phi, sparse=FALSE, verbose=FALSE, control, ...) {

   #########################################################################

   ###### data setup

   ### check argument specifications

   if (!is.element(method, c("FE","ML","REML")))
      stop("Unknown 'method' specified.")

   if (any(!is.element(struct, c("CS","HCS","UN","AR","HAR","UNHO","ID","DIAG"))))
      stop("Unknown 'struct' specified.")

   if (length(struct) == 1)
      struct <- c(struct, struct)

   na.act <- getOption("na.action")

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop("Unknown 'na.action' specified under options().")

   if (missing(random))
      random <- NULL

   if (missing(R))
      R <- NULL

   if (missing(sigma2))
      sigma2 <- NULL

   if (missing(tau2))
      tau2 <- NULL

   if (missing(rho))
      rho <- NULL

   if (missing(gamma2))
      gamma2 <- NULL

   if (missing(phi))
      phi <- NULL

   if (missing(control))
      control <- list()

   ### get ... argument and check for extra/superfluous arguments

   ddd <- list(...)

   .chkdots(ddd, c("tdist", "outlist"))

   ### handle 'tdist' argument from ...

   if (is.logical(ddd$tdist) && !ddd$tdist)
      test <- "z"
   if (is.logical(ddd$tdist) && ddd$tdist)
      test <- "t"

   if (!is.element(test, c("z","t","knha","adhoc")))
      stop("Invalid option selected for 'test' argument.")

   ### deal with Rscale argument (either character, logical, or integer)

   if (is.character(Rscale))
      Rscale <- match.arg(Rscale, c("none", "cor", "cor0", "cov0"))

   if (is.logical(Rscale))
      Rscale <- ifelse(Rscale, "cor", "none")

   if (is.numeric(Rscale)) {
      Rscale <- round(Rscale)
      if (Rscale > 3 | Rscale < 0)
         stop("Unknown 'Rscale' value specified.")
      Rscale <- switch(as.character(Rscale), "0"="none", "1"="cor", "2"="cor0", "3"="cov0")
   }

   #########################################################################

   if (verbose > 1)
      message("Extracting yi/V values ...")

   ### check if data argument has been specified

   if (missing(data))
      data <- NULL

   if (is.null(data)) {
      data <- sys.frame(sys.parent())
   } else {
      if (!is.data.frame(data))
         data <- data.frame(data)
   }

   mf <- match.call()

   ### extract yi, V, W, ni, slab, subset, and mods values, possibly from the data frame specified via data (arguments not specified are NULL)

   mf.yi     <- mf[[match("yi", names(mf))]]
   mf.V      <- mf[[match("V",  names(mf))]]
   mf.W      <- mf[[match("W",  names(mf))]]
   mf.ni     <- mf[[match("ni", names(mf))]] ### not yet possible to specify this
   mf.slab   <- mf[[match("slab",   names(mf))]]
   mf.subset <- mf[[match("subset", names(mf))]]
   mf.mods   <- mf[[match("mods",   names(mf))]]
   yi     <- eval(mf.yi,     data, enclos=sys.frame(sys.parent())) ### NULL if user does not specify this
   V      <- eval(mf.V,      data, enclos=sys.frame(sys.parent())) ### NULL if user does not specify this
   W      <- eval(mf.W,      data, enclos=sys.frame(sys.parent())) ### NULL if user does not specify this
   ni     <- eval(mf.ni,     data, enclos=sys.frame(sys.parent())) ### NULL if user does not specify this
   slab   <- eval(mf.slab,   data, enclos=sys.frame(sys.parent())) ### NULL if user does not specify this
   subset <- eval(mf.subset, data, enclos=sys.frame(sys.parent())) ### NULL if user does not specify this
   mods   <- eval(mf.mods,   data, enclos=sys.frame(sys.parent())) ### NULL if user does not specify this

   ### if yi is a formula, extract yi and X (this overrides anything specified via the mods argument further below)

   #is.formula <- FALSE

   if (inherits(yi, "formula")) {
      options(na.action = "na.pass")                   ### set na.action to na.pass, so that NAs are not filtered out (we'll do that later)
      mods <- model.matrix(yi, data=data)              ### extract model matrix (now mods is no longer a formula, so [a] further below is skipped)
      attr(mods, "assign") <- NULL                     ### strip assign attribute (not needed at the moment)
      yi <- model.response(model.frame(yi, data=data)) ### extract yi values from model frame
      options(na.action = na.act)                      ### set na.action back to na.act
      names(yi) <- NULL                                ### strip names (1:k) from yi (so res$yi is the same whether yi is a formula or not)
      intercept <- FALSE                               ### set to FALSE since formula now controls whether the intercept is included or not
      #is.formula <- TRUE                               ### note: code further below ([b]) actually checks whether intercept is included or not
   }

   ### in case user passed a matrix to yi, convert it to a vector

   if (is.matrix(yi))
      yi <- as.vector(yi)

   ### number of outcomes before subsetting

   k <- length(yi)
   k.all <- k

   ### set default measure argument

   measure <- "GEN"

   if (!is.null(attr(yi, "measure"))) ### take 'measure' from yi (if it is there)
      measure <- attr(yi, "measure")

   ### add measure attribute (back) to the yi vector

   attr(yi, "measure") <- measure

   ### some checks on V (and turn V into a diagonal matrix if it is a column/row vector)

   if (is.null(V))
      stop("Need to specify 'V' argument.")

   if (is.list(V)) {

      if (any(!sapply(V, .is.square)))
         stop("All list elements in 'V' must be square matrices.")

      ### need to do this first, since is.vector(V) is TRUE for lists and so that further code works

      if (sparse) {
         V <- bdiag(V)
      } else {
         V <- bldiag(V)
      }

   }

   ### check if user constrained V to 0

   if (is.vector(V) && length(V) == 1 && V == 0) {
      V0 <- TRUE
   } else {
      V0 <- FALSE
   }

   ### turn V into a diagonal matrix if it is a column/row vector
   ### note: if V is a scalar (e.g., V=0), then this will turn V into a kxk
   ### matrix with the value of V along the diagonal

   if (is.vector(V) || nrow(V) == 1L || ncol(V) == 1L)
      V <- diag(as.vector(V), nrow=k, ncol=k)

   if (is.data.frame(V))
      V <- as.matrix(V)

   ### remove row and column names (important for isSymmetric() function)
   ### (but only do this if V has row/column names to avoid making an unnecessary copy)

   if (!is.null(dimnames(V)))
      V <- unname(V)

   ### check whether V is square and symmetric

   if (!.is.square(V))
      stop("'V' must be a square matrix.")

   if (!isSymmetric(V)) ### note: copy of V is made when doing this
      stop("'V' must be a symmetric matrix.")

   ### check length of yi and V

   if (nrow(V) != k)
      stop("Length of 'yi' and length/dimensions of 'V' are not the same.")

   ### force V to be sparse when sparse=TRUE (and V is not yet sparse)

   if (sparse && inherits(V, "matrix"))
      V <- Matrix(V, sparse=TRUE)

   ### process W if it was specified

   if (!is.null(W)) {

      ### turn W into a diagonal matrix if it is a column/row vector
      ### in general, turn W into A (arbitrary weight matrix)

      if (is.vector(W) || nrow(W) == 1L || ncol(W) == 1L) {

         W <- as.vector(W)

         ### allow easy setting of W to a single value

         if (length(W) == 1L)
            W <- rep(W, k)

         A <- diag(W, nrow=length(W), ncol=length(W))

      } else {

         A <- W

      }

      if (is.data.frame(A))
         A <- as.matrix(A)

      ### remove row and column names (important for isSymmetric() function)
      ### (but only do this if A has row/column names to avoid making an unnecessary copy)

      if (!is.null(dimnames(A)))
         A <- unname(A)

      ### check whether A is square and symmetric

      if (!.is.square(A))
         stop("'W' must be a square matrix.")

      if (!isSymmetric(A))
         stop("'W' must be a symmetric matrix.")

      ### check length of yi and A

      if (nrow(A) != k)
         stop("Length of 'yi' and length/dimensions of 'W' are not the same.")

      ### force A to be sparse when sparse=TRUE (and A is not yet sparse)

      if (sparse && inherits(A, "matrix"))
         A <- Matrix(A, sparse=TRUE)

   } else {

      A <- NULL

   }

   ### if ni has not been specified (and hence is NULL) but is an attribute of yi, get it
   ### note: currently ni argument removed, so this is the only way to pass ni to the function

   if (is.null(ni) && !is.null(attr(yi, "ni")))
      ni <- attr(yi, "ni")

   ### check length of yi and ni
   ### if there is a mismatch, then ni cannot be trusted, so set it to NULL

   if (!is.null(ni) && length(ni) != k)
      ni <- NULL

   ### if ni is now available, add it (back) as an attribute to yi
   ### this is currently pointless, but may be useful if function has an ni argument

   #if (!is.null(ni))
   #   attr(yi, "ni") <- ni

   #########################################################################

   if (verbose > 1)
      message("Creating model matrix ...")

   ### convert mods formula to X matrix and set intercept equal to FALSE
   ### skipped if formula has already been specified via yi argument, since mods is then no longer a formula (see [a])

   if (inherits(mods, "formula")) {
      options(na.action = "na.pass")        ### set na.action to na.pass, so that NAs are not filtered out (we'll do that later)
      mods <- model.matrix(mods, data=data) ### extract model matrix
      attr(mods, "assign") <- NULL          ### strip assign attribute (not needed at the moment)
      options(na.action = na.act)           ### set na.action back to na.act
      intercept <- FALSE                    ### set to FALSE since formula now controls whether the intercept is included or not
      #is.formula <- TRUE                    ### note: code further below ([b]) actually checks whether intercept is included or not
   }

   ### turn a row vector for mods into a column vector

   if (is.vector(mods))
      mods <- cbind(mods)

   ### turn a mods data frame into a matrix

   if (is.data.frame(mods))
      mods <- as.matrix(mods)

   ### check if model matrix contains character variables

   if (is.character(mods))
      stop("Model matrix contains character variables.")

   ### check if mods matrix has the right number of rows

   if (!is.null(mods) && (nrow(mods) != k))
      stop("Number of rows of the model matrix does not match length of the outcome vector.")

   #########################################################################
   #########################################################################
   #########################################################################

   ### process random argument

   if (method != "FE" && !is.null(random)) {

      if (verbose > 1)
         message("Processing 'random' argument ...")

      ### make sure random argument is always a list (so lapply() below works)

      if (!is.list(random))
         random <- list(random)

      ### figure out if a formulas has a slash (as in ~ 1 | study/id)

      has.slash <- sapply(random, function(f) grepl("/", paste0(f, collapse="")))

      ### substitute + for | in all formulas (so that model.frame() below works)

      random.plus <- lapply(random, function(f) formula(sub("\\|", "+", paste0(f, collapse=""))))

      ### get all model frames corresponding to the formulas in the random argument
      ### mf.r <- lapply(random, get_all_vars, data=data)
      ### note: get_all_vars() does not carry out any functions calls within the formula
      ###       so use model.frame(), which allows for things like 'random = ~ factor(arm) | study'
      ###       need to use na.pass so that NAs are passed through and not omitted (check for NAs is done below)

      mf.r <- lapply(random.plus, model.frame, data=data, na.action=na.pass)

      ### count number of columns in each model frame

      mf.r.ncols <- sapply(mf.r, ncol)

      ### for formulas with slashes, create interaction terms

      for (j in seq_along(has.slash)) {

         if (!has.slash[j])
            next

         ### need to go backwards; otherwise, with 3 or more terms (e.g., ~ 1 | var1/var2/var3), the third term would be an
         ### interaction between var1, var1:var2, and var3; by going backwards, we get var1, var1:var2, and var1:var2:var3

         for (p in mf.r.ncols[j]:1) {
            mf.r[[j]][,p] <- interaction(mf.r[[j]][1:p], drop=TRUE, lex.order=TRUE, sep = "/")
            colnames(mf.r[[j]])[p] <- paste(colnames(mf.r[[j]])[1:p], collapse="/")
         }

      }

      ### create list where model frames with multiple columns based on slashes are flattened out

      if (any(has.slash)) {

         if (length(mf.r) == 1) {

            ### if formula only has one element of the form ~ 1 | var1/var2/..., create a list of the data frames (each with one column)

            mf.r <- lapply(seq(ncol(mf.r[[1]])), function(x) mf.r[[1]][x])

         } else {

            ### if there are non-slash elements, then this flattens things out (obviously ...)

            mf.r <- unlist(mapply(function(mf, sl) if (sl) lapply(seq(mf), function(x) mf[x]) else list(mf), mf.r, has.slash), recursive=FALSE, use.names=FALSE)

         }

         ### recount number of columns in each model frame

         mf.r.ncols <- sapply(mf.r, ncol)

      }

      #return(mf.r)

      ### make sure that each model frame has no more than 2 columns

      if (any(mf.r.ncols > 2))
         stop("No more than two elements allowed in each formula of the 'random' argument.")

      ### make sure that there are only up to 2 model frames with 2 columns

      if (sum(mf.r.ncols == 2) > 2)
         stop("Only up to two '~ inner | outer' formulas allowed in the 'random' argument.")

      ### separate mf.r into mf.s (~ 1 | id), mf.g (~ inner | outer), and mf.h (~ inner | outer) parts

      mf.s <- mf.r[which(mf.r.ncols == 1)]      ### if there is no (~ 1 | factor) term, this is list() ([] so that we get a list of data frames)
      mf.g <- mf.r[[which(mf.r.ncols == 2)[1]]] ### if there is no 1st (~ inner | outer) terms, this is NULL ([[]] so that we get a data frame, not a list)
      mf.h <- mf.r[[which(mf.r.ncols == 2)[2]]] ### if there is no 2nd (~ inner | outer) terms, this is NULL ([[]] so that we get a data frame, not a list)

      ### if there is no (~ 1 | factor) term, then mf.s is list(), so turn that into NULL

      if (length(mf.s) == 0)
         mf.s <- NULL

      ### does the random argument include at least one (~ 1 | id) term?

      withS <- !is.null(mf.s)

      ### does the random argument include (~ inner | outer) terms?

      withG <- !is.null(mf.g)
      withH <- !is.null(mf.h)

      ### count number of rows in each model frame

      mf.r.nrows <- sapply(mf.r, nrow)

      ### make sure that rows in each model frame match the length of the data

      if (any(mf.r.nrows != k))
         stop("Length of variables specified via the 'random' argument does not match length of the data.")

      ### need this for profile(); with things like 'random = ~ factor(arm) | study', 'mf.r' contains variables 'factor(arm)' and 'study'
      ### but the former won't work when using the same formula for the refitting (same when using interaction() in the random formula)
      ### careful: with ~ 1 | interaction(var1, var2), mf.r will have 2 columns, but is actually a 'one variable' term
      ###      and with ~ interaction(var1, var2) | var3, mf.r will have 3 columns, but is actually a 'two variable' term
      ### mf.r.ncols above is correct even in these cases (since it is based on the model.frame() results), but need
      ### to be careful that this doesn't screw up anything in other functions

      mf.r <- lapply(random.plus, get_all_vars, data=data)

   } else {

      ### set defaults for some elements when method="FE"

      mf.r  <- NULL
      mf.s  <- NULL
      mf.g  <- NULL
      mf.h  <- NULL
      withS <- FALSE
      withG <- FALSE
      withH <- FALSE

   }

   #return(list(mf.r, mf.s, mf.g, mf.h))

   ### note: checks on NAs in mf.s, mf.g, and mf.h after subsetting (since NAs may be removed by subsetting)

   #########################################################################
   #########################################################################
   #########################################################################

   ### generate study labels if none are specified (or none can be found in yi argument)

   if (verbose > 1)
      message("Generating/extracting study labels ...")

   ### study ids (1:k sequence before subsetting)

   ids <- seq_len(k)

   ### if slab has not been specified but is an attribute of yi, get it

   if (is.null(slab)) {

      if (!is.null(attr(yi, "slab")))
         slab <- attr(yi, "slab")

      ### check length of yi and slab (only if slab is not NULL)
      ### if there is a mismatch, then slab cannot be trusted, so set it to NULL

      if (is.null(slab) && length(slab) != k)
         slab <- NULL

   }

   if (is.null(slab)) {

      slab.null <- TRUE
      slab      <- ids

   } else {

      if (anyNA(slab))
         stop("NAs in study labels.")

      if (length(slab) != k)
         stop("Study labels not of same length as data.")

      slab.null <- FALSE

   }

   ### if a subset of studies is specified

   if (!is.null(subset)) {

      if (verbose > 1)
         message("Subsetting ...")

      yi   <- yi[subset]
      V    <- V[subset,subset,drop=FALSE]
      A    <- A[subset,subset,drop=FALSE]
      ni   <- ni[subset]
      mods <- mods[subset,,drop=FALSE]
      slab <- slab[subset]
      mf.r <- lapply(mf.r, function(x) x[subset,,drop=FALSE])
      mf.s <- lapply(mf.s, function(x) x[subset,,drop=FALSE])
      mf.g <- mf.g[subset,,drop=FALSE]
      mf.h <- mf.h[subset,,drop=FALSE]
      ids  <- ids[subset]
      k    <- length(yi)

      attr(yi, "measure") <- measure ### add measure attribute back
      attr(yi, "ni")      <- ni      ### add ni attribute back

   }

   ### check if study labels are unique; if not, make them unique

   if (anyDuplicated(slab))
      slab <- .make.unique(slab)

   ### add slab attribute back

   attr(yi, "slab") <- slab

   ### get the sampling variances from the diagonal of V

   vi <- diag(V)

   ### check for non-positive sampling variances (and set negative values to 0)

   if (any(vi <= 0, na.rm=TRUE)) {
      allvipos <- FALSE
      if (!V0)
         warning("There are outcomes with non-positive sampling variances.")
      vi.neg <- vi < 0
      if (any(vi.neg, na.rm=TRUE)) {
         V[vi.neg,] <- 0 ### note: entire row set to 0 (so covariances are also 0)
         V[,vi.neg] <- 0 ### note: entire col set to 0 (so covariances are also 0)
         vi[vi.neg] <- 0
         warning("Negative sampling variances constrained to zero.")
      }
   } else {
      allvipos <- TRUE
   }

   ### save full data (including potential NAs in yi/vi/V/W/ni/mods)

   yi.f    <- yi
   vi.f    <- vi
   V.f     <- V
   W.f     <- A
   ni.f    <- ni
   mods.f  <- mods
   mf.r.f  <- mf.r ### needed for cooks.distance()
   #mf.g.f <- mf.g ### copied further below
   #mf.h.f <- mf.h ### copied further below
   #mf.s.f <- mf.s ### copied further below

   k.f <- k ### total number of observed outcomes including all NAs

   #########################################################################
   #########################################################################
   #########################################################################

   ### stuff that need to be done after subsetting

   if (withS) {

      if (verbose > 1)
         message(paste0("Processing '", paste0("~ 1 | ", sapply(mf.s, names), collapse=", "), "' term(s) ..."))

      ### get variables names in mf.s

      s.names <- sapply(mf.s, names) ### one name per term

      ### turn each variable in mf.s into a factor (and turn each column vector into just a vector)
      ### if a variable was a factor to begin with, this drops any unused levels, but order of existing levels is preserved

      mf.s <- lapply(mf.s, function(x) factor(x[[1]]))

      ### check if there are any NAs anywhere in mf.s

      if (any(sapply(mf.s, anyNA)))
         stop("No NAs allowed in variables specified in the 'random' argument.")

      ### how many (~ 1 | id) terms does the random argument include? (0 if none, but if withS is TRUE, must be at least 1)

      sigma2s <- length(mf.s)

      ### set default value(s) for sigma2 argument if it is unspecified

      if (is.null(sigma2))
         sigma2 <- rep(NA_real_, sigma2s)

      ### allow quickly setting all sigma2 values to a fixed value

      if (length(sigma2) == 1)
         sigma2 <- rep(sigma2, sigma2s)

      ### check if sigma2 is of the correct length

      if (length(sigma2) != sigma2s)
         stop(paste("Length of 'sigma2' argument (", length(sigma2), ") does not match actual number of variance components (", sigma2s, ").", sep=""))

      ### checks on any fixed values of sigma2 argument

      if (any(sigma2 < 0, na.rm=TRUE))
         stop("Specified value(s) of 'sigma2' must be non-negative.")

      ### get number of levels of each variable in mf.s (vector with one value per term)

      s.nlevels <- sapply(mf.s, nlevels)

      ### get levels of each variable in mf.s (list with levels for each variable)

      s.levels <- lapply(mf.s, levels)

      ### checks on R (note: do this after subsetting, so user can filter out ids with no info in R)

      if (is.null(R)) {

         withR <- FALSE
         Rfix  <- rep(FALSE, sigma2s)

      } else {

         if (verbose > 1)
            message("Processing 'R' argument ...")

         withR <- TRUE

         ### make sure R is always a list (so lapply() below works)

         if (is.data.frame(R) || !is.list(R))
            R <- list(R)

         ### check if R list has no names at all or some names are missing
         ### (if only some elements of R have names, then names(R) is "" for the unnamed elements, so use nchar()==0 to check for that)

         if (is.null(names(R)) || any(nchar(names(R)) == 0))
            stop("Argument 'R' must be a *named* list.")

         ### remove elements in R that are NULL (not sure why this is needed; why would anybody ever do this?)
         ### maybe this had something to do with functions that repeatedly call rma.mv(); so leave this be for now

         R <- R[!sapply(R, is.null)]

         ### turn all elements in R into matrices (this would fail with a NULL element)

         R <- lapply(R, as.matrix)

         ### match up R matrices based on the s.names (and correct names of R)
         ### so if a particular ~ 1 | id term has a matching id=R element, the corresponding R element is that R matrix
         ###    if a particular ~ 1 | id term does not have a matching id=R element, the corresponding R element is NULL

         R <- R[s.names]

         ### NULL elements in R would have no name, so this makes sure that all R elements have the correct s.names

         names(R) <- s.names

         ### check for which components an R matrix has been specified

         Rfix <- !sapply(R, is.null)

         ### Rfix could be all FALSE (if user has used id names in R that are not actually in 'random')
         ### so only do the rest below if that is *not* the case

         if (any(Rfix)) {

            ### check if given R matrices are square and symmetric

            if (any(!sapply(R[Rfix], .is.square)))
               stop("Elements of 'R' must be square matrices.")
            if (any(!sapply(R[Rfix], function(x) isSymmetric(unname(x)))))
               stop("Elements of 'R' must be symmetric matrices.")

            for (j in seq_along(R)) {

               if (!Rfix[j])
                  next

               ### even if isSymmetric() is TRUE, there may still be minor numerical differences between the lower and upper triangular
               ### parts that could lead to isSymmetric() being FALSE once we do any potentially rescaling of the R matrices further
               ### below; this ensures strict symmetry to avoid this issue
               #R[[j]][lower.tri(R[[j]])] <- t(R[[j]])[lower.tri(R[[j]])]
               R[[j]] <- symmpart(R[[j]])

               ### if rownames are missing, copy colnames to rownames and vice-versa

               if (is.null(rownames(R[[j]])))
                  rownames(R[[j]]) <- colnames(R[[j]])
               if (is.null(colnames(R[[j]])))
                  colnames(R[[j]]) <- rownames(R[[j]])

               ### if colnames are still missing at this point, R element did not have dimension names to begin with

               if (is.null(colnames(R[[j]])))
                  stop("Elements of 'R' must have dimension names.")

            }

            ### if user specifies the entire (k x k) correlation matrix, this removes the duplicate rows/columns

            #R[Rfix] <- lapply(R[Rfix], unique, MARGIN=1)
            #R[Rfix] <- lapply(R[Rfix], unique, MARGIN=2)

            ### no, the user can specify an entire (k x k) matrix; the problem is repeated dimension names
            ### so let's filter out rows/columns with the same dimension names

            R[Rfix] <- lapply(R[Rfix], function(x) x[!duplicated(rownames(x)), !duplicated(colnames(x)), drop=FALSE])

            ### after the two commands above, this should always be FALSE, but leave for now just in case

            if (any(sapply(R[Rfix], function(x) length(colnames(x)) != length(unique(colnames(x))))))
               stop("Each element of 'R' must have unique dimension names.")

            ### check for R being positive definite
            ### skipped: even if R is not positive definite, the marginal var-cov matrix can still be; so just check for pd during optimization

            #if (any(sapply(R[Rfix], function(x) any(eigen(x, symmetric=TRUE, only.values=TRUE)$values <= .Machine$double.eps)))) ### any eigenvalue below double.eps is essentially 0
            #   stop("Matrix in R is not positive definite.")

            for (j in seq_along(R)) {

               if (!Rfix[j])
                  next

               ### check if there are NAs in a matrix specified via R

               if (anyNA(R[[j]]))
                  stop("No missing values allowed in matrices specified via 'R'.")

               ### check if there are levels in s.levels which are not in R (if yes, issue an error and stop)

               if (any(!is.element(s.levels[[j]], colnames(R[[j]]))))
                  stop(paste0("There are levels in '", s.names[j], "' for which there are no rows/columns in the corresponding 'R' matrix."))

               ### check if there are levels in R which are not in s.levels (if yes, issue a warning)

               if (any(!is.element(colnames(R[[j]]), s.levels[[j]])))
                  warning(paste0("There are rows/columns in the 'R' matrix for '", s.names[j], "' for which there are no data."))

            }

         } else {

            warning("Argument 'R' specified, but list name(s) not in 'random'.")

            withR <- FALSE
            Rfix  <- rep(FALSE, sigma2s)
            R     <- NULL

         }

      }

   } else {

      ### need one fixed sigma2 value for optimization function

      sigma2s <- 1
      sigma2  <- 0

      s.nlevels <- NULL
      s.levels  <- NULL
      s.names   <- NULL

      withR   <- FALSE
      Rfix    <- FALSE
      R       <- NULL

   }

   mf.s.f <- mf.s

   ### copy s.nlevels and s.levels

   s.nlevels.f <- s.nlevels
   s.levels.f  <- s.levels

   #########################################################################

   ### stuff that need to be done after subsetting

   if (withG) {

      tmp <- .process.G.aftersub(verbose, mf.g, struct[1], tau2, rho, isG=TRUE, k, sparse)

      mf.g      <- tmp$mf.g
      g.names   <- tmp$g.names
      g.nlevels <- tmp$g.nlevels
      g.levels  <- tmp$g.levels
      tau2s     <- tmp$tau2s
      rhos      <- tmp$rhos
      tau2      <- tmp$tau2
      rho       <- tmp$rho
      Z.G1      <- tmp$Z.G1
      Z.G2      <- tmp$Z.G2

   } else {

      ### need one fixed tau2 and rho value for optimization function

      tau2s <- 1
      rhos  <- 1
      tau2  <- 0
      rho   <- 0

      ### need Z.G1 and Z.G2 to exist further below and for optimization function

      Z.G1 <- NULL
      Z.G2 <- NULL
      g.nlevels <- NULL
      g.levels  <- NULL

      g.names <- NULL

   }

   mf.g.f <- mf.g ### needed for predict()

   #########################################################################

   ### stuff that need to be done after subsetting

   if (withH) {

      tmp <- .process.G.aftersub(verbose, mf.h, struct[2], gamma2, phi, isG=FALSE, k, sparse)

      mf.h      <- tmp$mf.g
      h.names   <- tmp$g.names
      h.nlevels <- tmp$g.nlevels
      h.levels  <- tmp$g.levels
      gamma2s   <- tmp$tau2s
      phis      <- tmp$rhos
      gamma2    <- tmp$tau2
      phi       <- tmp$rho
      Z.H1      <- tmp$Z.G1
      Z.H2      <- tmp$Z.G2

   } else {

      ### need one fixed gamma2 and phi value for optimization function

      gamma2s <- 1
      phis    <- 1
      gamma2  <- 0
      phi     <- 0

      ### need Z.H1 and Z.H2 to exist further below and for optimization function

      Z.H1 <- NULL
      Z.H2 <- NULL
      h.nlevels <- NULL
      h.levels  <- NULL

      h.names <- NULL

   }

   mf.h.f <- mf.h ### needed for predict()

   #########################################################################
   #########################################################################
   #########################################################################

   ### check for NAs and act accordingly

   has.na <- is.na(yi) | (if (is.null(mods)) FALSE else apply(is.na(mods), 1, any)) | .anyNAv(V) | (if (is.null(A)) FALSE else apply(is.na(A), 1, any))
   not.na <- !has.na

   if (any(has.na)) {

      if (verbose > 1)
         message("Handling NAs ...")

      if (na.act == "na.omit" || na.act == "na.exclude" || na.act == "na.pass") {

         yi   <- yi[not.na]
         V    <- V[not.na,not.na,drop=FALSE]
         A    <- A[not.na,not.na,drop=FALSE]
         vi   <- vi[not.na]
         ni   <- ni[not.na]
         mods <- mods[not.na,,drop=FALSE]
         mf.r <- lapply(mf.r, function(x) x[not.na,,drop=FALSE])
         mf.s <- lapply(mf.s, function(x) x[not.na]) ### note: mf.s is a list of vectors at this point
         mf.g <- mf.g[not.na,,drop=FALSE]
         mf.h <- mf.h[not.na,,drop=FALSE]
         Z.G1 <- Z.G1[not.na,,drop=FALSE]
         Z.G2 <- Z.G2[not.na,,drop=FALSE]
         Z.H1 <- Z.H1[not.na,,drop=FALSE]
         Z.H2 <- Z.H2[not.na,,drop=FALSE]
         k    <- length(yi)
         warning("Rows with NAs omitted from model fitting.")

         attr(yi, "measure") <- measure ### add measure attribute back
         attr(yi, "ni")      <- ni      ### add ni attribute back

         ### note: slab is always of the same length as the full yi vector (after subsetting), so missings are not removed and slab is not added back to yi

      }

      if (na.act == "na.fail")
         stop("Missing values in data.")

   }

   ### more than one study left?

   if (k <= 1)
      stop("Processing terminated since k <= 1.")

   ### check for V being positive definite (this should also cover non-positive variances)
   ### skipped: even if V is not positive definite, the marginal var-cov matrix can still be; so just check for pd during the optimization
   ### but at least issue a warning, since a fixed-effects model can then not be fitted and there is otherwise no indication why

   if (!V0 && any(eigen(V, symmetric=TRUE, only.values=TRUE)$values <= .Machine$double.eps)) ### any eigenvalue below double.eps is essentially 0
      warning("'V' appears to be not positive definite.")

   ### check ratio of largest to smallest sampling variance
   ### note: need to exclude some special cases (0/0 = Nan, max(vi)/0 = Inf)

 #  vimaxmin <- max(vi) / min(vi)

#   if (!is.nan(vimaxmin) && !is.infinite(vimaxmin) && vimaxmin >= 1e+07)
 #     warning("Ratio of largest to smallest sampling variance extremely large. Cannot obtain stable results.")

   ### make sure that there is at least one column in X ([b])

   if (is.null(mods) && !intercept) {
      warning("Must either include an intercept and/or moderators in model.\n  Coerced intercept into the model.")
      intercept <- TRUE
   }

   ### add vector of 1s to the X matrix for the intercept (if intercept=TRUE)

   if (intercept) {
      X   <- cbind(intrcpt=rep(1,k), mods)
      X.f <- cbind(intrcpt=rep(1,k.f), mods.f)
   } else {
      X   <- mods
      X.f <- mods.f
   }

   ### drop redundant predictors
   ### note: need to save coef.na for functions that modify the data/model and then refit the model (regtest() and the
   ### various function that leave out an observation); so we can check if there are redundant/dropped predictors then

   tmp <- lm(yi ~ X - 1)
   coef.na <- is.na(coef(tmp))
   if (any(coef.na)) {
      warning("Redundant predictors dropped from the model.")
      X   <- X[,!coef.na,drop=FALSE]
      X.f <- X.f[,!coef.na,drop=FALSE]
   }

   ### check whether intercept is included and if yes, move it to the first column (NAs already removed, so na.rm=TRUE for any() not necessary)

   is.int <- apply(X, 2, .is.intercept)
   if (any(is.int)) {
      int.incl <- TRUE
      int.indx <- which(is.int, arr.ind=TRUE)
      X        <- cbind(intrcpt=1,   X[,-int.indx, drop=FALSE]) ### this removes any duplicate intercepts
      X.f      <- cbind(intrcpt=1, X.f[,-int.indx, drop=FALSE]) ### this removes any duplicate intercepts
      #if (is.formula)
         intercept <- TRUE ### set intercept appropriately so that the predict() function works
   } else {
      int.incl <- FALSE
   }

   p <- NCOL(X) ### number of columns in X (including the intercept if it is included)

   ### check whether this is an intercept-only model

   if ((p == 1L) && .is.intercept(X)) {
      int.only <- TRUE
   } else {
      int.only <- FALSE
   }

   ### check if there are too many parameters for given k (currently skipped)

   ### set/check 'btt' argument

   btt <- .set.btt(btt, p, int.incl)
   m <- length(btt) ### number of betas to test (m = p if all betas are tested)

   #########################################################################
   #########################################################################
   #########################################################################

   ### stuff that need to be done after subsetting and filtering out NAs

   if (withS) {

      ### redo: turn each variable in mf.s into a factor (reevaluates the levels present, but order of existing levels is preserved)

      mf.s <- lapply(mf.s, factor)

      ### redo: get number of levels of each variable in mf.s (vector with one value per term)

      s.nlevels <- sapply(mf.s, nlevels)

      ### redo: get levels of each variable in mf.s

      s.levels <- lapply(mf.s, levels)

      ### for any single-level factor with unfixed sigma2, fix the sigma2 value to 0

      if (any(is.na(sigma2) & s.nlevels == 1)) {
         sigma2[is.na(sigma2) & s.nlevels == 1] <- 0
         warning("Single-level factor(s) found in 'random' argument. Corresponding 'sigma2' value(s) fixed to 0.")
      }

      ### create model matrix for each element in mf.s

      Z.S <- vector(mode="list", length=sigma2s)

      for (j in seq_len(sigma2s)) {
         if (s.nlevels[j] == 1) {
            Z.S[[j]] <- cbind(rep(1,k))
         } else {
            if (sparse) {
               #Z.S[[j]] <- Matrix(model.matrix(~ mf.s[[j]] - 1), sparse=TRUE, dimnames=list(NULL, NULL)) ### cannot use this for factors with a single level
               Z.S[[j]] <- sparse.model.matrix(~ mf.s[[j]] - 1)
            } else {
               Z.S[[j]] <- model.matrix(~ mf.s[[j]] - 1) ### cannot use this for factors with a single level
            }
         }
         attr(Z.S[[j]], "assign")    <- NULL
         attr(Z.S[[j]], "contrasts") <- NULL
      }

   } else {

      Z.S <- NULL

   }

   #########################################################################

   ### stuff that need to be done after subsetting and filtering out NAs

   if (withR) {

      ### R may contain levels that are not in ids (that's fine; just filter them out)
      ### also, R may not be in the order that Z.S is in, so this fixes that up

      for (j in seq_along(R)) {
         if (!Rfix[j])
            next
         R[[j]] <- R[[j]][s.levels[[j]], s.levels[[j]]]
      }

      ### TODO: allow Rscale to be a vector so that different Rs can be scaled differently

      ### force each element of R to be a correlation matrix

      if (Rscale=="cor" || Rscale=="cor0")
         R[Rfix] <- lapply(R[Rfix], cov2cor)

      ### rescale R so that entries are 0 to (max(R) - min(R)) / (1 - min(R))
      ### this preserves the ultrametric properties of R and makes levels split at the root uncorrelated

      if (Rscale=="cor0")
         R[Rfix] <- lapply(R[Rfix], function(x) (x - min(x)) / (1 - min(x)))

      ### rescale R so that min(R) is zero (this is for the case that R is covariance matrix)

      if (Rscale=="cov0")
         R[Rfix] <- lapply(R[Rfix], function(x) (x - min(x)))

   }

   #########################################################################

   ### create (kxk) indicator/correlation matrices for random intercepts

   if (withS) {

      D.S <- vector(mode="list", length=sigma2s)

      for (j in seq_len(sigma2s)) {
         if (Rfix[j]) {
            if (sparse) {
               D.S[[j]] <- Z.S[[j]] %*% Matrix(R[[j]], sparse=TRUE) %*% t(Z.S[[j]])
            } else {
               D.S[[j]] <- Z.S[[j]] %*% R[[j]] %*% t(Z.S[[j]])
            }
            # D.S[[j]] <- as.matrix(nearPD(D.S[[j]])$mat)
            ### this avoids that the full matrix becomes non-positive definite but adding
            ### a tiny amount to the diagonal of D.S[[j]] is easier and works just as well
            ### TODO: consider doing something like this by default
         } else {
            D.S[[j]] <- tcrossprod(Z.S[[j]])
         }
      }

   } else {

      D.S <- NULL

   }

   #########################################################################

   ### stuff that need to be done after subsetting and filtering out NAs

   if (withG) {

      tmp <- .process.G.afterrmna(mf.g, g.nlevels, g.levels, struct[1], tau2, rho, Z.G1, Z.G2, isG=TRUE)

      mf.g <- tmp$mf.g

      g.nlevels       <- tmp$g.nlevels
      g.nlevels.f     <- tmp$g.nlevels.f
      g.levels        <- tmp$g.levels
      g.levels.f      <- tmp$g.levels.f
      g.levels.r      <- tmp$g.levels.r
      g.levels.k      <- tmp$g.levels.k
      g.levels.comb.k <- tmp$g.levels.comb.k

      tau2 <- tmp$tau2
      rho  <- tmp$rho
      G    <- tmp$G

   } else {

      g.nlevels.f     <- NULL
      g.levels.f      <- NULL
      g.levels.r      <- NULL
      g.levels.k      <- NULL
      g.levels.comb.k <- NULL

      G <- NULL

   }

   #########################################################################

   ### stuff that need to be done after subsetting and filtering out NAs

   if (withH) {

      tmp <- .process.G.afterrmna(mf.h, h.nlevels, h.levels, struct[2], gamma2, phi, Z.H1, Z.H2, isG=FALSE)

      mf.h <- tmp$mf.g

      h.nlevels       <- tmp$g.nlevels
      h.nlevels.f     <- tmp$g.nlevels.f
      h.levels        <- tmp$g.levels
      h.levels.f      <- tmp$g.levels.f
      h.levels.r      <- tmp$g.levels.r
      h.levels.k      <- tmp$g.levels.k
      h.levels.comb.k <- tmp$g.levels.comb.k

      gamma2 <- tmp$tau2
      phi    <- tmp$rho
      H      <- tmp$G

   } else {

      h.nlevels.f     <- NULL
      h.levels.f      <- NULL
      h.levels.r      <- NULL
      h.levels.k      <- NULL
      h.levels.comb.k <- NULL

      H <- NULL

   }

   #########################################################################

   #return(list(Z.S=Z.S, sigma2=sigma2, Z.G1=Z.G1, Z.G2=Z.G2, tau2=tau2, rho=rho, G=G, Z.H1=Z.H1, Z.H2=Z.H2, gamma2=gamma2, phi=phi, H=H, Rfix=Rfix, R=R))

   #########################################################################
   #########################################################################
   #########################################################################

   Y <- as.matrix(yi)

   ### initial values for variance components (need to do something better here in the future; see rma.mv2() and rma.bv() for some general ideas)

   if (verbose > 1)
      message("Extracting/computing initial values ...")

   if (verbose > 1) {
      U <- try(chol(chol2inv(chol(V))), silent=FALSE)
   } else {
      U <- try(suppressWarnings(chol(chol2inv(chol(V)))), silent=TRUE)
   }

   if (inherits(U, "try-error")) {

      total <- sigma(lm(Y ~ X - 1))^2

      #sigma2.init <- rep(.001, sigma2s)
      #tau2.init   <- rep(.001, tau2s)
      #gamma2.init <- rep(.001, gamma2s)
      #rho.init    <- rep(.50,  rhos)
      #phi.init    <- rep(.50,  phis)

      QE <- NA
      QEp <- NA

   } else {

      sX <- U %*% X
      sY <- U %*% Y
      beta.FE <- try(solve(crossprod(sX), crossprod(sX, sY)), silent=TRUE)

      if (inherits(beta.FE, "try-error"))
         stop("Cannot compute initial values.")

      ### TODO: consider a better way to set initial values
      #total <- max(.001*(sigma2s + tau2s + gamma2s), var(c(Y - X %*% res.FE$beta)) - 1/mean(1/diag(V)))
      #total <- max(.001*(sigma2s + tau2s + gamma2s), var(as.vector(sY - sX %*% beta)) - 1/mean(1/diag(V)))
      total  <- max(.001*(sigma2s + tau2s + gamma2s), var(as.vector(Y) - as.vector(X %*% beta.FE)) - 1/mean(1/diag(V)))

      QE <- sum(as.vector(sY - sX %*% beta.FE)^2)

      ### QEp calculated below

   }

   sigma2.init <- rep(total / (sigma2s + tau2s + gamma2s), sigma2s)
   tau2.init   <- rep(total / (sigma2s + tau2s + gamma2s), tau2s)
   gamma2.init <- rep(total / (sigma2s + tau2s + gamma2s), gamma2s)
   rho.init    <- rep(.50, rhos)
   phi.init    <- rep(.50, phis)

   #########################################################################

   ### set default control parameters

   con <- list(verbose = FALSE,
               optimizer = "nlminb",      # optimizer to use ("optim", "nlminb", "uobyqa", "newuoa", "bobyqa", "nloptr", "nlm", "hjk", "nmk", "ucminf")
               optmethod = "BFGS",        # argument 'method' for optim() ("Nelder-Mead" and "BFGS" are sensible options)
               sigma2.init = sigma2.init, # initial value(s) for sigma2
               tau2.init = tau2.init,     # initial value(s) for tau2
               rho.init = rho.init,       # initial value(s) for rho
               gamma2.init = gamma2.init, # initial value(s) for gamma2
               phi.init = phi.init,       # initial value(s) for phi
               REMLf = TRUE,              # full REML likelihood (including all constants)
               tol = 1e-07,               # lower bound for eigenvalues to determine if var-cov matrix is positive definite
               cholesky = ifelse(struct=="UN", TRUE, FALSE), # by default, use Cholesky factorization for G and H matrix when struct="UN" (struct has 2 elements)
               posdefify = FALSE,         # to force G and H matrix to become positive definite
               hessian = FALSE,           # to compute Hessian
               hessianCtrl=list(r=8),     # arguments passed on to 'method.args' of hessian()
               vctransf = FALSE)          # if FALSE, Hessian is computed for the untransformed (raw) variance components
                                          # if TRUE,  Hessian is computed for the transformed components (log and r-to-z space)

   ### replace defaults with any user-defined values

   con.pos <- pmatch(names(control), names(con))
   con[c(na.omit(con.pos))] <- control[!is.na(con.pos)]

   if (verbose)
      con$verbose <- verbose

   verbose <- con$verbose

   ### checks on initial values set by the user (the initial values computed by the function are replaced by the user defined ones at this point)

   if (withS && any(con$sigma2.init <= 0))
      stop("Values of 'sigma2.init' must be positive.")
   if (withG && any(con$tau2.init <= 0))
      stop("Values of 'tau2.init' must be positive.")
   if (withG && any(con$rho.init <= -1 | con$rho.init >= 1))
      stop("Values of 'rho.init' must be in (-1,1).")
   if (withH && any(con$gamma2.init <= 0))
      stop("Values of 'gamma2.init' must be positive.")
   if (withH && any(con$phi.init <= -1 | con$phi.init >= 1))
      stop("Values of 'phi.init' must be in (-1,1).")

   ### in case user manually sets con$cholesky and specifies only a single value

   if (length(con$cholesky) == 1)
      con$cholesky <- rep(con$cholesky, 2)

   ### use of Cholesky factorization only applicable for models with "UN" structure ("UNHO" may also be possible, but that still requires a fix; see below)

   if (!withG) ### in case user sets cholesky=TRUE and struct="UN" even though there is no 1st 'inner | outer' term
      con$cholesky[1] <- FALSE

   if (con$cholesky[1] && struct[1] != "UN")
      con$cholesky[1] <- FALSE

   if (!withH) ### in case user sets cholesky=TRUE and struct="UN" even though there is no 2nd 'inner | outer' term
      con$cholesky[2] <- FALSE

   if (con$cholesky[2] && struct[2] != "UN")
      con$cholesky[2] <- FALSE

   ### copy initial values back (in case they were replaced by user-defined values); those values are
   ### then shown in the 'Variance Components in Model' table that is given when verbose=TRUE; cannot
   ### replace any fixed values, since that can lead to -Inf/+Inf below when transforming the initial
   ### values and then optim() throws an error and chol(G) and/or chol(H) is then likely to fail

   #sigma2.init <- ifelse(is.na(sigma2), con$sigma2.init, sigma2)
   #tau2.init   <- ifelse(is.na(tau2), con$tau2.init, tau2)
   #rho.init    <- ifelse(is.na(rho), con$rho.init, rho)
   sigma2.init <- con$sigma2.init
   tau2.init   <- con$tau2.init
   rho.init    <- con$rho.init
   gamma2.init <- con$gamma2.init
   phi.init    <- con$phi.init

   ### plug in fixed values for sigma2, tau2, rho, gamma2, and phi and transform initial values

   con$sigma2.init <- log(sigma2.init)

   if (con$cholesky[1]) {
      G <- .con.vcov.UN(tau2.init, rho.init)
      G <- try(chol(G), silent=TRUE)
      if (inherits(G, "try-error"))
         stop("Cannot take Choleski decomposition of initial 'G' matrix.")
      con$tau2.init <- diag(G)        ### note: con$tau2.init and con$rho.init are the 'choled' values of the initial G matrix, so con$rho.init really
      con$rho.init <- G[upper.tri(G)] ### contains the 'choled' covariances; and these values are also passed on the .ll.rma.mv as the initial values
   } else {
      con$tau2.init <- log(tau2.init)
      con$rho.init  <- transf.rtoz(rho.init)
   }

   if (con$cholesky[2]) {
      H <- .con.vcov.UN(gamma2.init, phi.init)
      H <- try(chol(H), silent=TRUE)
      if (inherits(H, "try-error"))
         stop("Cannot take Choleski decomposition of initial 'H' matrix.")
      con$gamma2.init <- diag(H)      ### note: con$gamma2.init and con$phi.init are the 'choled' values of the initial H matrix, so con$phi.init really
      con$phi.init <- H[upper.tri(H)] ### contains the 'choled' covariances; and these values are also passed on the .ll.rma.mv as the initial values
   } else {
      con$gamma2.init <- log(gamma2.init)
      con$phi.init  <- transf.rtoz(phi.init)
   }

   optimizer  <- match.arg(con$optimizer, c("optim","nlminb","uobyqa","newuoa","bobyqa","nloptr","nlm","hjk","nmk","ucminf"))
   optmethod  <- match.arg(con$optmethod, c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent"))
   tol        <- con$tol
   posdefify  <- con$posdefify
   cholesky   <- con$cholesky
   optcontrol <- control[is.na(con.pos)] ### get arguments that are control arguments for optimizer

   if (length(optcontrol) == 0)
      optcontrol <- list()

   reml <- ifelse(method=="REML", TRUE, FALSE)

   ### set NLOPT_LN_BOBYQA as the default algorithm for nloptr optimizer
   ### and by default use a relative convergence criterion of 1e-8 on the function value

   if (optimizer=="nloptr" && !is.element("algorithm", names(optcontrol)))
      optcontrol$algorithm <- "NLOPT_LN_BOBYQA"

   if (optimizer=="nloptr" && !is.element("ftol_rel", names(optcontrol)))
      optcontrol$ftol_rel <- 1e-8

   #return(list(con=con, optimizer=optimizer, optmethod=optmethod, tol=tol, posdefify=posdefify, optcontrol=optcontrol))

   if (is.element(optimizer, c("uobyqa","newuoa","bobyqa"))) {
      if (!requireNamespace("minqa", quietly=TRUE))
         stop("Please install the 'minqa' package to use this optimizer.")
   }

   if (optimizer == "nloptr") {
      if (!requireNamespace("nloptr", quietly=TRUE))
         stop("Please install the 'nloptr' package to use this optimizer.")
   }

   if (is.element(optimizer, c("hjk","nmk"))) {
      if (!requireNamespace("dfoptim", quietly=TRUE))
         stop("Please install the 'dfoptim' package to use this optimizer.")
   }

   if (optimizer == "ucminf") {
      if (!requireNamespace("ucminf", quietly=TRUE))
         stop("Please install the 'ucminf' package to use this optimizer.")
   }

   ### check if length of sigma2.init, tau2.init, rho.init, gamma2.init, and phi.init matches number of variance components
   ### note: if a particular component is not included, reset (transformed) initial values (in case the user still specifies multiple initial values)

   if (withS) {
      if (length(con$sigma2.init) != sigma2s)
         stop(paste("Length of 'sigma2.init' argument (", length(con$sigma2.init), ") does not match actual number of variance components (", sigma2s, ").", sep=""))
   } else {
      con$sigma2.init <- 0
   }

   if (withG) {
      if (length(con$tau2.init) != tau2s)
         stop(paste("Length of 'tau2.init' argument (", length(con$tau2.init), ") does not match actual number of variance components (", tau2s, ").", sep=""))
   } else {
      con$tau2.init <- 0
   }

   if (withG) {
      if (length(con$rho.init) != rhos)
         stop(paste("Length of 'rho.init' argument (", length(con$rho.init), ") does not match actual number of correlations (", rhos, ").", sep=""))
   } else {
      con$rho.init <- 0
   }

   if (withH) {
      if (length(con$gamma2.init) != gamma2s)
         stop(paste("Length of 'gamma2.init' argument (", length(con$gamma2.init), ") does not match actual number of variance components (", gamma2s, ").", sep=""))
   } else {
      con$gamma2.init <- 0
   }

   if (withH) {
      if (length(con$phi.init) != phis)
         stop(paste("Length of 'phi.init' argument (", length(con$phi.init), ") does not match actual number of correlations (", phis, ").", sep=""))
   } else {
      con$phi.init <- 0
   }

   #########################################################################

   ### check whether model matrix is of full rank

   if (any(eigen(crossprod(X), symmetric=TRUE, only.values=TRUE)$values <= tol))
      stop("Model matrix not of full rank. Cannot fit model.")

   ### which variance components are fixed? (TRUE/FALSE or NA if not applicable = not included)

   if (withS) {
      sigma2.fix <- !is.na(sigma2)
   } else {
      sigma2.fix <- NA
   }
   if (withG) {
      tau2.fix <- !is.na(tau2)
      rho.fix  <- !is.na(rho)
   } else {
      tau2.fix <- NA
      rho.fix  <- NA
   }
   if (withH) {
      gamma2.fix <- !is.na(gamma2)
      phi.fix    <- !is.na(phi)
   } else {
      gamma2.fix <- NA
      phi.fix    <- NA
   }

   vc.fix <- list(sigma2=sigma2.fix, tau2=tau2.fix, rho=rho.fix, gamma2=gamma2.fix, phi=phi.fix)

   ### show which variance components are included in the model, their initial value, and their specified value (NA if not specified)

   if (verbose) {
      cat("\nVariance Components in Model:")
      if (!withS && !withG && !withH) {
         cat(" none\n\n")
      } else {
         cat("\n\n")
         vcs <- rbind(c("sigma2" = if (withS) round(sigma2.init, digits=digits) else NA,
                        "tau2"   = if (withG) round(tau2.init, digits=digits) else NA,
                        "rho"    = if (withG) round(rho.init, digits=digits) else NA,
                        "gamma2" = if (withH) round(gamma2.init, digits=digits) else NA,
                        "phi"    = if (withH) round(phi.init, digits=digits) else NA),
                        round(c(   if (withS) sigma2 else NA,
                                   if (withG) tau2 else NA,
                                   if (withG) rho else NA,
                                   if (withH) gamma2 else NA,
                                   if (withH) phi else NA), digits=digits))
         vcs <- data.frame(vcs)
         rownames(vcs) <- c("initial", "specified")
         vcs <- rbind(included=ifelse(c(rep(withS, sigma2s), rep(withG, tau2s), rep(withG, rhos), rep(withH, gamma2s), rep(withH, phis)), "Yes", "No"), fixed=unlist(vc.fix), vcs)
         print(vcs, na.print="")
         cat("\n")
      }
   }

   level <- ifelse(level > 1, (100-level)/100, ifelse(level > .5, 1-level, level))

   #return(list(sigma2s, tau2s, rhos, gamma2s, phis))

   #########################################################################
   #########################################################################
   #########################################################################

   ###### model fitting, test statistics, and confidence intervals

   if (verbose > 1)
      message("Model fitting ...")

   ### estimate sigma2, tau2, rho, gamma2, and phi as needed

   if (optimizer=="optim") {
      par.arg <- "par"
      ctrl.arg <- ", control=optcontrol"
   }
   if (optimizer=="nlminb") {
      par.arg <- "start"
      ctrl.arg <- ", control=optcontrol"
   }
   if (is.element(optimizer, c("uobyqa","newuoa","bobyqa"))) {
      par.arg <- "par"
      optimizer <- paste0("minqa::", optimizer) ### need to use this since loading nloptr masks bobyqa() and newuoa() functions
      ctrl.arg <- ", control=optcontrol"
   }
   if (optimizer=="nloptr") {
      par.arg <- "x0"
      optimizer <- paste0("nloptr::nloptr") ### need to use this due to requireNamespace()
      ctrl.arg <- ", opts=optcontrol"
   }
   if (optimizer=="nlm") {
      par.arg <- "p" ### because of this, must use argument name pX for p (number of columns in X matrix)
      ctrl.arg <- paste(names(optcontrol), unlist(optcontrol), sep="=", collapse=", ")
      if (nchar(ctrl.arg) != 0)
         ctrl.arg <- paste0(", ", ctrl.arg)
   }
   if (is.element(optimizer, c("hjk","nmk"))) {
      par.arg <- "par"
      optimizer <- paste0("dfoptim::", optimizer) ### need to use this so that the optimizers can be found
      ctrl.arg <- ", control=optcontrol"
   }
   if (optimizer=="ucminf") {
      par.arg <- "par"
      optimizer <- paste0("ucminf::ucminf") ### need to use this due to requireNamespace()
      ctrl.arg <- ", control=optcontrol"
   }

   if (method != "FE" && !is.null(random)) {

      ### if at least one parameter needs to be estimated

      if (any(is.na(c(sigma2, tau2, rho, gamma2, phi)))) {

         optcall <- paste(optimizer, "(", par.arg, "=c(con$sigma2.init, con$tau2.init, con$rho.init, con$gamma2.init, con$phi.init),
            .ll.rma.mv, reml=reml, ", ifelse(optimizer=="optim", "method=optmethod, ", ""), "Y=Y, M=V, A=NULL, X.fit=X, k=k, pX=p,
            D.S=D.S, Z.G1=Z.G1, Z.G2=Z.G2, Z.H1=Z.H1, Z.H2=Z.H2,
            sigma2.val=sigma2, tau2.val=tau2, rho.val=rho, gamma2.val=gamma2, phi.val=phi,
            sigma2s=sigma2s, tau2s=tau2s, rhos=rhos, gamma2s=gamma2s, phis=phis,
            withS=withS, withG=withG, withH=withH,
            struct=struct, g.levels.r=g.levels.r, h.levels.r=h.levels.r,
            sparse=sparse, cholesky=cholesky, posdefify=posdefify, vctransf=TRUE,
            verbose=verbose, digits=digits, REMLf=con$REMLf, dofit=FALSE", ctrl.arg, ")\n", sep="")

         #return(optcall)
         opt.res <- try(eval(parse(text=optcall)), silent=!verbose)
         #return(opt.res)

         if (inherits(opt.res, "try-error"))
            stop("Error during optimization.")

         ### convergence checks

         if (is.element(optimizer, c("optim","nlminb","dfoptim::hjk","dfoptim::nmk")) && opt.res$convergence != 0)
            stop(paste0("Optimizer (", optimizer, ") did not achieve convergence (convergence = ", opt.res$convergence, ")."))

         if (is.element(optimizer, c("minqa::uobyqa","minqa::newuoa","minqa::bobyqa")) && opt.res$ierr != 0)
            stop(paste0("Optimizer (", optimizer, ") did not achieve convergence (ierr = ", opt.res$ierr, ")."))

         if (optimizer=="nloptr::nloptr" && !(opt.res$status >= 1 && opt.res$status <= 4))
            stop(paste0("Optimizer (", optimizer, ") did not achieve convergence (status = ", opt.res$status, ")."))

         if (optimizer=="ucminf::ucminf" && !(opt.res$convergence == 1 || opt.res$convergence == 2))
            stop(paste0("Optimizer (", optimizer, ") did not achieve convergence (convergence = ", opt.res$convergence, ")."))

         if (verbose > 1) {
            cat("\n")
            print(opt.res)
         }

         ### copy estimated values to 'par' so code below works

         if (optimizer=="nloptr::nloptr")
            opt.res$par <- opt.res$solution
         if (optimizer=="nlm")
            opt.res$par <- opt.res$estimate

         if (p == k) {

            ### when fitting a saturated model (with REML estimation), estimated values of variance components can remain stuck
            ### at their initial values; this ensures that the values are fixed to zero (unless values were fixed by the user)

            sigma2[is.na(sigma2)] <- 0
            tau2[is.na(tau2)]     <- 0
            rho[is.na(rho)]       <- 0
            gamma2[is.na(gamma2)] <- 0
            phi[is.na(phi)]       <- 0

         }

      } else {

         ### if all parameter are fixed to known values, can skip optimization

         opt.res <- list(par=c(sigma2, tau2, rho, gamma2, phi))

      }

      ### save these for Hessian computation

      sigma2.val <- sigma2
      tau2.val   <- tau2
      rho.val    <- rho
      gamma2.val <- gamma2
      phi.val    <- phi

   } else {

      opt.res <- list(par=c(0,0,0,0,0))

   }

   #########################################################################

   ### do the final model fit with estimated variance components

   fitcall <- .ll.rma.mv(opt.res$par, reml=reml, Y=Y, M=V, A=A, X.fit=X, k=k, pX=p,
      D.S=D.S, Z.G1=Z.G1, Z.G2=Z.G2, Z.H1=Z.H1, Z.H2=Z.H2,
      sigma2.val=sigma2, tau2.val=tau2, rho.val=rho, gamma2.val=gamma2, phi.val=phi,
      sigma2s=sigma2s, tau2s=tau2s, rhos=rhos, gamma2s=gamma2s, phis=phis,
      withS=withS, withG=withG, withH=withH,
      struct=struct, g.levels.r=g.levels.r, h.levels.r=h.levels.r,
      sparse=sparse, cholesky=cholesky, posdefify=posdefify, vctransf=TRUE,
      verbose=FALSE, digits=digits, REMLf=con$REMLf, dofit=TRUE)

   ### extract elements

   beta <- as.matrix(fitcall$beta)
   vb   <- as.matrix(fitcall$vb)

   if (withS)
      sigma2 <- fitcall$sigma2

   if (withG) {
      G <- as.matrix(fitcall$G)
      colnames(G) <- rownames(G) <- g.levels.f[[1]]
      tau2 <- fitcall$tau2
      rho  <- fitcall$rho
   }

   if (withH) {
      H <- as.matrix(fitcall$H)
      colnames(H) <- rownames(H) <- h.levels.f[[1]]
      gamma2 <- fitcall$gamma2
      phi    <- fitcall$phi
   }

   M <- fitcall$M

   ### remove row and column names of M
   ### (but only do this if M has row/column names)

   if (!is.null(dimnames(M)))
      M <- unname(M)

   #print(M[1:8,1:8])

   ### QM calculation

   QM <- try(as.vector(t(beta)[btt] %*% chol2inv(chol(vb[btt,btt])) %*% beta[btt]), silent=TRUE)

   if (inherits(QM, "try-error"))
      QM <- NA

   rownames(beta) <- rownames(vb) <- colnames(vb) <- colnames(X)

   se <- sqrt(diag(vb))
   names(se) <- NULL
   zval <- c(beta/se)

   if (is.element(test, c("t"))) {
      dfs <- k-p
      QM  <- QM / m
      if (dfs > 0) {
         QMp  <- pf(QM, df1=m, df2=dfs, lower.tail=FALSE)
         pval <- 2*pt(abs(zval), df=dfs, lower.tail=FALSE)
         crit <- qt(level/2, df=dfs, lower.tail=FALSE)
      } else {
         QMp  <- NaN
         pval <- NaN
         crit <- NaN
      }
   } else {
      dfs  <- NA
      QMp  <- pchisq(QM, df=m, lower.tail=FALSE)
      pval <- 2*pnorm(abs(zval), lower.tail=FALSE)
      crit <- qnorm(level/2, lower.tail=FALSE)
   }

   ci.lb <- c(beta - crit * se)
   ci.ub <- c(beta + crit * se)

   #########################################################################

   ### heterogeneity test (Wald-type test of the extra coefficients in the saturated model)

   if (verbose > 1)
      message("Heterogeneity testing ...")

   QE.df <- k-p

   if (QE.df > 0L) {

      if (!is.na(QE)) {

         ### if V is not positive definite, FE model fit will fail; then QE is NA
         ### otherwise compute the RSS (which is equal to the Q/QE-test statistic)

         QEp <- pchisq(QE, df=QE.df, lower.tail=FALSE)

      }

   } else {

      ### if the user fits a saturated model, then fit must be perfect and QE = 0 and QEp = 1

      QE  <- 0
      QEp <- 1

   }

   ### log-likelihood under a saturated model with ML estimation

   ll.QE <- -1/2 * (k) * log(2*base::pi) - 1/2 * determinant(V, logarithm=TRUE)$modulus

   #########################################################################

   ###### compute Hessian

   hessian <- NA

   if (con$hessian) {

      if (verbose > 1)
         message("Computing Hessian ...")

      if (!requireNamespace("numDeriv", quietly=TRUE))
         stop("Please install the 'numDeriv' package for Hessian computation.")

      hessian <- try(numDeriv::hessian(func=.ll.rma.mv, x = if (con$vctransf) opt.res$par else c(sigma2, tau2, rho, gamma2, phi),
         method.args=con$hessianCtrl, reml=reml, Y=Y, M=V, A=NULL, X.fit=X, k=k, pX=p,
         D.S=D.S, Z.G1=Z.G1, Z.G2=Z.G2, Z.H1=Z.H1, Z.H2=Z.H2,
         sigma2.val=sigma2.val, tau2.val=tau2.val, rho.val=rho.val, gamma2.val=gamma2.val, phi.val=phi.val,
         sigma2s=sigma2s, tau2s=tau2s, rhos=rhos, gamma2s=gamma2s, phis=phis,
         withS=withS, withG=withG, withH=withH,
         struct=struct, g.levels.r=g.levels.r, h.levels.r=h.levels.r,
         sparse=sparse, cholesky=ifelse(c(con$vctransf,con$vctransf) & cholesky, TRUE, FALSE), posdefify=posdefify, vctransf=con$vctransf,
         verbose=verbose, digits=digits, REMLf=con$REMLf), silent=TRUE)

      if (inherits(hessian, "try-error")) {

         warning("Error when trying to compute Hessian.")

      } else {

         ### row/column names

         colnames(hessian) <- seq_len(ncol(hessian)) ### need to do this, so the subsetting of colnames below works

         if (sigma2s == 1) {
            colnames(hessian)[1] <- "sigma^2"
         } else {
            colnames(hessian)[1:sigma2s] <- paste("sigma^2.", seq_len(sigma2s), sep="")
         }
         if (tau2s == 1) {
            colnames(hessian)[sigma2s+1] <- "tau^2"
         } else {
            colnames(hessian)[(sigma2s+1):(sigma2s+tau2s)] <- paste("tau^2.", seq_len(tau2s), sep="")
         }
         if (rhos == 1) {
            colnames(hessian)[sigma2s+tau2s+1] <- "rho"
         } else {
            colnames(hessian)[(sigma2s+tau2s+1):(sigma2s+tau2s+rhos)] <- paste("rho.", outer(seq_len(g.nlevels.f[1]), seq_len(g.nlevels.f), paste, sep=".")[upper.tri(matrix(NA,nrow=g.nlevels.f,ncol=g.nlevels.f))], sep="")
         }
         if (gamma2s == 1) {
            colnames(hessian)[sigma2s+tau2s+rhos+1] <- "gamma^2"
         } else {
            colnames(hessian)[(sigma2s+tau2s+rhos+1):(sigma2s+tau2s+rhos+gamma2s)] <- paste("gamma^2.", seq_len(gamma2s), sep="")
         }
         if (phis == 1) {
            colnames(hessian)[sigma2s+tau2s+rhos+gamma2s+1] <- "phi"
         } else {
            colnames(hessian)[(sigma2s+tau2s+rhos+gamma2s+1):(sigma2s+tau2s+rhos+gamma2s+phis)] <- paste("phi.", outer(seq_len(h.nlevels.f[1]), seq_len(h.nlevels.f), paste, sep=".")[upper.tri(matrix(NA,nrow=h.nlevels.f,ncol=h.nlevels.f))], sep="")
         }

         rownames(hessian) <- colnames(hessian)

         ### select correct rows/columns from Hessian depending on components in the model

         #if (withS && withG && withH)
            #hessian <- hessian[1:nrow(hessian),1:ncol(hessian), drop=FALSE]
         if (withS && withG && !withH)
            hessian <- hessian[1:(nrow(hessian)-2),1:(ncol(hessian)-2), drop=FALSE]
         if (withS && !withG && !withH)
            hessian <- hessian[1:(nrow(hessian)-4),1:(ncol(hessian)-4), drop=FALSE]
         if (!withS && withG && withH)
            hessian <- hessian[2:nrow(hessian),2:ncol(hessian), drop=FALSE]
         if (!withS && withG && !withH)
            hessian <- hessian[2:(nrow(hessian)-2),2:(ncol(hessian)-2), drop=FALSE]
         if (!withS && !withG && !withH)
            hessian <- NA

      }

   }

   #########################################################################

   ###### fit statistics

   if (verbose > 1)
      message("Computing fit statistics and log likelihood ...")

   ### note: this only counts *estimated* variance components and correlations for the total number of parameters

   parms <- p + ifelse(withS, sum(ifelse(sigma2.fix,0,1)), 0) +
                ifelse(withG, sum(ifelse(tau2.fix,0,1)), 0) +
                ifelse(withG, sum(ifelse(rho.fix,0,1)), 0) +
                ifelse(withH, sum(ifelse(gamma2.fix,0,1)), 0) +
                ifelse(withH, sum(ifelse(phi.fix,0,1)), 0)

   ### note: this counts all variance components and correlations for the total number of parameters, even if they were fixed by the user or function
   #parms <- p + ifelse(withS, sigma2s, 0) + ifelse(withG, tau2s, 0) + ifelse(withG, rhos, 0) + ifelse(withH, gamma2s, 0) + ifelse(withH, phis, 0)

   ll.ML   <- fitcall$llvals[1]
   ll.REML <- fitcall$llvals[2]

   dev.ML    <- -2 * (ll.ML - ll.QE)
   AIC.ML    <- -2 * ll.ML   + 2*parms
   BIC.ML    <- -2 * ll.ML   +   parms * log(k)
   AICc.ML   <- -2 * ll.ML   + 2*parms * max(k, parms+2) / (max(k, parms+2) - parms - 1)
   dev.REML  <- -2 * (ll.REML - 0) ### saturated model has ll = 0 when using the full REML likelihood
   AIC.REML  <- -2 * ll.REML + 2*parms
   BIC.REML  <- -2 * ll.REML +   parms * log(k-p)
   AICc.REML <- -2 * ll.REML + 2*parms * max(k-p, parms+2) / (max(k-p, parms+2) - parms - 1)

   fit.stats <- matrix(c(ll.ML, dev.ML, AIC.ML, BIC.ML, AICc.ML, ll.REML, dev.REML, AIC.REML, BIC.REML, AICc.REML), ncol=2, byrow=FALSE)
   dimnames(fit.stats) <- list(c("ll","dev","AIC","BIC","AICc"), c("ML","REML"))
   fit.stats <- data.frame(fit.stats)

   #########################################################################

   ###### prepare output

   if (verbose > 1)
      message("Preparing output ...")

   p.eff <- p
   k.eff <- k

   weighted <- TRUE

   if (is.null(ddd$outlist)) {

      res <- list(b=beta, beta=beta, se=se, zval=zval, pval=pval, ci.lb=ci.lb, ci.ub=ci.ub, vb=vb,
                  sigma2=sigma2, tau2=tau2, rho=rho, gamma2=gamma2, phi=phi,
                  k=k, k.f=k.f, k.eff=k.eff, k.all=k.all, p=p, p.eff=p.eff, parms=parms, m=m,
                  QE=QE, QEp=QEp, QM=QM, QMp=QMp,
                  int.only=int.only, int.incl=int.incl, allvipos=allvipos, coef.na=coef.na,
                  yi=yi, vi=vi, V=V, W=A, X=X, yi.f=yi.f, vi.f=vi.f, V.f=V.f, X.f=X.f, W.f=W.f, ni=ni, ni.f=ni.f, M=M, G=G, H=H, hessian=hessian,
                  ids=ids, not.na=not.na, subset=subset, slab=slab, slab.null=slab.null,
                  measure=measure, method=method, weighted=weighted, test=test, dfs=dfs, btt=btt, intercept=intercept, digits=digits, level=level, sparse=sparse, control=control,
                  fit.stats=fit.stats,
                  vc.fix=vc.fix,
                  withS=withS, withG=withG, withH=withH, withR=withR,
                  sigma2s=sigma2s, tau2s=tau2s, rhos=rhos, gamma2s=gamma2s, phis=phis,
                  s.names=s.names, g.names=g.names, h.names=h.names,
                  s.levels=s.levels, s.levels.f=s.levels.f,
                  s.nlevels=s.nlevels, s.nlevels.f=s.nlevels.f,
                  g.nlevels.f=g.nlevels.f, g.nlevels=g.nlevels,
                  h.nlevels.f=h.nlevels.f, h.nlevels=h.nlevels,
                  g.levels.f=g.levels.f, g.levels.k=g.levels.k, g.levels.comb.k=g.levels.comb.k,
                  h.levels.f=h.levels.f, h.levels.k=h.levels.k, h.levels.comb.k=h.levels.comb.k,
                  struct=struct, Rfix=Rfix, R=R, Rscale=Rscale,
                  mf.r=mf.r, mf.r.f=mf.r.f, mf.g=mf.g, mf.g.f=mf.g.f, mf.h=mf.h, mf.h.f=mf.h.f, Z.S=Z.S, Z.G1=Z.G1, Z.G2=Z.G2, Z.H1=Z.H1, Z.H2=Z.H2,
                  random=random, version="2.1.0", call=mf)

   }

   if (!is.null(ddd$outlist)) {
      if (ddd$outlist == "minimal") {
         res <- list(b=beta, beta=beta, se=se, zval=zval, pval=pval, ci.lb=ci.lb, ci.ub=ci.ub, vb=vb, int.only=int.only, digits=digits, k=k, k.eff=k.eff, k.all=k.all, p=p, p.eff=p.eff, parms=parms, m=m, sigma2=sigma2, tau2=tau2, rho=rho, gamma2=gamma2, phi=phi, withR=withR, withS=withS, withG=withG, withH=withH, s.nlevels=s.nlevels, vc.fix=vc.fix, s.names=s.names, Rfix=Rfix, g.names=g.names, g.nlevels=g.nlevels, g.nlevels.f=g.nlevels.f, g.levels.k=g.levels.k, g.levels.f=g.levels.f, g.levels.comb.k=g.levels.comb.k, h.names=h.names, h.nlevels=h.nlevels, h.levels.k=h.levels.k, h.levels.f=h.levels.f, h.nlevels.f=h.nlevels.f, h.levels.comb.k=h.levels.comb.k, struct=struct, method=method, fit.stats=fit.stats, QE=QE, QEp=QEp, QM=QM, QMp=QMp, btt=btt, test=test, dfs=dfs)
      } else {
         res <- eval(parse(text=paste0("list(", ddd$outlist, ")")))
      }
   }

   class(res) <- c("rma.mv", "rma")
   return(res)

}

.set.btt <- function(btt, p, int.incl) {
  
  if (missing(btt) || is.null(btt)) {
    
    if (p > 1) {                        ### if the model matrix has more than one column
      if (int.incl) {
        btt <- seq.int(from=2, to=p)     ### and the model has an intercept term, test all coefficients except the intercept
      } else {
        btt <- seq_len(p)                ### and the model does not have an intercept term, test all coefficients
      }
    } else {
      btt <- 1                         ### if the model matrix has a single column, test that single coefficient
    }
    
  } else {
    
    ### round, take unique values, and sort
    btt <- sort(unique(round(btt)))
    
    ### check for mix of positive and negative values
    if (any(btt < 0) && any(btt > 0))
      stop("Cannot mix positive and negative 'btt' values.")
    
    ### keep/remove from 1:p vector as specified
    btt <- seq_len(p)[btt]
    
    ### (1:5)[5:6] yields c(5, NA) so remove NAs if this happens
    btt <- btt[!is.na(btt)]
    
    ### make sure that at least one valid value is left
    if (length(btt) == 0L)
      stop("Non-existent coefficients specified via 'btt'.")
    
  }
  
  return(btt)
  
}



.chkdots <- function(ddd, okargs) {
  
  for (i in seq_along(okargs))
    ddd[okargs[i]] <- NULL
  
  if (length(ddd) > 0)
    warning(paste0("Extra argument", ifelse(length(ddd) > 1, "s ", " "), "(", paste0("'", names(ddd), "'", collapse=", "), ") disregarded."), call.=FALSE)
  
}


.invcalc <- function(X, W, k) {
  
  sWX <- sqrt(W) %*% X
  res.qrs <- qr.solve(sWX, diag(k))
  #res.qrs <- try(qr.solve(sWX, diag(k)), silent=TRUE)
  #if (inherits(res.qrs, "try-error"))
  #   stop("Cannot compute QR decomposition.")
  return(tcrossprod(res.qrs))
  
}


.process.G.aftersub <- function(verbose, mf.g, struct, tau2, rho, isG, k, sparse) {
  
  if (verbose > 1)
    message(paste0("Processing '", paste0("~ ", names(mf.g)[1], " | ", names(mf.g)[2]), "' term ..."))
  
  ### get variables names in mf.g
  
  g.names <- names(mf.g) ### two names for inner and outer factor
  
  if (is.element(struct, c("CS","HCS","UN","ID","DIAG","UNHO")) && !is.factor(mf.g[[1]]) && !is.character(mf.g[[1]]))
    stop("Inner variable in (~ inner | outer) must be a factor or character variable.", call.=FALSE)
  
  ### turn each variable in mf.g into a factor (and turn the list into a data frame with 2 columns)
  ### if a variable was a factor to begin with, this drops any unused levels, but order of existing levels is preserved
  
  mf.g <- data.frame(inner=factor(mf.g[[1]]), outer=factor(mf.g[[2]]))
  
  ### check if there are any NAs anywhere in mf.g
  
  if (anyNA(mf.g))
    stop("No NAs allowed in variables specified in the 'random' argument.", call.=FALSE)
  
  ### get number of levels of each variable in mf.g (vector with two values, for the inner and outer factor)
  
  g.nlevels <- c(nlevels(mf.g[[1]]), nlevels(mf.g[[2]]))
  
  ### get levels of each variable in mf.g
  
  g.levels <- list(levels(mf.g[[1]]), levels(mf.g[[2]]))
  
  ### determine appropriate number of tau2 and rho values (care: this is done *after* subsetting)
  ### care: if g.nlevels[1] is 1, then technically there is no correlation, but we need one rho
  ### for the optimization function (this rho is fixed to 0 further in the rma.mv() function)
  
  if (struct == "CS" || struct == "ID") {
    tau2s <- 1
    rhos  <- 1
  }
  if (struct == "HCS" || struct == "DIAG") {
    tau2s <- g.nlevels[1]
    rhos  <- 1
  }
  if (struct == "UN") {
    tau2s <- g.nlevels[1]
    rhos  <- ifelse(g.nlevels[1] > 1, g.nlevels[1]*(g.nlevels[1]-1)/2, 1)
  }
  if (struct == "AR") {
    tau2s <- 1
    rhos  <- 1
  }
  if (struct == "HAR") {
    tau2s <- g.nlevels[1]
    rhos  <- 1
  }
  if (struct == "UNHO") {
    tau2s <- 1
    rhos  <- ifelse(g.nlevels[1] > 1, g.nlevels[1]*(g.nlevels[1]-1)/2, 1)
  }
  
  ### set default value(s) for tau2 if it is unspecified
  
  if (is.null(tau2))
    tau2 <- rep(NA_real_, tau2s)
  
  ### set default value(s) for rho argument if it is unspecified
  
  if (is.null(rho))
    rho <- rep(NA_real_, rhos)
  
  ### allow quickly setting all tau2 values to a fixed value
  
  if (length(tau2) == 1)
    tau2 <- rep(tau2, tau2s)
  
  ### allow quickly setting all rho values to a fixed value
  
  if (length(rho) == 1)
    rho <- rep(rho, rhos)
  
  ### check if tau2 and rho are of correct length
  
  if (length(tau2) != tau2s)
    stop(paste0("Length of ", ifelse(isG, 'tau2', 'gamma2'), " argument (", length(tau2), ") does not match actual number of variance components (", tau2s, ")."), call.=FALSE)
  if (length(rho) != rhos)
    stop(paste0("Length of ", ifelse(isG, 'rho', 'phi'), " argument (", length(rho), ") does not match actual number of correlations (", rhos, ")."), call.=FALSE)
  
  ### checks on any fixed values of tau2 and rho arguments
  
  if (any(tau2 < 0, na.rm=TRUE))
    stop(paste0("Specified value(s) of ", ifelse(isG, 'tau2', 'gamma2'), " must be non-negative."), call.=FALSE)
  if (any(rho > 1 | rho < -1, na.rm=TRUE))
    stop(paste0("Specified value(s) of ", ifelse(isG, 'rho', 'phi'), " must be in [-1,1]."), call.=FALSE)
  
  ### create model matrix for inner and outer factors of mf.g
  
  if (g.nlevels[1] == 1) {
    Z.G1 <- cbind(rep(1,k))
  } else {
    if (sparse) {
      #Z.G1 <- Matrix(model.matrix(~ mf.g[[1]] - 1), sparse=TRUE, dimnames=list(NULL, NULL))
      Z.G1 <- sparse.model.matrix(~ mf.g[[1]] - 1)
    } else {
      Z.G1 <- model.matrix(~ mf.g[[1]] - 1)
    }
  }
  if (g.nlevels[2] == 1) {
    Z.G2 <- cbind(rep(1,k))
  } else {
    if (sparse) {
      #Z.G2 <- Matrix(model.matrix(~ mf.g[[2]] - 1), sparse=TRUE, dimnames=list(NULL, NULL))
      Z.G2 <- sparse.model.matrix(~ mf.g[[2]] - 1)
    } else {
      Z.G2 <- model.matrix(~ mf.g[[2]] - 1)
    }
  }
  
  attr(Z.G1, "assign")    <- NULL
  attr(Z.G1, "contrasts") <- NULL
  attr(Z.G2, "assign")    <- NULL
  attr(Z.G2, "contrasts") <- NULL
  
  return(list(mf.g=mf.g, g.names=g.names, g.nlevels=g.nlevels, g.levels=g.levels, tau2s=tau2s, rhos=rhos, tau2=tau2, rho=rho, Z.G1=Z.G1, Z.G2=Z.G2))
  
}

.process.G.afterrmna <- function(mf.g, g.nlevels, g.levels, struct, tau2, rho, Z.G1, Z.G2, isG) {
  
  ### copy g.nlevels and g.levels
  
  g.nlevels.f <- g.nlevels
  g.levels.f  <- g.levels
  
  ### redo: turn each variable in mf.g into a factor (reevaluates the levels present, but order of existing levels is preserved)
  
  mf.g <- data.frame(inner=factor(mf.g[[1]]), outer=factor(mf.g[[2]]))
  
  ### redo: get number of levels of each variable in mf.g (vector with two values, for the inner and outer factor)
  
  g.nlevels <- c(nlevels(mf.g[[1]]), nlevels(mf.g[[2]]))
  
  ### redo: get levels of each variable in mf.g
  
  g.levels <- list(levels(mf.g[[1]]), levels(mf.g[[2]]))
  
  ### determine which levels of the inner factor were removed
  
  g.levels.r <- !is.element(g.levels.f[[1]], g.levels[[1]])
  
  ### warn if any levels were removed
  
  if (any(g.levels.r))
    warning("One or more levels of inner factor removed due to NAs.", call.=FALSE)
  
  ### for "ID" and "DIAG", fix rho to 0
  
  if (is.element(struct, c("ID","DIAG")))
    rho <- 0
  
  ### if there is only a single arm for "CS","HCS","AR","HAR" (either to begin with or after removing NAs), then fix rho to 0
  
  if (g.nlevels[1] == 1 && is.element(struct, c("CS","HCS","AR","HAR")) && is.na(rho)) {
    rho <- 0
    warning(paste0("Inner factor has only a single level, so fixed value of ", ifelse(isG, 'rho', 'phi'), " to 0."), call.=FALSE)
  }
  
  ### k per level of the inner factor
  
  g.levels.k <- table(factor(mf.g[[1]], levels=g.levels.f[[1]]))
  
  ### for "HCS","UN","DIAG","HAR": if a particular level of the inner factor only occurs once, then set corresponding tau2 value to 0 (if not already fixed)
  ### note: no longer done; variance component should still be (weakly) identifiable
  
  #if (is.element(struct, c("HCS","UN","DIAG","HAR"))) {
  #   if (any(is.na(tau2) & g.levels.k == 1)) {
  #      tau2[is.na(tau2) & g.levels.k == 1] <- 0
  #      warning("Inner factor has k=1 for one or more levels. Corresponding 'tau2' value(s) fixed to 0.")
  #   }
  #}
  
  ### create matrix where each row (= study) indicates how often each arm occurred
  ### then turn this into a list (with each element equal to a row (= study))
  
  g.levels.comb.k <- crossprod(Z.G2, Z.G1)
  g.levels.comb.k <- split(g.levels.comb.k, seq_len(nrow(g.levels.comb.k)))
  
  ### check if each study has only a single arm (could be different arms!)
  ### for "CS","HCS","AR","HAR", if yes, then must fix rho to 0 (if not already fixed)
  
  if (all(unlist(lapply(g.levels.comb.k, sum)) == 1)) {
    if (is.element(struct, c("CS","HCS","AR","HAR")) && is.na(rho)) {
      rho <- 0
      warning(paste0("Each level of the outer factor contains only a single level of the inner factor, so fixed value of ", ifelse(isG, 'rho', 'phi'), " to 0."), call.=FALSE)
    }
  }
  
  ### create matrix for each element (= study) that indicates which combinations occurred
  ### sum up all matrices (numbers indicate in how many studies each combination occurred)
  ### take upper triangle part that corresponds to the arm combinations (in order of rho)
  
  g.levels.comb.k <- lapply(g.levels.comb.k, function(x) outer(x,x, FUN="&"))
  g.levels.comb.k <- Reduce("+", g.levels.comb.k)
  g.levels.comb.k <- g.levels.comb.k[upper.tri(g.levels.comb.k)]
  
  ### UN/UNHO: if a particular combination of arms never occurs in any of the studies, then must fix the corresponding rho to 0 (if not already fixed)
  ### this also takes care of the case where each study has only a single arm
  
  if (is.element(struct, c("UN","UNHO")) && any(g.levels.comb.k == 0 & is.na(rho))) {
    rho[g.levels.comb.k == 0] <- 0
    warning(paste0("Some combinations of the levels of the inner factor never occurred. Corresponding ", ifelse(isG, 'rho', 'phi'), " value(s) fixed to 0."), call.=FALSE)
  }
  
  ### if there was only a single arm for "UN/UNHO" to begin with, then fix rho to 0
  ### (technically there is then no rho at all to begin with, but rhos was still set to 1 earlier for the optimization routine)
  ### (if there is a single arm after removing NAs, then this is dealt with below by setting tau2 and rho values to 0)
  
  if (is.element(struct, c("UN","UNHO")) && g.nlevels.f[1] == 1 && is.na(rho)) {
    rho <- 0
    warning(paste0("Inner factor has only a single level, so fixed value of ", ifelse(isG, 'rho', 'phi'), " to 0."), call.=FALSE)
  }
  
  ### construct G matrix for the various structures
  
  if (struct == "CS") {
    G <- matrix(rho*tau2, nrow=g.nlevels.f[1], ncol=g.nlevels.f[1])
    diag(G) <- tau2
  }
  
  if (struct == "HCS") {
    G <- matrix(rho, nrow=g.nlevels.f[1], ncol=g.nlevels.f[1])
    diag(G) <- 1
    G <- diag(sqrt(tau2), nrow=g.nlevels.f[1], ncol=g.nlevels.f[1]) %*% G %*% diag(sqrt(tau2), nrow=g.nlevels.f[1], ncol=g.nlevels.f[1])
    diag(G) <- tau2
  }
  
  if (struct == "UN") {
    G <- .con.vcov.UN(tau2, rho)
  }
  
  if (struct == "ID" || struct == "DIAG" ) {
    G <- diag(tau2, nrow=g.nlevels.f[1], ncol=g.nlevels.f[1])
  }
  
  if (struct == "UNHO") {
    G <- matrix(NA_real_, nrow=g.nlevels.f[1], ncol=g.nlevels.f[1])
    G[upper.tri(G)] <- rho
    G[lower.tri(G)] <- t(G)[lower.tri(G)]
    diag(G) <- 1
    G <- diag(sqrt(rep(tau2, g.nlevels.f[1])), nrow=g.nlevels.f[1], ncol=g.nlevels.f[1]) %*% G %*% diag(sqrt(rep(tau2, g.nlevels.f[1])), nrow=g.nlevels.f[1], ncol=g.nlevels.f[1])
    diag(G) <- tau2
  }
  
  if (struct == "AR") {
    if (is.na(rho)) {
      G <- matrix(NA_real_, nrow=g.nlevels.f[1], ncol=g.nlevels.f[1])
    } else {
      ### is g.nlevels.f[1] == 1 even possible here?
      if (g.nlevels.f[1] > 1) {
        G <- toeplitz(ARMAacf(ar=rho, lag.max=g.nlevels.f[1]-1))
      } else {
        G <- diag(1)
      }
    }
    G <- diag(sqrt(rep(tau2, g.nlevels.f[1])), nrow=g.nlevels.f[1], ncol=g.nlevels.f[1]) %*% G %*% diag(sqrt(rep(tau2, g.nlevels.f[1])), nrow=g.nlevels.f[1], ncol=g.nlevels.f[1])
    diag(G) <- tau2
  }
  
  if (struct == "HAR") {
    if (is.na(rho)) {
      G <- matrix(NA_real_, nrow=g.nlevels.f[1], ncol=g.nlevels.f[1])
    } else {
      ### is g.nlevels.f[1] == 1 even possible here?
      if (g.nlevels.f[1] > 1) {
        G <- toeplitz(ARMAacf(ar=rho, lag.max=g.nlevels.f[1]-1))
      } else {
        G <- diag(1)
      }
    }
    G <- diag(sqrt(tau2), nrow=g.nlevels.f[1], ncol=g.nlevels.f[1]) %*% G %*% diag(sqrt(tau2), nrow=g.nlevels.f[1], ncol=g.nlevels.f[1])
    diag(G) <- tau2
  }
  
  ### for "CS","AR","ID" set tau2 value to 0 for any levels that were removed
  
  if (any(g.levels.r) && is.element(struct, c("CS","AR","ID"))) {
    G[g.levels.r,] <- 0
    G[,g.levels.r] <- 0
  }
  
  ### for "HCS","HAR","DIAG" set tau2 value(s) to 0 for any levels that were removed
  
  if (any(g.levels.r) && is.element(struct, c("HCS","HAR","DIAG"))) {
    G[g.levels.r,] <- 0
    G[,g.levels.r] <- 0
    tau2[g.levels.r] <- 0
    warning(paste0("Fixed ", ifelse(isG, 'tau2', 'gamma2'), " to 0 for removed level(s)."), call.=FALSE)
  }
  
  ### for "UN", set tau2 value(s) and corresponding rho(s) to 0 for any levels that were removed
  
  if (any(g.levels.r) && struct == "UN") {
    G[g.levels.r,] <- 0
    G[,g.levels.r] <- 0
    tau2[g.levels.r] <- 0
    rho <- G[upper.tri(G)]
    warning(paste0("Fixed ", ifelse(isG, 'tau2', 'gamma2'), " and corresponding ", ifelse(isG, 'rho', 'phi'), " value(s) to 0 for removed level(s)."), call.=FALSE)
  }
  
  ### for "UNHO", set rho(s) to 0 corresponding to any levels that were removed
  
  if (any(g.levels.r) && struct == "UNHO") {
    G[g.levels.r,] <- 0
    G[,g.levels.r] <- 0
    diag(G) <- tau2 ### don't really need this
    rho <- G[upper.tri(G)]
    warning(paste0("Fixed ", ifelse(isG, 'rho', 'phi'), " value(s) to 0 corresponding to removed level(s)."), call.=FALSE)
  }
  
  ### special handling for the bivariate model:
  ### if tau2 (for "CS","AR","UNHO") or either tau2.1 or tau2.2 (for "HCS","UN","HAR") is fixed to 0, then rho must be fixed to 0
  
  if (g.nlevels.f[1] == 2) {
    if (is.element(struct, c("CS","AR","UNHO")) && !is.na(tau2) && tau2 == 0)
      rho <- 0
    if (is.element(struct, c("HCS","UN","HAR")) && ((!is.na(tau2[1]) && tau2[1] == 0) || (!is.na(tau2[2]) && tau2[2] == 0)))
      rho <- 0
  }
  
  #return(list(G=G, tau2=tau2, rho=rho, Z.G1=Z.G1, Z.G2=Z.G2))
  
  return(list(mf.g=mf.g, g.nlevels=g.nlevels, g.nlevels.f=g.nlevels.f, g.levels=g.levels, g.levels.f=g.levels.f, g.levels.r=g.levels.r, g.levels.k=g.levels.k, g.levels.comb.k=g.levels.comb.k, tau2=tau2, rho=rho, G=G))
  
}

### function to construct var-cov matrix for struct="UN" given vector of variances and correlations

.con.vcov.UN <- function(vars, cors) {
  dims <- length(vars)
  G <- matrix(1, nrow=dims, ncol=dims)
  G[upper.tri(G)] <- cors
  G[lower.tri(G)] <- t(G)[lower.tri(G)]
  H <- diag(sqrt(vars), nrow=dims, ncol=dims)
  return(H %*% G %*% H)
}

### function to construct var-cov matrix for struct="UN" given vector of 'choled' variances and covariances

.con.vcov.UN.chol <- function(vars, covs) {
  dims <- length(vars)
  G <- matrix(0, nrow=dims, ncol=dims)
  G[upper.tri(G)] <- covs
  diag(G) <- vars
  return(crossprod(G))
}

### function to construct var-cov matrix (G or H) for '~ inner | outer' terms

.con.E <- function(v, r, v.val, r.val, Z1, levels.r, struct, cholesky, vctransf, posdefify, sparse) {
  
  ### if cholesky=TRUE, back-transformation/substitution is done below; otherwise, back-transform and replace fixed values
  if (!cholesky) {
    if (vctransf) {
      ### variance is optimized in log-space, so exponentiate
      v <- ifelse(is.na(v.val), exp(v), v.val)
      ### correlation is optimized in r-to-z space, so use transf.ztor
      r  <- ifelse(is.na(r.val), transf.ztor(r), r.val)
    } else {
      ### for Hessian computation, can choose to leave as is
      v <- ifelse(is.na(v.val), v, v.val)
      r <- ifelse(is.na(r.val), r, r.val)
      v[v < 0] <- 0
      r[r >  1] <-  1
      r[r < -1] <- -1
    }
    v <- ifelse(v <= .Machine$double.eps*10, 0, v) ### don't do this with Cholesky factorization, since variance can be negative
  }
  
  ncol.Z1 <- ncol(Z1)
  
  if (struct == "CS") {
    E <- matrix(r*v, nrow=ncol.Z1, ncol=ncol.Z1)
    diag(E) <- v
  }
  
  if (struct == "HCS") {
    E <- matrix(r, nrow=ncol.Z1, ncol=ncol.Z1)
    diag(E) <- 1
    E <- diag(sqrt(v), nrow=ncol.Z1, ncol=ncol.Z1) %*% E %*% diag(sqrt(v), nrow=ncol.Z1, ncol=ncol.Z1)
    diag(E) <- v
  }
  
  if (struct == "UN") {
    if (cholesky) {
      E <- .con.vcov.UN.chol(v, r)
      v <- diag(E)                  ### need this, so correct values are shown when verbose=TRUE
      r <- cov2cor(E)[upper.tri(E)] ### need this, so correct values are shown when verbose=TRUE
      v[!is.na(v.val)] <- v.val[!is.na(v.val)] ### replace any fixed values
      r[!is.na(r.val)] <- r.val[!is.na(r.val)] ### replace any fixed values
    }
    E <- .con.vcov.UN(v, r)
    if (posdefify) {
      E <- as.matrix(nearPD(E)$mat) ### nearPD() in Matrix package
      v <- diag(E)                  ### need this, so correct values are shown when verbose=TRUE
      r <- cov2cor(E)[upper.tri(E)] ### need this, so correct values are shown when verbose=TRUE
    }
  }
  
  if (struct == "ID" || struct == "DIAG") {
    E <- diag(v, nrow=ncol.Z1, ncol=ncol.Z1)
  }
  
  if (struct == "UNHO") {
    E <- matrix(NA_real_, nrow=ncol.Z1, ncol=ncol.Z1)
    E[upper.tri(E)] <- r
    E[lower.tri(E)] <- t(E)[lower.tri(E)]
    diag(E) <- 1
    E <- diag(sqrt(rep(v, ncol.Z1)), nrow=ncol.Z1, ncol=ncol.Z1) %*% E %*% diag(sqrt(rep(v, ncol.Z1)), nrow=ncol.Z1, ncol=ncol.Z1)
    if (posdefify) {
      E <- as.matrix(nearPD(E, keepDiag=TRUE)$mat) ### nearPD() in Matrix package
      v <- E[1,1]                                  ### need this, so correct values are shown when verbose=TRUE
      r <- cov2cor(E)[upper.tri(E)]                ### need this, so correct values are shown when verbose=TRUE
    }
  }
  
  if (struct == "AR") {
    if (ncol.Z1 > 1) {
      E <- toeplitz(ARMAacf(ar=r, lag.max=ncol.Z1-1))
    } else {
      E <- diag(1)
    }
    E <- diag(sqrt(rep(v, ncol.Z1)), nrow=ncol.Z1, ncol=ncol.Z1) %*% E %*% diag(sqrt(rep(v, ncol.Z1)), nrow=ncol.Z1, ncol=ncol.Z1)
    diag(E) <- v
  }
  
  if (struct == "HAR") {
    if (ncol.Z1 > 1) {
      E <- toeplitz(ARMAacf(ar=r, lag.max=ncol.Z1-1))
    } else {
      E <- diag(1)
    }
    E <- diag(sqrt(v), nrow=ncol.Z1, ncol=ncol.Z1) %*% E %*% diag(sqrt(v), nrow=ncol.Z1, ncol=ncol.Z1)
    diag(E) <- v
  }
  
  ### set variance and corresponding correlation value(s) to 0 for any levels that were removed
  
  if (any(levels.r)) {
    E[levels.r,] <- 0
    E[,levels.r] <- 0
  }
  
  if (sparse)
    E <- Matrix(E, sparse=TRUE)
  
  return(list(v=v, r=r, E=E))
  
}

############################################################################

### -1 times the log likelihood (regular or restricted) for rma.mv models

.ll.rma.mv <- function(par, reml, Y, M, A, X.fit, k, pX, # note: X.fit due to hessian(); pX due to nlm()
                       D.S, Z.G1, Z.G2, Z.H1, Z.H2,
                       sigma2.val, tau2.val, rho.val, gamma2.val, phi.val,
                       sigma2s, tau2s, rhos, gamma2s, phis,
                       withS, withG, withH,
                       struct, g.levels.r, h.levels.r,
                       sparse, cholesky, posdefify, vctransf,
                       verbose, digits, REMLf, dofit=FALSE) {
  
  ### only NA values in sigma2.val, tau2.val, rho.val, gamma2.val, phi.val should be estimated; otherwise, replace with fixed values
  
  if (withS) {
    
    if (vctransf) {
      ### sigma2 is optimized in log-space, so exponentiate
      sigma2 <- ifelse(is.na(sigma2.val), exp(par[seq_len(sigma2s)]), sigma2.val)
    } else {
      ### for Hessian computation, can choose to leave as is
      sigma2 <- ifelse(is.na(sigma2.val), par[seq_len(sigma2s)], sigma2.val)
      sigma2[sigma2 < 0] <- 0
    }
    
    #if (any(is.nan(sigma2)))
    #   return(Inf)
    
    ### set really small sigma2 values equal to 0 (anything below .Machine$double.eps*10 is essentially 0)
    sigma2 <- ifelse(sigma2 <= .Machine$double.eps*10, 0, sigma2)
    
    for (j in seq_len(sigma2s)) {
      M <- M + sigma2[j] * D.S[[j]]
    }
    
  }
  
  if (withG) {
    
    resG <- .con.E(v=par[(sigma2s+1):(sigma2s+tau2s)], r=par[(sigma2s+tau2s+1):(sigma2s+tau2s+rhos)],
                   v.val=tau2.val, r.val=rho.val, Z1=Z.G1, levels.r=g.levels.r,
                   struct=struct[1], cholesky=cholesky[1], vctransf=vctransf, posdefify=posdefify, sparse=sparse)
    tau2 <- resG$v
    rho  <- resG$r
    G    <- resG$E
    
    M <- M + (Z.G1 %*% G %*% t(Z.G1)) * tcrossprod(Z.G2)
    
  }
  
  if (withH) {
    
    resH <- .con.E(v=par[(sigma2s+tau2s+rhos+1):(sigma2s+tau2s+rhos+gamma2s)], r=par[(sigma2s+tau2s+rhos+gamma2s+1):(sigma2s+tau2s+rhos+gamma2s+phis)],
                   v.val=gamma2.val, r.val=phi.val, Z1=Z.H1, levels.r=h.levels.r,
                   struct=struct[2], cholesky=cholesky[2], vctransf=vctransf, posdefify=posdefify, sparse=sparse)
    gamma2 <- resH$v
    phi    <- resH$r
    H      <- resH$E
    
    M <- M + (Z.H1 %*% H %*% t(Z.H1)) * tcrossprod(Z.H2)
    
  }
  
  ### note: if M is sparse, then using nearPD() could blow up
  
  if (posdefify)
    M <- as.matrix(nearPD(M)$mat)
  
  if (verbose > 1) {
    W <- try(chol2inv(chol(M)), silent=FALSE)
  } else {
    W <- try(suppressWarnings(chol2inv(chol(M))), silent=TRUE)
  }
  
  ### note: need W for REML llval computation
  
  if (inherits(W, "try-error")) {
    
    ### if M is not positive-definite, set the (restricted) log likelihood to -Inf
    ### this idea is based on: http://stats.stackexchange.com/q/11368/1934 (this is crude, but should
    ### move the parameter estimates away from values that create the non-positive-definite M matrix)
    
    if (dofit) {
      stop("Final variance-covariance matrix not positive definite.")
    } else {
      llval <- -Inf
    }
    
  } else {
    
    if (verbose > 1) {
      U <- try(chol(W), silent=FALSE)
    } else {
      U <- try(suppressWarnings(chol(W)), silent=TRUE)
    }
    
    ### Y ~ N(Xbeta, M), so UY ~ N(UXbeta, UMU) where UMU = I
    ### return(U %*% M %*% U)
    
    if (inherits(U, "try-error")) {
      
      if (dofit) {
        stop("Cannot fit model based on estimated marginal variance-covariance matrix.")
      } else {
        llval <- -Inf
      }
      
    } else {
      
      if (!dofit || is.null(A)) {
        
        sX   <- U %*% X.fit
        sY   <- U %*% Y
        beta <- solve(crossprod(sX), crossprod(sX, sY))
        RSS  <- sum(as.vector(sY - sX %*% beta)^2)
        if (dofit)
          vb <- matrix(solve(crossprod(sX)), nrow=pX, ncol=pX)
        
      } else {
        
        stXAX <- chol2inv(chol(as.matrix(t(X.fit) %*% A %*% X.fit)))
        #stXAX <- tcrossprod(qr.solve(sX, diag(k)))
        beta  <- matrix(stXAX %*% crossprod(X.fit,A) %*% Y, ncol=1)
        RSS   <- as.vector(t(Y - X.fit %*% beta) %*% W %*% (Y - X.fit %*% beta))
        vb    <- matrix(stXAX %*% t(X.fit) %*% A %*% M %*% A %*% X.fit %*% stXAX, nrow=pX, ncol=pX)
        
      }
      
      llvals <- c(NA_real_, NA_real_)
      
      if (dofit || !reml)
        llvals[1]  <- -1/2 * (k) * log(2*base::pi) - 1/2 * determinant(M, logarithm=TRUE)$modulus - 1/2 * RSS
      
      if (dofit || reml)
        llvals[2]  <- -1/2 * (k-pX) * log(2*base::pi) + ifelse(REMLf, 1/2 * determinant(crossprod(X.fit), logarithm=TRUE)$modulus, 0) +
        -1/2 * determinant(M, logarithm=TRUE)$modulus - 1/2 * determinant(crossprod(X.fit,W) %*% X.fit, logarithm=TRUE)$modulus - 1/2 * RSS
      
      if (dofit) {
        
        res <- list(beta=beta, vb=vb, M=M, llvals=llvals)
        
        if (withS)
          res$sigma2 <- sigma2
        
        if (withG) {
          res$G <- G
          res$tau2 <- tau2
          res$rho <- rho
        }
        
        if (withH) {
          res$H <- H
          res$gamma2 <- gamma2
          res$phi <- phi
        }
        
        return(res)
        
      } else {
        
        llval <- ifelse(reml, llvals[2], llvals[1])
        
      }
      
    }
    
  }
  
  if ((vctransf && verbose) || (!vctransf && (verbose > 1))) {
    if (withS)
      cat("sigma2 =", ifelse(is.na(sigma2), NA, paste(formatC(sigma2, digits=digits, format="f", flag=" "), " ", sep="")), "  ", sep="")
    if (withG) {
      cat("tau2 =",   ifelse(is.na(tau2),   NA, paste(formatC(tau2,   digits=digits, format="f", flag=" "), " ", sep="")), "  ", sep="")
      cat("rho =",    ifelse(is.na(rho),    NA, paste(formatC(rho,    digits=digits, format="f", flag=" "), " ", sep="")), "  ", sep="")
    }
    if (withH) {
      cat("gamma2 =", ifelse(is.na(gamma2), NA, paste(formatC(gamma2, digits=digits, format="f", flag=" "), " ", sep="")), "  ", sep="")
      cat("phi =",    ifelse(is.na(phi),    NA, paste(formatC(phi,    digits=digits, format="f", flag=" "), " ", sep="")), "  ", sep="")
    }
    cat("  ll = ", ifelse(is.na(llval), NA, formatC(llval, digits=digits, format="f", flag=" ")), sep="", "\n")
  }
  
  return(-1 * c(llval))
  
}

.is.square <- function(X) NROW(X) == NCOL(X)

.anyNAv <- function(x) {
  k <- nrow(x)
  not.na <- not.na.diag <- !is.na(diag(x))
  for (i in (1:k)[not.na.diag]) {
    not.na[i] <- !anyNA(x[i, (1:k)[not.na.diag]])
  }
  return(!not.na)
}

.is.intercept <- function(x, eps=1e-08) return(all(abs(x - 1) < eps))


.make.unique <- function(x) {
  
  x <- as.character(x)
  ux <- unique(x)
  
  for (i in 1:length(ux)) {
    xiTF <- x == ux[i]
    xi <- x[xiTF]
    if (length(xi) == 1L)
      next
    x[xiTF] <- paste(xi, seq_along(xi), sep=".")
  }
  
  return(x)
  
}