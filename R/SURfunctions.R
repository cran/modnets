#' Fit SUR models with or without constraints
#'
#' A wrapper for the \code{\link[systemfit:systemfit]{systemfit::systemfit}}
#' function that will construct formulas for all equations based on specified
#' moderators. This function was NOT designed for user-level functionality, but
#' rather exists to be embedded within \code{\link{fitNetwork}}. The purpose for
#' making it available to the user is for allowing the exact fitted model to be
#' highly customizable.
#'
#' See the \code{systemfit} package for details on customizing
#' \code{\link[systemfit:systemfit]{systemfit::systemfit}} objects. Constraints
#' can be applied via the \code{varMods} argument, which is intended to
#' facilitate the output of the \code{\link{varSelect}} and
#' \code{\link{resample}} functions. These objects can be further edited to
#' apply constraints not specified by these automated functions. Moreover, there
#' are a variety of additional arguments that can be supplied to the
#' \code{\link[systemfit:systemfit]{systemfit::systemfit}} function if desired.
#'
#' If the variable selection results from \code{\link{resample}} are intended to
#' be used as input for the \code{varMods} argument, then these results must be
#' fed into the \code{\link{modSelect}} function.
#'
#' @param data Dataframe or matrix containing idiographic temporal data.
#' @param varMods Output of \code{\link{varSelect}} or \code{\link{modSelect}}.
#'   The latter must be applied to \code{\link{resample}} results in order for
#'   it to work as input for this argument.
#' @param mod Character string. Only applies if output from
#'   \code{\link{varSelect}} or \code{\link{modSelect}} is used to constrain the
#'   model, and cross-validation \code{"CV"} was set as the criterion for
#'   model/variable selection. Options include \code{"min"}, which uses the
#'   lambda value that minimizes the objective function, or \code{"1se"} which
#'   uses the lambda value at 1 standard error above the value that minimizes
#'   the objective function.
#' @param maxiter Numeric. The maximum number of iterations to attempt before
#'   stopping the function.
#' @param m Character string or numeric value to specify the moderator (if any).
#' @param type Indicates the type of model to use, either \code{"g"} for
#'   gaussian, or \code{"c"} for categorical (i.e., binary, at present). This
#'   argument should not be edited by the user, as the appropriate input will
#'   automatically be detected.
#' @param center Logical. Determines whether to mean-center the variables.
#' @param scale Logical. Determines whether to standardize the variables.
#' @param exogenous Logical. See \code{\link{fitNetwork}} function for details.
#' @param covs something
#' @param sur Logical. Provides input to the \code{method} argument of the
#'   \code{\link[systemfit:systemfit]{systemfit::systemfit}} function. If
#'   \code{TRUE}, then the \code{method} will be \code{"SUR"}. If \code{FALSE},
#'   then the \code{method} will be \code{"OLS"}. These two methods only differ
#'   when constraints are applied. When a saturated model is fit, both methods
#'   produce the same results.
#' @param consec A logical vector that identifies which values to include in
#'   accordance with the \code{beepno} and \code{dayno} arguments in the
#'   \code{\link{fitNetwork}} function.
#' @param ... Additional arguments.
#'
#' @return A SUR model, as fit with the
#'   \code{\link[systemfit:systemfit]{systemfit::systemfit}} function.
#' @export
#'
#' @seealso \code{\link{SURnet}, \link{fitNetwork},
#'   \link[systemfit:systemfit]{systemfit::systemfit}}
SURfit <- function(data, varMods = NULL, mod = "min", maxiter = 100, m = NULL,
                   type = "g", center = TRUE, scale = FALSE, exogenous = TRUE,
                   covs = NULL, sur = TRUE, consec = NULL, ...){
  if(!is(varMods, 'list')){type <- varMods; varMods <- NULL}
  eqs <- surEqs(data = data, varMods = varMods, mod = match.arg(mod, c("min", "1se")),
                m = m, exogenous = exogenous, covs = covs)
  dat <- lagMat(data = data, type = type, m = m, covariates = covs, center = center,
                scale = scale, exogenous = exogenous, consec = consec)
  fit <- systemfit::systemfit(formula = eqs, method = ifelse(sur, "SUR", "OLS"),
                              data = cbind.data.frame(dat$Y, dat$X),
                              maxiter = maxiter, ...) # FULLFIX
  for(i in 1:length(fit$eq)){attr(fit$eq[[i]], "family") <- "gaussian"}
  names(fit$eq) <- colnames(dat$Y)
  if(!is.null(m)){
    attr(fit, "exogenous") <- exogenous
    attr(fit, "moderators") <- colnames(data)[m]
  }
  if(!is.null(covs)){attr(fit, "covariates") <- covs}
  fit
}

#' Creates temporal and contemporaneous network of SUR results
#'
#' A method for converting outputs from the
#' \code{\link[systemfit:systemfit]{systemfit::systemfit}} function into
#' temporal and contemporaneous networks. Intended as an internal function of
#' \code{\link{fitNetwork}}. Not intended for use by the user. The only purpose
#' of making it available is to allow for extreme customization, and the
#' capacity to convert any
#' \code{\link[systemfit:systemfit]{systemfit::systemfit}} output into a pair of
#' network models compatible with the \code{modnets} package.
#'
#' @param fit Output from \code{\link{SURfit}}
#' @param dat A list containing elements \code{"Y"} and \code{"X"} elements, to
#'   reflect the outcome and predictor matrices. These are lagged data matrices,
#'   and can be automatically created through the internal
#'   \code{modnets:::lagMat} function. These transformed matrices must be
#'   supplied in conjunction with the \code{\link{SURfit}} output in order to
#'   construct network models.
#' @param s Character string indicating which type of residual covariance matrix
#'   to compute for SUR models. Options include \code{"res", "dfres", "sigma"}.
#'   \code{"sigma"} uses the residual covariance matrix as computed by the
#'   \code{\link[systemfit:systemfit]{systemfit::systemfit}} function.
#'   \code{"res"} and \code{"dfres"} compute the matrix based directly on the
#'   residual values. \code{"dfres"} is the sample estimator that uses \code{N -
#'   1} in the denominator, while \code{"res"} just uses \code{N}.
#' @param m Character string or numeric value to specify the moderator (if any).
#' @param threshold See corresponding argument of \code{\link{fitNetwork}}
#' @param mval Numeric. See corresponding argument of \code{\link{fitNetwork}}
#' @param medges Numeric. See corresponding argument of \code{\link{fitNetwork}}
#' @param pcor See corresponding argument of \code{\link{fitNetwork}}
#'
#' @return Temporal and contemporaneous networks
#' @export
#'
#' @seealso \code{\link{SURfit}, \link{fitNetwork},
#'   \link[systemfit:systemfit]{systemfit::systemfit}}
SURnet <- function(fit, dat, s = "sigma", m = NULL, threshold = FALSE,
                   mval = NULL, medges = 1, pcor = "none"){
  y <- dat$Y
  p <- ncol(y)
  fitobj <- fit$eq
  yhat <- predict(fit)
  attr(fitobj, "rank") <- fit$rank
  ynames <- gsub("[.]y$", "", colnames(y))
  beta <- getCoefs(fit = fit, mat = "beta", data = dat)
  pvals <- getCoefs(fit = fit, mat = "pvals", data = dat)
  mods <- lapply(1:nrow(beta), function(z){
    model <- as.matrix(beta[z, ], ncol = 1)
    deviance <- sum((y[, z] - yhat[, z])^2)
    s <- sqrt(deviance/nrow(y))
    LL_model <- sum(dnorm(y[, z], mean = yhat[, z], sd = s, log = TRUE))
    k <- nrow(model) + 1
    aic <- (2 * k) - (2 * LL_model)
    bic <- (log(nrow(y)) * k) - (2 * LL_model)
    out <- list(deviance = deviance, LL_model = LL_model, AIC = aic, BIC = bic, model = model)
    return(out)
  })
  mods0 <- lapply(mods, '[[', "model")
  names(mods) <- names(mods0) <- names(fitobj) <- colnames(y)
  s <- match.arg(tolower(s), choices = c("sigma", "res", "dfres"))
  pcor <- ifelse(is.logical(pcor), "none", pcor)
  call <- list(type = rep("g", ncol(dat$Y)), moderators = m, mval = mval,
               lags = 1, residMat = s, threshold = threshold, pcor = pcor)
  call$mval <- mval <- ifelse(length(m) != 1, list(NULL), ifelse(
    attr(fit, "exogenous"), list(mval), list(NULL)))[[1]]
  if(!is.null(m)){
    mname <- call$moderators <- attr(fit, "moderators")
    exogenous <- call$exogenous <- attr(fitobj, "exogenous") <- attr(fit, "exogenous")
    intnames <- colnames(beta)[grep(":", colnames(beta))]
    beta2 <- beta[, intnames, drop = FALSE]
    pvals2 <- pvals[, intnames, drop = FALSE]
    interactions <- list(beta = beta2, pvals = pvals2)
    modEdges <- ifelse(pvals2 <= ifelse(!is.numeric(threshold), .05, threshold), 1, 0)
    #if(any(pvals2 == 0)){modEdges <- modEdges * ifelse(pvals2 == 0, 0, 1)}
    modEdges <- modEdges + 1
    if(length(m) == 1){
      if(is.character(m)){
        m <- which(colnames(dat$X) %in% m)
      } else if("covariates" %in% names(attributes(fit))){
        m <- m - sum(attr(fit, "covariates") < m)
      }
      ints <- lapply(mods0, function(z) rownames(z)[z[, 1] != 0][grep(":", rownames(z)[z[, 1] != 0])])
      inds0 <- unlist(lapply(ints, function(z) gsub(paste0(mname, ":|:", mname), "", z)))
      inds1 <- cbind(y = rep(names(mods), sapply(ints, length)), ints = unname(unlist(ints)))
      vars <- lapply(fitobj, vcov)
      interactions$coefvars <- data.frame(
        Y = inds1[, 1], X = unname(inds0), Z = mname, Int = inds1[, 2],
        t(sapply(seq_len(nrow(inds1)), function(i){
          vb1 <- vars[[inds1[i, 1]]][inds0[i], inds0[i]]
          vb3 <- vars[[inds1[i, 1]]][inds1[i, 2], inds1[i, 2]]
          vb1b3 <- vars[[inds1[i, 1]]][inds0[i], inds1[i, 2]]
          return(c(varX = vb1, varInt = vb3, varCov = vb1b3))
        }))
      )
      vars0 <- as.matrix(interactions$coefvars[, 1:4])
      vars1 <- interactions$coefvars[, -c(1:4)]
      if(!exogenous){
        modEdges0 <- matrix(1, p, p)
        modEdges0[, -m] <- modEdges
        modEdges0[, m] <- unname(apply(modEdges, 1, function(z) ifelse(any(z == 2), medges, 1)))
        rownames(modEdges0) <- rownames(modEdges)
        mcols0 <- character(p)
        mcols0[m] <- paste(rep(attr(fit, "moderators"), 2), collapse = ":")
        mcols0[-m] <- colnames(modEdges)
        colnames(modEdges0) <- mcols0
        modEdges <- modEdges0
      } else if(!is.null(mval)){
        margSE <- function(x, vars){
          as.numeric(sqrt(vars[1] + ((x^2) * vars[2]) + (2 * x * vars[3])))}
        ses <- getCoefs(fit = fit, mat = "ses", data = dat)
        for(i in 1:nrow(vars1)){
          B3Z <- mval * beta[vars0[i, "Y"], vars0[i, "Int"]]
          beta[vars0[i, "Y"], vars0[i, "X"]] <- beta[vars0[i, "Y"], vars0[i, "X"]] + B3Z
          ses[vars0[i, "Y"], vars0[i, "X"]] <- margSE(mval, vars1[i, ])
        }
        dfs <- matrix(rep(sapply(fitobj, '[[', "df.residual"), each = ncol(beta)),
                      nrow = nrow(beta), ncol = ncol(beta), byrow = TRUE)
        pvals <- (2 * pt(abs(beta/ses), df = dfs, lower.tail = FALSE))
        if(any(is.na(pvals))){pvals[is.na(pvals)] <- 1}
      }
    }
  } else {
    modEdges <- matrix(1, p, p)
  }
  b <- beta[, ynames, drop = FALSE]
  kappa <- solve(getCoefs(fit = fit, mat = s))
  PCC <- getCoefs(fit = fit, mat = "pcor")
  PDC <- b/(sqrt(diag(solve(kappa)) %o% diag(kappa) + b^2))
  pvals3 <- psych::corr.p(r = PCC, n = nrow(y), adjust = pcor)[[4]]
  pvals4 <- matrix(0, ncol(PCC), ncol(PCC))
  pvals4[upper.tri(pvals4)] <- pvals3[upper.tri(pvals3)]
  pvals4 <- as.matrix(Matrix::forceSymmetric(pvals4))
  dimnames(PCC) <- dimnames(kappa) <- dimnames(pvals4) <- dimnames(PDC)
  getEdgeColors <- function(adjMat){
    obj <- sign(as.vector(adjMat))
    colMat <- rep(NA, length(obj))
    if(any(obj == 1)){colMat[obj == 1] <- "darkgreen"}
    if(any(obj == 0)){colMat[obj == 0] <- "darkgrey"}
    if(any(obj == -1)){colMat[obj == -1] <- "red"}
    colMat <- matrix(colMat, ncol = ncol(adjMat), nrow = nrow(adjMat))
    colnames(colMat) <- paste0(colnames(adjMat), ifelse(
      !any(grepl("lag", colnames(adjMat))), ".lag1.", ""))
    rownames(colMat) <- rownames(adjMat)
    return(colMat)
  }
  if(threshold != FALSE){
    if(!is.character(threshold)){
      if(isTRUE(threshold)){threshold <- .05}
      b <- b * ifelse(pvals[, ynames, drop = FALSE] <= threshold, 1, 0)
      PDC <- PDC * ifelse(pvals[, ynames, drop = FALSE] <= threshold, 1, 0)
      if(!is.null(m)){interactions$beta <- beta2 * ifelse(pvals2 <= threshold, 1, 0)}
    }
    kdiag <- diag(kappa)
    kappa <- kappa * ifelse(pvals4 <= ifelse(!is.numeric(threshold), .05, threshold), 1, 0)
    PCC <- PCC * ifelse(pvals4 <= ifelse(!is.numeric(threshold), .05, threshold), 1, 0)
    diag(kappa) <- kdiag
  }
  temporal <- list(adjMat = b, edgeColors = getEdgeColors(b), modEdges = modEdges,
                   PDC = list(adjMat = PDC, edgeColors = getEdgeColors(PDC)),
                   coefs = list(beta = beta, pvals = pvals))
  if(length(m) != 1){temporal$modEdges <- NULL}
  colnames(temporal$adjMat) <- colnames(temporal$PDC$adjMat) <- colnames(temporal$edgeColors)
  contemporaneous <- list(adjMat = PCC, edgeColors = getEdgeColors(PCC), pvals = pvals4, kappa = kappa)
  colnames(contemporaneous$adjMat) <- colnames(contemporaneous$edgeColors)
  surNet <- list(call = call, temporal = temporal, contemporaneous = contemporaneous)
  if(!is.null(m)){surNet$interactions <- interactions}
  if(length(m) == 1 & ifelse(is.null(m), FALSE, exogenous)){
    madj <- modEdges0 <- matrix(0, p + 1, p + 1)
    madj[-m, -m] <- surNet$temporal$adjMat
    dimnames(madj) <- rep(list(1:(p + 1)), 2)
    rownames(madj)[-m] <- rownames(surNet$temporal$adjMat)
    rownames(madj)[m] <- paste0(mname, ".y")
    colnames(madj)[-m] <- colnames(surNet$temporal$adjMat)
    colnames(madj)[m] <- paste0(mname, ".lag1.")
    if(threshold != FALSE & !is.character(threshold)){
      madj[-m, m] <- beta[, mname] * ifelse(pvals[, mname] <= threshold, 1, 0)
    } else {
      madj[-m, m] <- beta[, mname]
    }
    modEdges0[-m, -m] <- modEdges
    modEdges0[-m, m] <- unname(apply(
      modEdges, 1, function(z){ifelse(any(z == 2), medges, 1)}))
    dimnames(modEdges0) <- dimnames(madj)
    shape <- rep("circle", p + 1)
    shape[m] <- "square"
    surNet$mnet <- list(adjMat = madj, edgeColors = getEdgeColors(madj),
                        modEdges = modEdges0, shape = shape)
  }
  surNet <- append(surNet, list(mods = mods, data = dat))
  return(surNet)
}

##### surEqs: create regression equations for fitting SUR models
surEqs <- function(data, varMods = NULL, mod = "min", m = NULL,
                   exogenous = TRUE, covs = NULL){
  if(!is.null(varMods)){
    if(attr(varMods, "criterion") != "CV"){mod <- "min"}
    mod <- match.arg(mod, c("min", "1se"))
    mod <- ifelse(mod == "min" | attr(varMods, "method") == "regsubsets", "mod0", "mod1se")
    x <- lapply(varMods, '[[', mod)
    y <- names(x)
    if("covs" %in% names(attributes(varMods))){covs <- attr(varMods, "covs")}
    if(!is.null(covs) & !is(data, 'list')){
      covs <- colnames(data)[covs]
      x <- lapply(x, function(z) c(z[!grepl(":", z)], covs, z[grepl(":", z)]))
    }
  } else {
    if(!is(data, 'list')){
      data <- lagMat(data = data, m = m, exogenous = exogenous,
                     covariates = covs, checkType = TRUE)
    }
    y <- colnames(data$Y)
    x <- rep(list(colnames(data$X)), length(y))
  }
  eqs <- lapply(seq_along(y), function(z){
    as.formula(paste(y[z], "~", paste(x[[z]], collapse = " + ")))})
  return(eqs)
}

##### getCoefs: extract beta matrix and/or sigma matrix from systemfit model
getCoefs <- function(fit, mat = "beta", data = NULL){
  if("SURfit" %in% names(fit)){fit <- fit$SURfit}
  if(class(fit) == "systemfit"){
    mat <- match.arg(arg = mat, c(
      "beta", "pvals", "ses", "sigma", "res", "dfres", "cor", "pcor"))
    if(mat %in% c("beta", "pvals", "ses")){
      ynames <- c(); n <- list()
      for(i in 1:length(fit$eq)){
        ynames[i] <- as.character(fit$eq[[i]]$terms[[2]])
        n[[i]] <- names(coef(fit$eq[[i]]))
      }
      if(!is.null(data)){
        if(class(data) == "list"){data <- data$X}
        if(is.null(colnames(data))){colnames(data) <- paste0("X", 1:ncol(data))}
        bnames <- c("(Intercept)", colnames(data))
        N <- length(bnames)
      } else {
        N <- ifelse(any(grep(":", names(coef(fit)))), length(fit$eq) * 2, length(fit$eq) + 1)
        bnames <- n[[which.max(sapply(n, length))]]
      }
      b <- t(matrix(ifelse(mat == "beta", 0, 1), ncol = length(fit$eq), nrow = N))
      rownames(b) <- ynames
      colnames(b) <- bnames
      for(i in 1:nrow(b)){
        bb <- ifelse(
          mat == "beta", list(coef(fit$eq[[i]])), ifelse(
            mat == "pvals", list(summary(fit$eq[[i]])$coefficients[, 4]),
            list(summary(fit$eq[[i]])$coefficients[, 2])))[[1]]
        for(j in 1:length(n[[i]])){b[i, names(bb)[j]] <- bb[j]}
      }
    } else {
      if(mat %in% c("res", "dfres")){
        e <- as.matrix(residuals(fit))
        N <- ifelse(mat == "res", nrow(e), nrow(e) - 1)
        b <- (t(e) %*% e)/N
      } else if(mat == "pcor"){
        b <- tryCatch({-solve(cor(residuals(fit)))}, error = function(e){
          -corpcor::pseudoinverse(cor(residuals(fit)))})
        diag(b) <- -diag(b)
        delta <- (1/sqrt(diag(b)))
        b <- t(delta * b) * delta
        diag(b) <- 0
      } else {
        b <- fit$residCov
        if(mat == "cor"){b <- cov2cor(b)}
      }
    }
  } else {
    b <- do.call(cbind, sapply(fit$mods, '[', "model"))
    colnames(b) <- rownames(fit$adjMat)
  }
  return(b)
}

##### SURsampler: Sample data for fitting system of lagged SURs
SURsampler <- function(B = NULL, S, n, seed = NULL, beta, beta2 = NULL,
                       full = TRUE, mu = 0, cholesky = FALSE,
                       time = FALSE, allDat = TRUE){
  p <- ncol(S)
  if(!is.null(B)){
    if(is(B, 'list')){B <- do.call(rbind, B)}
    if(missing(beta)){
      if(ncol(B) == (p + 1)){
        beta <- B
      } else {
        beta <- B[, 1:(p + 1)]
        beta2 <- B[, -c(1:(p + 1)), drop = FALSE]
      }
    }
  }
  if(is(beta, 'list')){beta <- do.call(rbind, beta)}
  if(dim(beta)[2] != dim(beta)[1] + 1){
    if(dim(beta)[2] == dim(beta)[1]){
      beta <- cbind(0, beta)
    } else {stop("Invalid dimensions for beta matrix")}
  }
  if(!is.null(beta2)){
    if(is(beta2, 'list')){beta2 <- do.call(rbind, beta2)}
    stopifnot(dim(beta2)[2] == dim(beta2)[1] - 1)
  }
  if(length(mu) == 1){mu <- rep(mu, p)}
  if(!is.null(seed)){set.seed(seed)}
  tx <- Sys.time()
  if(!cholesky){
    R <- mvtnorm::rmvnorm(n = n, mean = mu, sigma = S)
    #R <- MASS::mvrnorm(n = n, mu = mu, Sigma = S)
  } else {
    R <- matrix(rnorm(p * n), ncol = p) %*% chol(S)
  }
  X <- matrix(NA, ncol = p, nrow = n + 1)
  if(!is.null(seed)){set.seed(seed)}
  X[1, ] <- rnorm(p)
  it <- 0
  if(is.null(beta2)){
    repeat{
      it <- it + 1
      for(i in 1:p){
        X[it + 1, i] <- beta[i, 1] + sum(X[it, ] * beta[i, -1]) + R[it, i]
      }
      if(it == n){break}
    }
    Y <- X[-1, ]
  } else {
    Y <- X
    X2 <- X[, 1] * X[, -1]
    X <- cbind(X, X2)
    beta <- cbind(beta, beta2)
    repeat{
      it <- it + 1
      for(i in 1:p){
        Y[it + 1, i] <- beta[i, 1] + sum(X[it, ] * beta[i, -1]) + R[it, i]
        X <- cbind(Y, (Y[, 1] * Y[, -1]))
      }
      if(it == n){break}
    }
    Y <- Y[-1, ]
  }
  tx <- Sys.time() - tx
  if(time){cat("Completed in", round(tx, 3), attr(tx, "units"), "\n")}
  colnames(Y) <- paste0("X", 1:ncol(Y), ".y")
  X <- X[-nrow(X), ]
  colnames(X) <- paste0("X", 1:ncol(X))
  if(!is.null(beta2)){
    ints <- (ncol(beta) - ncol(beta2)):ncol(X)
    for(i in ints){colnames(X)[i] <- paste0("X1:X", which(ints == i) + 1)}
  }
  dat <- list(Y = Y, X = X)
  if(full){
    dat <- list(Y = Y, X = X, full = data.frame(do.call(cbind, dat)))
    if(!is.null(beta2)){colnames(dat$full)[ncol(dat$Y) + ints] <- colnames(dat$X)[ints]}
  }
  if(time){attributes(dat)$time <- tx}
  if(!allDat){
    X <- dat$X[, 1:p]
    X <- rbind(X, dat$Y[nrow(dat$Y), ])
    dat <- data.frame(X)
  }
  dat
}
