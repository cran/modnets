#' Convert continuous variables into ordinal variables
#'
#' Allows for easy conversion of continuous variables into ordinal variables.
#'
#' If a moderator value is specified via the \code{m} argument, that variable
#' will automatically be relegated to the last column of the resultant dataframe
#' or matrix. It will also be renamed "M"
#'
#' @param data An \code{n x k} dataframe or matrix containing only numeric
#'   values. Can also be a numeric vector.
#' @param m The column number or name of the moderator variable, if applicable.
#'   Leave as \code{NULL} if there is no moderator, and set to \code{TRUE} if
#'   the moderator is the last column in the matrix or dataframe.
#' @param nLevels Number of levels for the ordinal variables.
#' @param thresholds List of length \code{k}, where each element is a numeric
#'   vector of length \code{(nLevels - 1)} containing the splitpoints for
#'   grouping each variable into ordered categories.
#' @param mthresh Vector of length \code{(nLevels - 1)} containing thresholds to
#'   group values of the moderator into ordered categories.
#' @param mord if \code{FALSE}, then the moderator will not be converted into an
#'   ordinal variable (if applicable).
#' @param minOrd The minimum number of unique values allowed for each variable.
#'
#' @return A dataframe or matrix containing the ordinalized data.
#' @export
#'
#' @examples
#' dat <- data.frame(sapply(1:5, function(z) rnorm(100)))
#' ord_dat <- ordinalize(dat)
#'
#' # Including a moderator, without converting the moderator into an ordinal variable
#' ord_dat <- ordinalize(dat, m = 5, mord = FALSE)
#'
#' colnames(dat)[5] <- 'M'
#' ord_dat <- ordinalize(dat, m = 'M', mord = FALSE)
#'
#' # Use thresholds to break each variable into quartiles
#' thresh <- lapply(dat, function(z) quantile(z, probs = c(.25, .5, .75)))
#' ord_dat <- ordinalize(dat, thresholds = thresh)
ordinalize <- function(data, m = NULL, nLevels = 5, thresholds = NULL,
                       mthresh = NULL, mord = TRUE, minOrd = 3){
  if(!is(data, 'matrix') & !is(data, 'data.frame') & is(data, 'numeric')){data <- matrix(data, ncol = 1)}
  if(nLevels < minOrd){minOrd <- nLevels}
  if(!is.null(m)){
    if(isTRUE(m)){m <- ncol(data)}
    if(is.character(m)){m <- which(colnames(data) == m)}
    m0 <- m
    m <- data[, m0]
    data <- data[, -m0]
    if(mord){
      tick <- 0
      while(length(unique(mord)) < minOrd){
        thresh <- switch(2 - is.null(mthresh), rnorm(nLevels - 1), unlist(mthresh))
        mord <- as.numeric(cut(m, sort(c(-Inf, thresh, Inf))))
        tick <- tick + 1
        if(tick == 10 & !is.null(mthresh)){mthresh <- NULL}
      }
      m <- mord
    }
  }
  for(vv in 1:ncol(data)){
    tick <- ord <- 0
    while(length(unique(ord)) < minOrd){
      thresh <- switch(2 - is.null(thresholds), rnorm(nLevels - 1), thresholds[[vv]])
      ord <- as.numeric(cut(data[, vv], sort(c(-Inf, thresh, Inf))))
      tick <- tick + 1
      if(tick == 10 & !is.null(thresholds)){thresholds <- NULL}
    }
    data[, vv] <- ord
  }
  if(!is.null(m)){data$M <- m}
  return(data)
}

##### getConsec
getConsec <- function(data, beepno = NULL, dayno = NULL, makeAtt = TRUE){
  stopifnot(!is.null(beepno) & !is.null(dayno))
  stopifnot(is.data.frame(data))
  beepday <- list(beepno, dayno)
  stopifnot(sum(sapply(c(1, nrow(data)), function(i){
    all(sapply(beepday, length) == i)})) == 1)
  if(all(sapply(beepday, length) == 1)){
    if(is.character(beepno)){beepno <- which(colnames(data) == beepno)}
    if(is.character(dayno)){dayno <- which(colnames(data) == dayno)}
    data0 <- data[, -c(beepno, dayno)]
    beepday <- as.list(data[, c(beepno, dayno)])
    data <- data0
  }
  #consec <- mgm:::beepday2consec(beepvar = beepday[[1]], dayvar = beepday[[2]])
  #out <- mgm:::lagData(data = data, lags = 1, consec = consec)[[3]][-1]
  consec <- makeConsec(beepvar = beepday[[1]], dayvar = beepday[[2]])
  out <- lagData(data = data, lags = 1, consec = consec)[[3]][-1]
  if(makeAtt){
    attr(data, 'samp_ind') <- which(out)
    out <- data
  }
  return(out)
}

##### lagMat: create lagged matrices for fitting SUR models
lagMat <- function(data, type = "g", m = NULL, covariates = NULL, center = TRUE,
                   scale = FALSE, exogenous = TRUE, lags = 1,
                   consec = NULL, checkType = FALSE){
  if(lags != 1){stop("Only lag = 1 currently supported")}
  if("samp_ind" %in% names(attributes(data))){samp_ind <- attr(data, "samp_ind")}
  data <- data.frame(data)
  vs <- ynames <- colnames(data)
  ynames <- paste0(ynames, ".y")
  binary <- c()
  if(class(type) == "list" | ifelse(length(type) != ncol(data), TRUE, FALSE)){
    if(is.list(type) & "moderators" %in% names(attributes(type))){
      vm <- which(vs %in% attr(type, "moderators"))
    }
    type <- rep("gaussian", ncol(data))
    binary <- unname(which(apply(data, 2, function(z) dim(table(z)) <= 2)))
    if(length(binary) > 0){type[binary] <- "binomial"}
    if(all(type == "binary")){stop("Can't have binary outcomes for lagged models")}
  } else if(any(!type %in% c("g", "gaussian"))){
    binary <- which(!type %in% c("g", "gaussian"))
  }
  names(type) <- ynames
  if(center | scale){
    if(length(binary) > 0){
      data[, -binary] <- apply(data[, -binary], 2, scale, center, scale)
    } else {
      data <- apply(data, 2, scale, center, scale)
    }
  }
  X <- data[-nrow(data), ]
  Y <- data[-1, ]
  colnames(Y) <- ynames
  makeMods <- function(X, m){
    vs <- colnames(X)
    mods <- list()
    if(length(vs) == 2){m <- m[1]}
    for(i in seq_along(m)){
      mform <- as.formula(paste0("~ . * ", vs[m[i]]))
      mods[[i]] <- model.matrix(mform, data.frame(X))
      mods[[i]] <- mods[[i]][, grepl(":", colnames(mods[[i]]))]
    }
    mods <- do.call(cbind, mods)
    if(any(duplicated(colnames(mods)))){mods <- mods[, !duplicated(colnames(mods))]}
    if(ncol(mods) == 1){colnames(mods) <- paste0(colnames(X), collapse = ":")}
    return(mods)
  }
  if(!is.null(covariates)){
    stopifnot(class(covariates) %in% c("numeric", "integer"))
    Y <- as.matrix(Y[, -covariates])
    if(!is.null(m)){
      if(exogenous & length(union(m, covariates)) < length(ynames)){
        ynames <- ynames[-union(m, covariates)]
      } else {
        ynames <- ynames[-covariates]
      }
      covx <- X[, covariates]
      covnames <- vs[covariates]
      X <- X[, -covariates]
      m <- which(colnames(X) %in% vs[m])
      mods <- makeMods(X, m)
      mnames <- colnames(mods)
      xnames <- colnames(X)
      X <- cbind(X, covx, mods)
      colnames(X) <- c(xnames, covnames, mnames)
      if(exogenous & length(m) < ncol(Y)){Y <- Y[, -m]}
    }
  } else if(!is.null(m)){
    if(exogenous & length(m) < ncol(Y)){
      Y <- Y[, -m]
      ynames <- ynames[-m]
    }
    mods <- makeMods(X, m)
    X <- as.matrix(cbind(X, mods))
  }
  if(!is(Y, 'matrix')){
    Y <- as.matrix(Y, ncol = 1)
    colnames(Y) <- ynames
  }
  if(is.null(m) & exists("vm", inherits = FALSE) & exogenous){Y <- Y[, -vm]}
  if(exists("samp_ind", inherits = FALSE)){
    Y <- Y[samp_ind, ]
    X <- X[samp_ind, ]
  }
  rownames(Y) <- rownames(X) <- NULL
  #full <- cbind.data.frame(Y, X)
  out <- list(Y = Y, X = X) #, full = full) # FULLFIX
  if("binomial" %in% type[match(colnames(out$Y), names(type))] & checkType){
    stop("Can't have binary outcomes for lagged models")
  }
  if(!is.null(consec)){
    if(exists('samp_ind', inherits = FALSE)){consec <- consec[samp_ind]}
    out$Y <- out$Y[consec, ]
    out$X <- out$X[consec, ]
  }
  return(out)
}

##### matrixDist: compute similarity between two matrices
matrixDist <- function(mat1, mat2 = NULL, ind = "correlation", directed = TRUE,
                       similarity = TRUE, distMat = FALSE){
  if(length(mat1) == 0){return(NA)}
  if(is.null(mat2) & is(mat1, 'list')){mat2 <- mat1[[2]]; mat1 <- mat1[[1]]}
  mat1 <- as.matrix(mat1); mat2 <- as.matrix(mat2)
  if(!all(dim(mat1) == dim(mat2))){stop("Matrices must have the same dimensions")}
  d <- as.vector(mat1) - as.vector(mat2)
  if(distMat != FALSE){
    distMat <- match.arg(distMat, choices = c(0, 1, 2))
    if(distMat == 0){return(matrix(d, ncol = ncol(mat1), nrow = nrow(mat1)))}
    if(distMat == 1){return(matrix(abs(d), ncol = ncol(mat1), nrow = nrow(mat1)))}
    if(distMat == 2){return(matrix((d^2), ncol = ncol(mat1), nrow = nrow(mat1)))}
  }
  ind <- match.arg(tolower(ind), c("cosine", "mse", "rmse", "ssd", "mae", "msd", "correlation"))
  if(!directed){
    #if(ind == "cosine"){ind <- "correlation"}
    mat1 <- mat1[lower.tri(mat1)]
    mat2 <- mat2[lower.tri(mat2)]
  } else {
    mat1 <- c(mat1)
    mat2 <- c(mat2)
  }
  if(ind == "cosine"){
    #f1 <- norm(x = mat1, type = "F")
    #f2 <- norm(x = mat2, type = "F")
    #top <- sum(diag(t(mat1) %*% mat2))
    f1 <- sqrt(sum(mat1^2))
    f2 <- sqrt(sum(mat2^2))
    top <- sum(mat1 * mat2)
    return(ifelse(similarity == TRUE, top/(f1 * f2), 1 - (top/(f1 * f2))))
  }
  d <- c(mat1) - c(mat2)
  if(ind == "correlation"){
    if(all(is.na(mat1)) || all(is.na(mat2))){return(NA)}
    if(sd(mat1, na.rm = TRUE) == 0 | sd(mat2, na.rm = TRUE) == 0){return(0)}
    return(cor(mat1, mat2, use = "pairwise"))
  }
  if(ind == "ssd"){return(sum(d^2))}
  if(ind == "mse"){return(sum(d^2)/length(d))}
  if(ind == "rmse"){return(sqrt(sum(d^2)/length(d)))}
  if(ind == "mae"){return((sum(abs(d))/length(d)))}
  if(ind == "msd"){return(sum(d/length(d)))}
}

##### detrender: detrend variables based on linear model
detrender <- function(data, timevar = NULL, vars = NULL,
                      rmTimevar = TRUE, verbose = TRUE){
  if(!is(data, 'data.frame')){data <- data.frame(data)}
  data_detrend <- data
  if(is.null(timevar)){
    timevar <- 'time'
    data_detrend$time <- 1:nrow(data_detrend)
  }
  if(is.null(vars)){vars <- setdiff(colnames(data_detrend), timevar)}
  vars <- setdiff(vars, timevar)
  for(i in seq_along(vars)){
    ff <- as.formula(paste0(vars[[i]], ' ~ ', timevar))
    fit <- lm(ff, data = data_detrend)
    if(anova(fit)$P[1] < .05){
      if(verbose){message(paste0('Detrending variable ', i))}
      data_detrend[[vars[i]]][!is.na(data_detrend[[vars[i]]])] <- residuals(fit)
    }
  }
  if(rmTimevar){data_detrend <- data_detrend[, setdiff(colnames(data_detrend), timevar)]}
  return(data_detrend)
}

##### capitalize: capitalize strings
capitalize <- function(x){
  unname(sapply(x, function(z){
    z1 <- substr(z, 1, 1)
    z2 <- substr(z, 2, nchar(z))
    lower <- which(letters == z1)
    if(length(lower) != 0){z1 <- LETTERS[lower]}
    return(paste0(z1, z2))
  }))
}

##### getEdgeColors
getEdgeColors <- function(adjMat){
  obj <- sign(as.vector(adjMat))
  colMat <- rep(NA, length(obj))
  if(any(obj == 1)){colMat[obj == 1] <- "darkgreen"}
  if(any(obj == 0)){colMat[obj == 0] <- "darkgrey"}
  if(any(obj == -1)){colMat[obj == -1] <- "red"}
  colMat <- matrix(colMat, ncol = ncol(adjMat), nrow = nrow(adjMat))
  if(all(adjMat %in% 0:1)){colMat[colMat != 'darkgrey'] <- 'darkgrey'}
  dimnames(colMat) <- dimnames(adjMat)
  colMat
}

##### mmat: get interaction matrices when there are multiple exogenous moderators
mmat <- function(fit, m = NULL){
  if(isTRUE(attr(fit, 'ggm'))){
    stopifnot(!is.null(fit$call$moderators))
    allmods <- fit$call$moderators
    p <- length(fit$mods)
    vs <- varnames <- names(fit$mods)
    exind <- function(fit, mm, ex){
      lapply(lapply(lapply(fit$mods, '[[', ifelse(ex == 'b', 'model', 'pvals')), function(i){
        i[grep(mm, rownames(i)), ]
      }), function(z){
        orignames <- names(z)
        names(z) <- gsub(mm, '', names(z))
        if(any(!names(z) %in% vs)){
          orignames <- orignames[which(names(z) %in% vs)]
          z <- z[which(names(z) %in% vs)]
        }
        attr(z, 'orignames') <- orignames
        return(z)
      })
    }
    if(is.null(m)){
      m <- allmods[1]
    } else if(!tolower(m) %in% tolower(allmods)){
      if(any(grepl(tolower(m), tolower(allmods)))){
        m <- allmods[grep(tolower(m), tolower(allmods))[1]]
      } else {
        stop(paste0('No interactions with ', m))
      }
    } else if(!m %in% allmods){
      m <- allmods[which(tolower(allmods) == tolower(m))]
    }
    m1 <- paste0(':', m); m2 <- paste0(m, ':')
    mm <- paste0(c(m1, m2), collapse = '|')
    betas <- exind(fit = fit, mm = mm, ex = 'b')
    pvals <- exind(fit = fit, mm = mm, ex = 'p')
    b1 <- p1 <- structure(matrix(0, p, p), dimnames = rep(list(vs), 2))
    for(i in 1:p){
      nvi <- attr(betas[[i]], 'orignames')
      b1[i, match(names(betas[[i]]), vs)] <- betas[[i]]
      p1[i, match(names(pvals[[i]]), vs)] <- pvals[[i]]
      if(any(grepl(m2, nvi))){
        varnames[which(paste0(m2, vs) %in% nvi)] <- paste0(m2, vs)[which(paste0(m2, vs) %in% nvi)]
      }
    }
    varnames[!grepl(m2, varnames)] <- paste0(varnames[!grepl(m2, varnames)], m1)
    b1 <- t(b1)
    p1 <- t(p1); diag(p1) <- 1
    rownames(b1) <- rownames(p1) <- varnames
    out <- list(betas = b1, pvals = p1)
    attr(out, 'moderator') <- m
  } else {
    stop('Not developed yet')
  }
  return(out)
}

##### margCIs: retrieve CIs for the effect of Z on the coefficient relating X to Y
margCIs <- function(mods, data = NULL, modname = NULL, alpha = .05, nsims = 500,
                    compare = NULL, tTests = FALSE, seed = 666){
  set.seed(seed)
  if(alpha == FALSE){alpha <- .05}
  if("adjMat" %in% names(mods)){mods <- mods$mods0}
  if(is.null(modname) & "dat" %in% names(mods)){modname <- colnames(mods$dat)[ncol(mods$dat)]}
  if(is.null(data) & "dat" %in% names(mods)){data <- mods$dat[, -which(colnames(mods$dat) == modname)]}
  if("models" %in% names(mods)){mods <- mods$models}
  if(!any(grepl(":", unlist(sapply(mods, function(z) names(coef(z))))))){stop("No interaction terms in the models")}
  p <- length(mods)
  vs <- colnames(data)
  if(is.null(compare)){
    xmin <- min(mods[[1]]$model[modname])
    xmax <- max(mods[[1]]$model[modname])
  } else if(length(compare) > 1){
    xmin <- compare[1]
    xmax <- compare[2]
  }
  msims <- min_sims <- max_sims <- ci_diff <- list()
  for(i in 1:p){
    msims[[i]] <- arm::sim(mods[[i]], nsims)
    intvars <- colnames(msims[[i]]@coef)[grep(":", colnames(msims[[i]]@coef))]
    vars <- which(colnames(msims[[i]]@coef) %in% vs)
    vars <- vars[colnames(msims[[i]]@coef)[vars] %in% gsub(":.*", "", intvars)]
    min_sims[[i]] <- max_sims[[i]] <- ci_diff[[i]] <- vector("list", length = length(vars))
    if(length(intvars) != 0){
      for(j in 1:length(vars)){
        intTerm <- paste0(colnames(msims[[i]]@coef)[vars[j]], ":", modname)
        min_sims[[i]][[j]] <- msims[[i]]@coef[, vars[j]] + xmin * msims[[i]]@coef[, intTerm]
        max_sims[[i]][[j]] <- msims[[i]]@coef[, vars[j]] + xmax * msims[[i]]@coef[, intTerm]
        dm <- max_sims[[i]][[j]] - min_sims[[i]][[j]]
        if(tTests){
          se <- sqrt(2 * (((var(min_sims[[i]][[j]]) + var(max_sims[[i]][[j]]))/2)/nsims))
          tval <- mean(dm)/se
          ci_diff[[i]][[j]] <- c(quantile(dm, alpha/2), quantile(dm, 1 - alpha/2), b = mean(dm), se = se,
                                 t = tval, p = 2 * pt(abs(tval), df = mods[[i]]$df.residual, lower.tail = F))
        } else {
          ci_diff[[i]][[j]] <- c(b = mean(dm), quantile(dm, alpha/2), quantile(dm, 1 - alpha/2))
        }
      }
      names(min_sims[[i]]) <- names(max_sims[[i]]) <- names(ci_diff[[i]]) <- colnames(msims[[i]]@coef)[vars]
    }
  }
  ci_diff <- lapply(ci_diff, function(z) do.call(rbind, z))
  names(ci_diff) <- vs
  attributes(ci_diff)$moderator <- modname
  if(!is.null(compare)){attributes(ci_diff)$compare <- compare[1:2]}
  return(ci_diff)
}

##### getInts: retrieve interactions effects; e.g., which ones apply to both variables
getInts <- function(x, allInts = FALSE, getNames = FALSE, ...){
  if("adjMat" %in% names(x)){x <- x$mods0}
  if("models" %in% names(x)){x <- margCIs(x, ...)}
  ones <- lapply(x, function(z) as.numeric(!(z[, 2] < 0 & z[, 3] > 0)))
  newx <- lapply(1:length(x), function(z) data.frame(x[[z]], ones[[z]]))
  for(i in 1:length(newx)){
    if(nrow(newx[[i]]) >= 1){
      if(any(newx[[i]] == 1)){
        newx[[i]] <- rownames(newx[[i]][newx[[i]][, 4] == 1, ])
      } else {
        newx[[i]] <- NA
      }
    } else {
      newx[[i]] <- NA
    }
  }
  ints <- matrix(0, length(x), length(x))
  for(i in 1:length(x)){ints[i, match(newx[[i]], names(x))] <- 1}
  dimnames(ints) <- rep(list(names(x)), 2)
  zints <- t(ints) * ints
  if(allInts == FALSE & any(zints == 1)){
    f1 <- function(x){
      n <- ncol(x) - 1
      vals <- function(z){
        z1 <- (z * (z - 1))/2 + z + 1
        z1:(z1 + z)
      }
      y <- list()
      for(i in 1:n){y[[i]] <- vals(i - 1)}
      lapply(y, function(z) x[upper.tri(x)][z])
    }
    f2 <- function(x){
      vs <- colnames(x)
      x2 <- f1(x)
      x3 <- list()
      for(i in 1:length(x2)){
        if(!any(x2[[i]] == 1)){
          x3[[i]] <- NA
        } else {
          x3[[i]] <- vs[which(x2[[i]] == 1)]
        }
      }
      x3 <- as.list(unlist(x3[!is.na(x3)]))
      x4 <- rep(2:(length(x2) + 1), sapply(x2, sum))
      lapply(1:length(x3), function(z) c(vs[x4[z]], x3[[z]]))
    }
    if(getNames){return(f2(zints))}
  }
  if(allInts){return(ints)} else {return(zints)}
}

##### condEffects: generate values of beta for X conditioned on certain values of Z
condEffects <- function(mods, x = NULL, xn = NULL, data = NULL, alpha = .05,
                        adjCI = FALSE, saveMods = TRUE){
  if("SURnet" %in% c(names(mods), names(attributes(mods)))){
    if(!"SURnet" %in% names(mods)){stop("Need 'SURfit' mods")}
    if(!"mnet" %in% names(attributes(mods))){stop("Must have only one exogenous moderator")}
    fitobj <- mods$SURfit$eq
    net <- mods$SURnet
    ynames <- names(net$mods)
    mname <- net$call$moderators
    vars <- net$interactions$coefvars
    vars0 <- as.matrix(vars[, 5:7])
    vars <- as.matrix(vars[, 1:4])
    dfs <- unname(sapply(fitobj, '[[', "df.residual")[match(vars[, 1], ynames)])
    xr <- range(net$data$X[, mname])
    if(length(xn) == 1){x <- seq(xr[1], xr[2], length.out = xn)}
    if(is.null(x)){x <- seq(xr[1], xr[2], length.out = 100)}
    mats <- lapply(seq_along(x), function(i){
      out <- SURnet(fit = mods$SURfit, dat = net$data, m = mname, mval = x[i])
      out <- out$temporal$adjMat
      colnames(out) <- gsub("[.]lag1[.]$", "", colnames(out))
      return(out)
    })
    margSE <- function(x, vars){sqrt(vars[1] + ((x^2) * vars[2]) + (2 * x * vars[3]))}
    dats <- lapply(seq_len(nrow(vars)), function(i){
      D <- data.frame(x, matrix(NA, ncol = 4, nrow = length(x)))
      D[, 2] <- sapply(mats, function(z) z[vars[i, 1], vars[i, 2]])
      D[, 3] <- margSE(x = x, vars = vars0[i, ])
      outlog <- capture.output({fdr <- interactionTest::fdrInteraction(
        D[, 2], D[, 3], df = dfs[i], level = 1 - alpha)})
      D[, 4] <- D[, 2] - (ifelse(adjCI, fdr, qnorm(1 - alpha/2)) * D[, 3])
      D[, 5] <- D[, 2] + (ifelse(adjCI, fdr, qnorm(1 - alpha/2)) * D[, 3])
      colnames(D) <- c("x", "y", "se", "lower", "upper")
      return(D)
    })
    yints <- unique(vars[, 1])
    Y <- lapply(seq_along(yints), function(i){
      whichy <- which(vars[, 1] == yints[i])
      out <- dats[whichy]
      names(out) <- vars[whichy, 2]
      return(out)
    })
    names(Y) <- gsub("[.]y$", "", yints)
    attributes(Y)[c("moderator", "alpha", "SURnet")] <- list(mname, alpha, TRUE)
    if(saveMods){attributes(Y)$mods <- net$data$X[, mname]}
  } else {
    if("adjMat" %in% names(mods)){mods <- mods$mods0}
    if(!any(grepl(":", unlist(sapply(mods$models, function(z) names(coef(z))))))){stop("No interaction terms in the models")}
    if(is.null(data)){data <- mods$dat[, -ncol(mods$dat)]}
    xr <- range(mods$dat[, ncol(mods$dat)])
    if(!is.null(xn)){
      stopifnot(length(xn) == 1)
      x <- seq(xr[1], xr[2], length.out = xn)
    }
    if(is.null(x)){
      x <- ifelse(dim(table(mods$dat[, ncol(mods$dat)])) <= 2, list(c(0, 1)), ifelse(nrow(data) > 100, 100, nrow(data)))[[1]]
      if(length(x) == 1){x <- seq(xr[1], xr[2], length.out = x)}
    }
    mats <- list()
    for(i in seq_along(x)){mats[[i]] <- t(modNet(mods, data, mval = x[i], nsims = 0)$nodewise$adjNW)}
    net0 <- modNet(mods, data)
    p <- ncol(net0$adjMat)
    vars <- net0$interactions$coefvars
    if("varMods" %in% names(attributes(mods$models))){
      inds <- net0$interactions$inds
      inds[, 2] <- match(gsub(":.*", "", inds[, 2]), colnames(data))
      n <- nrow(inds)
      df <- sapply(mods$models, function(z) z$df.residual)[inds[, 1]]
    } else {
      n <- (p * (p - 1))
      inds1 <- net0$interactions$inds
      inds2 <- cbind(inds1[, 2], inds1[, 1])
      inds <- rbind(inds1, inds2)
      inds <- inds[order(inds[, 2]), ]
      inds <- inds[order(inds[, 1]), ]
      inds3 <- cbind(inds[, 1], rep(c(1:(p - 1)), p))
      df <- rep(nrow(data) - (2 * p), n)
    }
    margSE <- function(x, vars){sqrt(vars[1] + ((x^2) * vars[2]) + (2 * x * vars[3]))}
    fdr <- c()
    dats <- list()
    for(i in 1:n){
      dats[[i]] <- data.frame(matrix(NA, ncol = 5, nrow = length(x)))
      dats[[i]][, 1] <- x
      dats[[i]][, 2] <- sapply(mats, function(z) z[inds[i, 1], inds[i, 2]])
      if(n < (p * (p - 1))){
        dats[[i]][, 3] <- margSE(x = x, vars = vars[[i]])
      } else {
        dats[[i]][, 3] <- margSE(x = x, vars = vars[[inds3[i, 1]]][[inds3[i, 2]]])
      }
      log <- capture.output({fdr[i] <- interactionTest::fdrInteraction(dats[[i]][, 2], dats[[i]][, 3], df = df[i], level = 1 - alpha)})
      dats[[i]][, 4] <- dats[[i]][, 2] - (ifelse(adjCI, fdr[i], qnorm(1 - alpha/2)) * dats[[i]][, 3])
      dats[[i]][, 5] <- dats[[i]][, 2] + (ifelse(adjCI, fdr[i], qnorm(1 - alpha/2)) * dats[[i]][, 3])
      colnames(dats[[i]]) <- c("x", "y", "se", "lower", "upper")
    }
    Y <- list()
    for(i in 1:p){
      if(i %in% inds[, 1]){
        Y[[i]] <- dats[which(inds[, 1] == i)]
        if(n < (p * (p - 1))){
          names(Y[[i]]) <- colnames(data)[inds[inds[, 1] == i, 2]]
        } else {
          names(Y[[i]]) <- names(vars[[i]])
        }
      } else {
        Y[[i]] <- NA
      }
    }
    names(Y) <- colnames(data)
    attributes(Y)$moderator <- colnames(mods$dat)[ncol(mods$dat)]
    attributes(Y)$alpha <- alpha
    if(saveMods){attributes(Y)$mods <- mods}
  }
  return(Y)
}

##### checkInclude
checkInclude <- function(x, which.net = "temporal"){
  if("SURnet" %in% c(names(x), names(attributes(x)))){
    if("SURnet" %in% names(x)){x <- x$SURnet}
    if(is.numeric(which.net)){which.net <- c("t", "c", "p")[which.net]}
    which.net <- match.arg(tolower(which.net), c("temporal", "contemporaneous", "pdc"))
    x <- x[[ifelse(which.net == "contemporaneous", "contemporaneous", "temporal")]]
    if(which.net == "pdc"){x <- x$PDC}
  }
  if("adjMat" %in% names(x)){x <- t(x$adjMat)}
  directed <- !isTRUE(all.equal(x, t(x), check.attributes = FALSE))
  weighted <- any(!x %in% c(0, 1))
  include <- c("Degree", "Strength", "OutDegree", "InDegree", "OutStrength",
               "InStrength", "Closeness", "Betweenness", "ExpectedInfluence",
               "OutExpectedInfluence", "InExpectedInfluence")
  if(directed){
    include <- c("OutStrength", "InStrength", "Closeness", "Betweenness",
                 "OutExpectedInfluence", "InExpectedInfluence")
  } else {
    include <- c("Strength", "Closeness", "Betweenness", "ExpectedInfluence")
  }
  if(!weighted){include <- gsub("Strength", "Degree", include)}
  return(include)
}

##### setupVAR: prepare data for analysis
setupVAR <- function(data, idvar = NULL, method = c("gvar", "lmer", "all"),
                     center = TRUE, scale = TRUE, vars = NULL,
                     centerWithin = TRUE, scaleWithin = FALSE){
  method <- match.arg(method)
  if(class(data) == "list"){
    if(!"data" %in% names(data)){stop("Must supply data frame")}
    data <- as.data.frame(data[["data"]])
  }
  if(is.null(idvar)){
    if(!any(grepl("ID", colnames(data)))){stop("Must supply 'idvar'")}
    idvar <- colnames(data)[grep("ID", colnames(data))]
  }
  data <- data.frame(data[, -which(colnames(data) == idvar)], ID = data[, idvar])
  if(is.null(vars)){vars <- colnames(data)[!colnames(data) %in% "ID"]}
  ids <- unique(data[, "ID"])
  binary <- apply(data[, -which(colnames(data) == "ID")], 2, function(z) length(table(z)) <= 2)
  binary <- ifelse(any(binary), list(names(which(binary))), list(NULL))[[1]]
  vars0 <- setdiff(vars, binary)
  if(center){data[, vars0] <- apply(data[, vars0], 2, scale, center, scale)}
  dataByID <- lapply(ids, function(z) data[data[, "ID"] == z, seq_along(vars)])
  N <- rep(seq_along(ids), (sapply(dataByID, nrow) - ifelse(method == "all", 0, 1)))
  dataMeans <- do.call(rbind, lapply(dataByID, colMeans))[N, ]
  colnames(dataMeans) <- paste0(colnames(dataMeans), ".m")
  if(method != "lmer"){
    data0 <- lapply(dataByID, function(z){if(centerWithin){
      z[, vars0] <- apply(z[, vars0], 2, scale, centerWithin, scaleWithin)}
      return(z)
    })
    if(method == "all"){return(data.frame(do.call(rbind, data0), dataMeans, ID = ids[N]))}
    Y <- do.call(rbind, lapply(data0, function(z) z[-1, ]))
    X <- do.call(rbind, lapply(data0, function(z) z[-nrow(z), ]))
  } else {
    Y <- do.call(rbind, lapply(dataByID, function(z) z[-1, ]))
    X <- do.call(rbind, lapply(dataByID, function(z){
      z <- z[-nrow(z), ]
      if(centerWithin){z[, vars0] <- apply(
        z[, vars0], 2, scale, centerWithin, scaleWithin)}
      return(z)
    }))
  }
  dat <- data.frame(Y, X, dataMeans, ID = ids[N])
  dat
}

##### trevSimulateVAR: core single-subject VAR sampler
trevSimulateVAR <- function(parms, means = 0, lags = 1, Nt = 100, init,
                            residuals = 0.1, burnin, m = NULL, mb1 = NULL,
                            mb2 = NULL, mcenter = TRUE, skewErr = FALSE){
  if(is.matrix(parms)){parms <- list(parms)} ##### START 3
  if(any(sapply(parms, function(x) length(unique(dim(x))) > 1))){
    stop("non-square graph detected.")
  }
  if(missing(burnin)){burnin <- min(round(Nt/2), 100)}
  Ni <- ncol(parms[[1]])
  if(length(means) == 1){means <- rep(means, Ni)}
  maxLag <- max(lags)
  if(length(residuals) == 1){
    residuals <- diag(residuals, Ni)
  } else if(length(residuals) == Ni){
    residuals <- diag(residuals)
  }
  if(!is.matrix(residuals) && ncol(residuals) != Ni && nrow(residuals) != Ni){
    stop("'residuals' is not a square matrix")
  }
  totTime <- Nt + burnin
  if(missing(init)){init <- matrix(0, maxLag, Ni)}
  Res <- matrix(NA, totTime, Ni)
  if(!is.null(m)){
    if(!is(m, "matrix")){
      if(mcenter){m <- m - mean(m)}
      m <- matrix(rep(m, ncol(Res)), ncol = ncol(Res))
    }
  }
  skewErr <- ifelse(is.numeric(skewErr), skewErr, ifelse(isTRUE(skewErr), 3, FALSE))
  Res[1:maxLag, ] <- init ##### STOP 3
  for(t in (maxLag + 1):(totTime)){
    if(!is.null(m)){
      x1 <- rowSums(do.call(cbind, lapply(seq_along(lags), function(i){
        parms[[i]] %*% (Res[t - lags[i], ] - means)})))
      x2 <- mb1 * m[t - 1, ]
      x3 <- rowSums(mb2 %*% ((Res[t - 1, ] - means) * m[t - 1, ]))
      Res[t, ] <- means + x1 + x2 + x3
    } else {
      Res[t, ] <- means + rowSums(do.call(cbind, lapply(seq_along(lags),function(i){
        parms[[i]] %*% (Res[t - lags[i], ] - means)})))
    }
    e <- switch(
      2 - is.numeric(skewErr),
      sn::rmsn(1, rep(0, Ni), residuals, rep(skewErr, Ni)),
      mvtnorm::rmvnorm(1, rep(0, Ni), residuals)
    )
    Res[t, ] <- Res[t, ] + e
  }
  out <- as.data.frame(Res[-(1:burnin), ])
  return(out)
}

##### simPcor: creates network matrices
simPcor <- function(Nvar, sparsity = 0.5, parRange = c(0.5, 1), constant = 1.5,
                    propPos = 0.5, precision = FALSE, finalRange = NULL){
  trueKappa <- matrix(0, Nvar, Nvar)
  kupper <- upper.tri(trueKappa)
  klower <- lower.tri(trueKappa)
  totEdges <- sum(kupper)
  trueKappa[kupper][sample(seq_len(totEdges), round((1 - sparsity) * totEdges))] <- 1
  vals <- sample(c(-1, 1), totEdges, TRUE, prob = c(propPos, 1 - propPos)) * runif(totEdges, min(parRange), max(parRange))
  trueKappa[kupper] <- trueKappa[kupper] * vals
  trueKappa[klower] <- t(trueKappa)[klower]
  diag(trueKappa) <- constant * rowSums(abs(trueKappa))
  diag(trueKappa) <- ifelse(diag(trueKappa) == 0, 1, diag(trueKappa))
  trueKappa <- trueKappa/diag(trueKappa)[row(trueKappa)]
  trueKappa <- (trueKappa + t(trueKappa))/2
  if(!precision){trueKappa <- as.matrix(qgraph::wi2net(trueKappa))}
  if(!is.null(finalRange)){
    if(max(abs(trueKappa)) > max(finalRange)){
      trueKappa <- trueKappa/(max(abs(abs(trueKappa)))/max(finalRange))
    }
    if(min(abs(trueKappa[trueKappa != 0])) < min(finalRange)){
      wmin <- abs(trueKappa[trueKappa != 0]) < min(finalRange)
      while(min(abs(trueKappa[trueKappa != 0])) < min(finalRange)){
        trueKappa[trueKappa != 0][wmin] <- trueKappa[trueKappa != 0][wmin] * constant
      }
    }
  }
  return(trueKappa)
}

##### pcor2: creates partial correlation matrix
pcor2 <- function(x){
  x <- -corpcor::pseudoinverse(x)
  diag(x) <- -diag(x)
  x <- cov2cor(x) - diag(ncol(x))
  return((x + t(x))/2)
}


### -------------------------- QGRAPH FUNCTIONS ---------------------------- ###

##### cor0: deal with sd == 0
cor0 <- function(x, y){
  if(all(is.na(x)) | all(is.na(y))){return(NA)}
  if(sd(x, na.rm = TRUE) == 0 | sd(y, na.rm = TRUE) == 0){return(0)}
  return(cor(x, y, use = 'pairwise.complete.obs'))
}

##### fnames: rename
fnames <- function(x, n = ''){
  if(is.null(names(x))){names(x) <- paste0(n, seq_along(x))}
  out <- ifelse(names(x) == '', paste0(n, seq_along(x)), names(x))
  return(out)
}

##### scaleNA: deal with NAs
scaleNA <- function(x){
  if(all(is.na(x))){
    out <- NA
  } else if(sd(x, na.rm = TRUE) != 0){
    out <- (x - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)
  } else {
    out <- rep(0, length(x))
  }
  return(out)
}

##### WS: clustering method
WS <- function(x, threshWS = 0){
  W <- qgraph::getWmat(x)
  if(is.list(W)){
    out <- lapply(W, WS, threshWS = threshWS)
  } else {
    thresh <- threshWS
    diag(W) <- 0
    A <- matrix(0, nrow = nrow(W), ncol = ncol(W))
    A[W > thresh] <- 1
    A[W < (-thresh)] <- -1
    diag(A) <- 0
    absA <- abs(A)
    abs_top <- diag(absA %*% absA %*% absA)
    top <- diag(A %*% A %*% A)
    k <- colSums(absA)
    bottom <- k * (k - 1)
    CW <- top/bottom
    absCW <- abs_top/bottom
    out <- data.frame(cbind(clustWS = absCW, signed_clustWS = CW))
  }
  return(out)
}

##### zhang: clustering method
zhang <- function(x){
  W <- qgraph::getWmat(x)
  if(is.list(W)){
    out <- lapply(W, zhang)
  } else {
    diag(W) <- 0
    absW <- abs(W)
    top <- diag(W %*% W %*% W)
    abs_top <- diag(absW %*% absW %*% absW)
    bottom <- colSums(absW)^2 - colSums(W^2)
    Z <- top/bottom
    absZ <- abs_top/bottom
    out <- data.frame(cbind(clustZhang = absZ, signed_clustZhang = Z))
  }
  return(out)
}

##### onnela: clustering method
onnela <- function(x, threshON = 0){
  W <- qgraph::getWmat(x)
  if(is.list(W)){
    out <- lapply(W, onnela, threshON = threshON)
  } else {
    thresh <- threshON
    diag(W) <- 0
    W[abs(W) < thresh] <- 0
    absW <- abs(W)
    absW13 <- absW^(1/3)
    W13 <- matrix(nrow = nrow(W), ncol = ncol(W))
    W13[W >= 0] <- absW13[W >= 0]
    W13[W < 0] <- (abs(W[W < 0])^(1/3)) * (-1)
    top <- diag(W13 %*% W13 %*% W13)
    abs_top <- diag(absW13 %*% absW13 %*% absW13)
    A <- matrix(0, nrow = nrow(W), ncol = ncol(W))
    A[abs(W) > thresh] <- 1
    k <- colSums(A)
    bottom <- k * (k - 1)
    on <- top/bottom
    abs_on <- abs_top/bottom
    out <- data.frame(cbind(clustOnnela = abs_on, signed_clustOnnela = on))
  }
  return(out)
}

### ---------------------------  FROM BOOTNET ------------------------------ ###

##### netsimulator is adapted from bootnet:::plot.netSimulator
netsimulator <- function(x, xvar = "factor(nCases)", yvar = c("sensitivity", "specificity", "correlation"),
                         xfacet = "measure", yfacet = ".", color = NULL, ylim = c(0, 1), print = TRUE,
                         xlab = "Number of cases", ylab, outlier.size = 0.5, boxplot.lwd = 0.5,
                         style = c("fancy", "basic"), ...){
  style <- match.arg(style)
  if(xvar != "factor(nCases)" && xlab == "Number of cases"){
    warning("argument 'xvar' is not 'factor(nCases)' while argument 'xlab' is still 'Number of cases'. X-axis label might be wrong.")
  }
  if(missing(ylab)){
    if(xfacet != "measure"){
      ylab <- paste(yvar, collapse = "; ")
    } else {
      ylab <- ""
    }
  }
  #Gathered <- tidyr::gather_(x, "measure", "value", yvar)
  #Gathered <- x %>% tidyr::gather_("measure", "value", yvar)
  xx <- do.call(rbind, replicate(length(yvar), x, FALSE))[, setdiff(colnames(x), yvar)]
  Gathered <- data.frame(xx, setNames(stack(x[, yvar, drop = FALSE])[, 2:1], c('measure', 'value')))
  Gathered$measure <- as.character(Gathered$measure)
  if(!is.null(color)){
    Gathered[[color]] <- as.factor(Gathered[[color]])
    AES <- ggplot2::aes_string(x = xvar, y = "value", fill = color)
  } else {
    AES <- ggplot2::aes_string(x = xvar, y = "value")
  }
  g <- ggplot2::ggplot(Gathered, AES) + ggplot2::facet_grid(paste0(yfacet, " ~ ", xfacet)) +
    ggplot2::geom_boxplot(outlier.size = outlier.size, lwd = boxplot.lwd, fatten = boxplot.lwd,
                          position = position_dodge2(preserve = "total"))
  if(style == "fancy"){
    g <- g + ggplot2::theme_bw() + ggplot2::scale_y_continuous(limits = ylim, breaks = seq(ylim[1], ylim[2], by = 0.1)) +
      ggplot2::ylab(ylab) + ggplot2::xlab(xlab) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
      geom_vline(xintercept = seq(1.5, length(unique(eval(parse(text = xvar), envir = Gathered))) - 0.5, 1), lwd = 0.5, colour = "black", alpha = 0.25) +
      theme(legend.position = "top")
  }
  if(print){
    print(g)
    invisible(g)
  } else {
    return(g)
  }
}

### ------------------------------- FROM MGM ------------------------------- ###

##### makeConsec: adapted from mgm:::beepday2consec
makeConsec <- function(beepvar, dayvar){
  if(!all(dayvar == round(dayvar))){stop('dayno has to be a vector of non-negative integers')}
  if(!all(beepvar == round(beepvar))){stop('beepno has to be a vector of non-negative integers')}
  if(length(beepvar) != length(dayvar)){stop('beepno has to be the same length as dayno')}
  n <- length(beepvar)
  sameday <- c(1, (dayvar[-1] == dayvar[-n])[-1])
  consec <- c(1, rep(NA, n - 1))
  counter <- 1
  for(i in 2:n){
    bdiff <- beepvar[i] - beepvar[i - 1]
    ddiff <- dayvar[i] - dayvar[i - 1]
    counter <- counter + ifelse(bdiff == 1 & ddiff == 0, 1, 2)
    consec[i] <- counter
  }
  return(consec)
}

##### lagData: adapted from mgm:::lagData
lagData <- function(data, lags, consec = NULL){
  data <- as.matrix(data)
  max_lag <- max(lags)
  lags_ext <- 1:max(lags)
  n <- nrow(data)
  p <- ncol(data)
  n_var <- nrow(data) - max(lags)
  n_lags <- length(lags_ext)
  data_response <- data
  if(!is.null(consec)){m_consec <- matrix(NA, nrow = n, ncol = n_lags)}
  l_data_lags <- list()
  lag_pos <- 1
  for(lag in lags){
    lagged_data <- matrix(NA, nrow = n, ncol = p)
    lagged_data[(lag + 1):n, ] <- data[-((n - lag + 1):n), ]
    lagged_data <- matrix(lagged_data, ncol = p, nrow = n)
    colnames(lagged_data) <- paste("V", 1:p, ".lag", lag, ".", sep = "")
    l_data_lags[[lag_pos]] <- lagged_data
    lag_pos <- lag_pos + 1
  }
  if(!is.null(consec)){
    for(lag in lags_ext){m_consec[(lag + 1):n, lag] <- consec[-((n - lag + 1):n)]}
    m_consec_check <- cbind(consec, m_consec)
    v_check <- apply(m_consec_check, 1, function(x){
      if(any(is.na(x))){
        FALSE
      } else {
        check_row <- x[1] - x[-1] == 1:length(x[-1])
        check_row_relevant <- check_row[lags_ext %in% lags]
        if(any(check_row_relevant == FALSE)){
          FALSE
        } else {
          TRUE
        }
      }
    })
  } else {
    v_check <- rep(TRUE, n)
    v_check[1:n_lags] <- FALSE
  }
  outlist <- list()
  outlist$data_response <- data_response
  outlist$l_data_lags <- l_data_lags
  outlist$included <- v_check
  return(outlist)
}

