#' Fit GVAR models with multilevel data
#'
#' Fits a graphical vector autoregressive model to data containing multiple time
#' points measured for multiple individuals.
#'
#' Uses a pseudo-mixed effects approach, wherein fixed effects are estimated and
#' random effects are approximated. See the work of Epskamp et al. (2018) for
#' more details on how these types of effects are estimated.
#'
#' @param data \code{n x k} dataframe or matrix
#' @param m Character vector or numeric vector indicating the moderator(s), if
#'   any. Can also specify \code{"all"} to make every variable serve as a
#'   moderator, or \code{0} to indicate that there are no moderators. If the
#'   length of \code{m} is \code{k - 1} or longer, then it will not be possible
#'   to have the moderators as exogenous variables. Thus, \code{exogenous} will
#'   automatically become \code{FALSE}.
#' @param selectFUN Choose a variable selection function. Can specify either
#'   \code{"varSelect"} or \code{"resample"} to use the corresponding functions.
#'   If you want to use the \code{\link{resample}} function though, then it is
#'   recommended to specify \code{selectFUN} as one of: \code{"stability",
#'   "split", "bootstrap"} in order to identify the specific method. If
#'   \code{selectFUN = "resample"}, then it is recommended to add the
#'   \code{sampMethod} argument to the call to \code{\link{mlGVAR}}.
#' @param subjectNets If \code{TRUE}, then subject-specific networks are fit for
#'   all subjects and returned in the final output. Otherwise, can specify a
#'   single value or a vector of values to represent which subjects to return
#'   individual networks for -- specifically, the SUR network. One caveat is
#'   that variable selection methods are not applied to these subject-specific
#'   networks. Further modeling could be done using the output, however.
#' @param idvar Character string to indicate which variable contains the
#'   participant identification numbers.
#' @param exogenous Logical. Indicates whether moderator variables should be
#'   treated as exogenous or not. If they are exogenous, they will not be
#'   modeled as outcomes/nodes in the network. If the number of moderators
#'   reaches \code{k - 1} or \code{k}, then \code{exogenous} will automatically
#'   be \code{FALSE}.
#' @param center Logical. Determines whether to mean-center the variables.
#' @param scale Logical. Determines whether to standardize the variables.
#' @param fixedType If logical, then any variable selection procedure specified
#'   by \code{selectFUN} will not be applied to the SUR network. Alternatively,
#'   a variable selection result, such as the output from either
#'   \code{\link{varSelect}} or \code{\link{modSelect}}, can be supplied to
#'   choose a specific constrained model in advance.
#' @param betweenType If logical, then any variable selection procedure
#'   specified by \code{selectFUN} will not be applied to the SUR network.
#'   Alternatively, a variable selection result, such as the output from either
#'   \code{\link{varSelect}} or \code{\link{modSelect}}, can be supplied to
#'   choose a specific constrained model in advance.
#' @param centerWithin Following the application of \code{center} and
#'   \code{scale}, this determines whether to center variables within individual
#'   subjects to create subject-centered values.
#' @param scaleWithin Following the application of \code{center} and
#'   \code{scale}, this determines whether to scale variables within individual
#'   subjects to create subject-standardized values.
#' @param rule Only applies to the between-subject network when a threshold is
#'   supplied. The \code{"AND"} rule will only preserve edges when both
#'   corresponding coefficients have p-values below the threshold, while the
#'   \code{"OR"} rule will preserve an edge so long as one of the two
#'   coefficients have a p-value below the supplied threshold.
#' @param threshold Logical or numeric. If \code{TRUE}, then a default value of
#'   .05 will be set. Indicates whether a threshold should be placed on the
#'   models at each iteration of the sampling. A significant choice by the
#'   researcher.
#' @param verbose Logical. Determines whether to output progress bars and
#'   messages in the console during the fitting process.
#' @param pcor See corresponding argument in the \code{\link{fitNetwork}}
#'   function
#' @param fixedArgs A named list of arguments for the variable selection
#'   function can be provided here, specifically those that are meant to be
#'   applied to the SUR network estimation.
#' @param betweenArgs A named list of arguments for the variable selection
#'   function can be provided for the between-subjects network.
#' @param bm Logical. Determines whether the same moderators are applied in the
#'   between-subjects network. By default, the value of \code{m} only applies to
#'   the SUR network. This allows one to decide whether or not to apply those
#'   moderators in the between-subject network.
#' @param beepno Character string or numeric value to indicate which variable
#'   (if any) encodes the survey number within a single day. Must be used in
#'   conjunction with \code{dayno} argument.
#' @param dayno Character string or numeric value to indiciate which variable
#'   (if any) encodes the survey number within a single day. Must be used in
#'   conjunction with \code{beepno} argument.
#' @param deleteMissing Logical. Determines whether to automatically perform
#'   listwise deletion if there are any missing values in the dataset.
#' @param ... Additional arguments.
#'
#' @return \code{mlGVAR} objects
#' @export
#' @references Epskamp, S., Waldorp, L. J., Mottus, R., & Borsboom, B. (2018).
#'   The gaussian graphical model in cross-sectional and time-series data.
#'   Multivariate Behavioral Research. 53, 453-580.
#'
#' @seealso \code{\link{mlGVARsim}, \link{lmerVAR}, \link{fitNetwork},
#'   \link{varSelect}, \link{resample}, \link{net}, \link{netInts},
#'   \link{plotNet}}
#'
#' @examples
#' \donttest{
#' fit1 <- mlGVAR(mlgvarDat, 'M')
#'
#' fit2 <- mlGVAR(mlgvarDat, 'M', bm = TRUE) # Fit the same moderator in the between-subjects network
#'
#' fit3 <- mlGVAR(mlgvarDat, 'M', selectFUN = 'varSelect')
#' }
mlGVAR <- function(data, m = NULL, selectFUN = NULL, subjectNets = FALSE, idvar = 'ID',
                   exogenous = TRUE, center = TRUE, scale = TRUE, fixedType = 'g',
                   betweenType = 'g', centerWithin = TRUE, scaleWithin = FALSE,
                   rule = 'OR', threshold = 'none', verbose = TRUE, pcor = FALSE,
                   fixedArgs = NULL, betweenArgs = NULL, bm = FALSE, beepno = NULL,
                   dayno = NULL, deleteMissing = TRUE, ...){
  t1 <- Sys.time()
  mnames <- mi <- m
  args0 <- tryCatch({list(...)}, error = function(e){list()})
  if(any(is.na(data))){
    if(deleteMissing){
      ww <- which(apply(data, 1, function(z) any(is.na(z))))
      data <- na.omit(data)
      warning(paste0(length(ww), ' rows deleted due to missingness'))
    } else {
      stop(paste0(length(ww), ' rows contain missing values'))
    }
  }
  if(!bm){
    bm <- list(moderators = NULL)
    betweenArgs <- switch(
      2 - is.null(betweenArgs), bm,
      if(!'moderators' %in% names(betweenArgs)){append(betweenArgs, bm)} else {betweenArgs})
  }
  if(is.character(threshold)){
    threshold <- sapply(match.arg(tolower(
      threshold), c('none', 'pcc', 'between', 'fixed', 'all'), several.ok = TRUE),
      function(z) switch(z, none = FALSE, all = TRUE, z), USE.NAMES = FALSE)
  }
  if(!idvar %in% colnames(data)){stop('Must supply idvar')}
  if(!is.null(m)){
    mnames <- switch(2 - is.character(m), m, colnames(data)[m])
  } else if(!is.null(selectFUN)){
    a0 <- names(args0)
    if(!'method' %in% a0){
      args0$method <- 'glmnet'
    } else if(args0$method == 'glinternet'){
      mnames <- m <- setdiff(colnames(data), idvar)
    }
    if(!'criterion' %in% a0){args0$criterion <- 'AIC'}
    if(any(c('resample', 'split') %in% selectFUN) & !'sampMethod' %in% a0){
      args0$sampMethod <- 'split'
    }
  }
  data <- data.frame(data[, -which(colnames(data) == idvar)], ID = data[, idvar])
  if(!is.null(beepno) | !is.null(dayno)){
    stopifnot(!is.null(beepno) & !is.null(dayno))
    stopifnot(length(beepno) == 1 & length(dayno) == 1)
    if(is.numeric(dayno)){dayno <- colnames(data)[dayno]}
    if(is.numeric(beepno)){beepno <- colnames(data)[beepno]}
    data0 <- data
    data <- data[, setdiff(colnames(data), c(beepno, dayno))]
  }
  vars <- setdiff(colnames(data), 'ID')
  dat <- setupVAR(data = data, idvar = 'ID', method = 'all', center = center,
                  scale = scale, vars = vars, centerWithin = centerWithin,
                  scaleWithin = scaleWithin)
  fixedDat <- dat[, vars]
  samp_ind <- as.numeric(cumsum(table(dat[, 'ID'])))
  if(!is.null(beepno) & !is.null(dayno)){
    dat <- cbind.data.frame(dat, data0[, c(beepno, dayno)])
    dat0 <- split(dat, dat$ID)
    consec <- vector('list', length(dat0))
    for(i in seq_along(dat0)){
      consec[[i]] <- which(!getConsec(data = dat0[[i]], beepno = beepno, dayno = dayno, makeAtt = FALSE))
      if(i > 1){
        consec[[i]] <- consec[[i]] + sum(sapply(1:(i - 1), function(z) nrow(dat0[[z]])))
      }
    }
    samp_ind <- union(unlist(consec), samp_ind)
  }
  attr(fixedDat, 'samp_ind') <- (1:nrow(dat))[-samp_ind]
  if(!is.null(m)){
    m <- which(vars %in% mnames)
    if(length(m) >= ncol(fixedDat) - 1){exogenous <- FALSE}
    if(!is.null(selectFUN)){args0$method <- 'glinternet'}
  }
  if(verbose){message('Estimating fixed networks')}
  fixedThresh <- ifelse(!is.character(threshold), threshold, ifelse(
    'fixed' %in% threshold, TRUE, ifelse('pcc' %in% threshold, 'PCC', FALSE)))
  fitNetArgs <- setdiff(formalArgs('fitNetwork'), '...')
  args1 <- list(data = fixedDat, moderators = m, type = fixedType, lags = 1,
                exogenous = exogenous, center = FALSE, scale = FALSE, pcor = pcor,
                rule = rule, threshold = fixedThresh, verbose = verbose)
  if(length(args0) > 0){args0 <- args0[setdiff(names(args0), names(args1))]}
  if(!is.null(fixedArgs)){
    fix1 <- intersect(names(fixedArgs), names(args1))
    if(length(fix1) > 0){args1 <- replace(args1, fix1, fixedArgs[fix1])}
    fixedArgs <- fixedArgs[setdiff(names(fixedArgs), names(args1))]
    args1 <- append(args1, fixedArgs[intersect(fitNetArgs, names(fixedArgs))])
  }
  if(is.null(selectFUN)){
    args1 <- append(args1, args0[intersect(fitNetArgs, names(args0))])
    fixedNets <- do.call(fitNetwork, args1)
  } else {
    if(isTRUE(selectFUN)){selectFUN <- 'varSelect'}
    if(length(selectFUN) == 2){betweenType <- isTRUE(selectFUN[2] == FALSE)}
    selectFUN <- match.arg(selectFUN[1], c('varSelect', 'resample', 'stability', 'bootstrap', 'split'))
    if(!selectFUN %in% c('varSelect', 'resample')){
      args0$sampMethod <- selectFUN
      selectFUN <- 'resample'
    }
    FUNargs <- setdiff(formalArgs(selectFUN), '...')
    FUN <- match.fun(selectFUN)
    args1 <- setNames(args1, gsub('moderators', 'm', names(args1)))
    args1.1 <- append(args1, args0)
    args1.1 <- args1.1[intersect(FUNargs, names(args1.1))]
    if('criterion' %in% names(args1.1)){
      if(length(args1.1$criterion) > 1){
        args1.1$criterion <- args1.1$criterion[1]}}
    if('gamma' %in% names(args1.1)){
      if(length(args1.1$gamma) > 1){args1.1$gamma <- args1.1$gamma[1]}
    }
    fixedType <- tryCatch({do.call(FUN, args1.1)}, error = function(e){TRUE})
    if(selectFUN == 'varSelect' | isTRUE(fixedType)){
      args1.2 <- append(replace(args1, c('type', 'verbose'), list(
        type = 'g', verbose = FALSE)), args0)
      if(!isTRUE(fixedType)){args1.2$type <- fixedType}
      names(args1.2)[names(args1.2) == 'm'] <- 'moderators'
      args1.2 <- args1.2[intersect(fitNetArgs, names(args1.2))]
      fixedNets <- tryCatch({do.call(fitNetwork, replace(args1.2, 'saveMods', FALSE))}, error = function(e){
        do.call(fitNetwork, replace(args1.2, c('type', 'saveMods'), list(type = 'g', saveMods = FALSE)))})
    } else {
      if('fit0' %in% names(fixedType)){fixedType$fit0 <- NULL}
      if(is(fixedType, 'list')){attr(fixedType[[2]], 'threshold') <- fixedThresh}
      args1.2 <- append(list(obj = fixedType, data = args1.1$data, fit = TRUE),
                        args0[intersect(c('select', 'thresh'), names(args0))])
      fixedNets <- tryCatch({do.call(modSelect, replace(args1.2, 'saveMods', FALSE))}, error = function(e){
        do.call(fitNetwork, setNames(args1, gsub('^m$', 'moderators', names(args1))))
      })
      if(verbose){cat('\n')}
    }
  }
  ids <- unique(dat[, 'ID'])
  meansDat <- data.frame(do.call(rbind, lapply(ids, function(z){
    z <- dat[dat[, 'ID'] == z, paste0(vars, '.m')]
    return(z[1, ])
  })))
  colnames(meansDat) <- vars
  rownames(meansDat) <- 1:nrow(meansDat)
  if(verbose){message('Estimating between-subjects network')}
  if(nrow(meansDat) <= ncol(meansDat)){
    warning('Not enough subjects to fit unregularized between-subjects network')
  }
  betThresh <- ifelse(!is.character(threshold), threshold, 'between' %in% threshold)
  args2 <- list(data = meansDat, moderators = m, type = betweenType, center = FALSE,
                scale = FALSE, exogenous = exogenous, lags = NULL, rule = rule,
                threshold = betThresh, pcor = pcor, verbose = verbose)
  if(!is.null(betweenArgs)){
    if('covariates' %in% names(betweenArgs) & !'moderators' %in% names(betweenArgs)){
      betweenArgs <- append(betweenArgs, list(moderators = NULL))
    }
    fix2 <- intersect(names(betweenArgs), names(args2))
    if(length(fix2) > 0){
      if('moderators' %in% fix2){
        if(is.null(betweenArgs$moderators) & !is.null(args2$moderators)){
          if(!'covariates' %in% names(betweenArgs) & exogenous){
            args2$data <- args2$data[, -args2$moderators]
          }
        }
      }
      args2 <- replace(args2, fix2, betweenArgs[fix2])
    }
    betweenArgs <- betweenArgs[setdiff(names(betweenArgs), names(args2))]
    args2 <- append(args2, betweenArgs[intersect(fitNetArgs, names(betweenArgs))])
  }
  if(is.null(selectFUN) | isTRUE(betweenType)){
    if(isTRUE(betweenType)){args2$type <- betweenType <- 'g'}
    args2 <- append(args2, args0[intersect(fitNetArgs, names(args0))])
    betNet <- do.call(fitNetwork, args2)
  } else {
    if(is.logical(betweenType)){
      betweenType <- 'g'
      selectFUN <- 'varSelect'
      FUNargs <- formalArgs('varSelect')
      FUN <- match.fun('varSelect')
    }
    args2 <- setNames(args2, gsub('moderators', 'm', names(args2)))
    args2.1 <- append(args2, args0)
    if(is.null(args2.1$m) & args2.1$method == 'glinternet'){args2.1$method <- 'glmnet'}
    args2.1 <- args2.1[intersect(FUNargs, names(args2.1))]
    if('criterion' %in% names(args2.1)){
      if(length(args2.1$criterion) > 1){
        args2.1$criterion <- args2.1$criterion[2]}}
    if('gamma' %in% names(args2.1)){
      if(length(args2.1$gamma) > 1){args2.1$gamma <- args2.1$gamma[2]}
    }
    betweenType <- tryCatch({do.call(FUN, args2.1)}, error = function(e){TRUE})
    if(selectFUN == 'varSelect' | isTRUE(betweenType)){
      args2.2 <- append(replace(args2, c('type', 'verbose'), list(
        type = 'g', verbose = FALSE)), args0)
      if(!isTRUE(betweenType)){args2.2$type <- betweenType}
      names(args2.2)[names(args2.2) == 'm'] <- 'moderators'
      args2.2 <- args2.2[intersect(fitNetArgs, names(args2.2))]
      betNet <- tryCatch({do.call(fitNetwork, replace(args2.2, 'saveMods', FALSE))}, error = function(e){
        do.call(fitNetwork, replace(args2.2, c('type', 'saveMods'), list(type = 'g', saveMods = FALSE)))})
    } else {
      if('fit0' %in% names(betweenType)){betweenType$fit0 <- NULL}
      if(is(betweenType, 'list')){attr(betweenType[[2]], 'threshold') <- betThresh}
      args2.2 <- append(list(obj = betweenType, data = args2.1$data, fit = TRUE),
                        args0[intersect(c('select', 'thresh'), names(args0))])
      betNet <- tryCatch({do.call(modSelect, replace(args2.2, 'saveMods', FALSE))}, error = function(e){
        do.call(fitNetwork, setNames(args2, gsub('^m$', 'moderators', names(args2))))
      })
    }
  }
  if(all(subjectNets != FALSE)){
    subjectSelect <- ifelse(all(subjectNets == 'select'), TRUE, FALSE)
    if(isTRUE(subjectNets) | isTRUE(subjectSelect)){subjectNets <- ids}
    if(any(!subjectNets %in% ids)){subjectNets <- intersect(subjectNets, ids)}
    if(verbose){
      message('Estimating subject-specific networks')
      pb <- txtProgressBar(max = length(subjectNets), style = 3)
    }
    indNets <- lapply(seq_along(subjectNets), function(i){
      dati <- dat[dat[, 'ID'] == subjectNets[i], vars]
      if(!is.null(m)){
        mi <- sapply(seq_along(m), function(j) length(unique(dati[, m[j]])) != 1)
        mi <- ifelse(all(!mi), list(NULL), list(m[mi]))[[1]]
        if(!identical(m, mi)){dati <- dati[, -setdiff(m, mi)]}
      }
      p0 <- ncol(dati)
      fiti <- nrow(dati) > (2 * p0) + 1
      if(!is.null(mi)){
        p1 <- sum(sapply(seq_along(mi), function(p) p0 - p))
        fiti <- ifelse(exogenous & length(mi) < ncol(dati) - 1,
                       nrow(dati) > p1 + (2 * p0) - (length(mi) - 1),
                       nrow(dati) > p1 + (2 * p0) + 1)
      }
      if(fiti){
        indi <- fitNetwork(data = dati, moderators = mi, threshold = fixedThresh,
                           exogenous = exogenous, center = center, type = 'g',
                           lags = 1, scale = ifelse(center, scale, FALSE),
                           saveMods = FALSE, ...)
      }
      if(verbose){setTxtProgressBar(pb, i)}
      if(fiti){return(indi)} else {return(list())}
    })
    names(indNets) <- paste0('subject', subjectNets)
    if(any(sapply(indNets, length) == 0)){
      inderrs <- subjectNets[sapply(indNets, length) == 0]
      indNets <- indNets[which(sapply(indNets, length) > 0)]
      if(length(indNets) > 0){
        message(paste0('Too few time points to estimate networks for subject',
                       ifelse(length(inderrs) == 1, ': ', 's: '),
                       paste(inderrs, collapse = ', ')))
      } else {
        subjectNets <- FALSE
        message('Too few time points to estimate subject-specific networks')
      }
    }
  }
  outcall <- list(fixedType = fixedType, betweenType = betweenType, m = mnames,
                  center = center, scale = scale, exogenous = exogenous,
                  threshold = threshold, centerWithin = centerWithin,
                  scaleWithin = scaleWithin, rule = rule)
  if(length(args0) > 0){outcall <- append(outcall, args0)}
  if(!is.null(fixedArgs) & length(fixedArgs) > 0){outcall <- append(outcall, list(fixedArgs = fixedArgs))}
  if(!is.null(betweenArgs) & length(betweenArgs) > 0){outcall <- append(outcall, list(betweenArgs = betweenArgs))}
  out <- list(call = outcall, fixedNets = fixedNets, betweenNet = betNet, ids = ids)
  if(!is.null(selectFUN)){
    out$ids <- NULL
    out$varMods <- list(fixedMods = fixedType, betweenMods = betweenType)
    out$call$fixedType <- fixedNets$call$type
    out$call$betweenType <- betNet$call$type
    out$call$selectFUN <- selectFUN
    out$ids <- ids
  }
  if(all(subjectNets != FALSE)){out$subjectNets <- indNets}
  out$netData <- list(samp_ind = samp_ind, fixedDat = fixedDat, meansDat = meansDat)
  out$dat <- dat
  attr(out, 'mlGVAR') <- TRUE
  class(out) <- c('list', 'mlGVAR')
  attr(out, 'time') <- t2 <- Sys.time() - t1
  if(verbose){cat('\nCOMPLETE:', round(t2, 2), attr(t2, 'units')); cat('\n')}
  return(out)
}

#' Main workhorse for simulating VAR and mlGVAR data
#'
#' Affords the generation of simulated data containing multiple timepoint
#' measurements for a number of subjects. Can simulate data with a single
#' moderator as well.
#'
#' Made to simulate data based on pre-specified parameters, possibly for power
#' simulations or other analyses. Output can be used to fit models with either
#' \code{\link{mlGVAR}} or \code{\link{lmerVAR}}.
#'
#' @param nTime Numeric value. The number of timepoints to simulate for each
#'   individual.
#' @param nPerson The number of subjects to create data for. Can set to \code{1}
#'   to just simulate a single graphical VAR network.
#' @param nNode The number of nodes/variables to simulate. Does not include a
#'   moderator if one is specified.
#' @param m Logical. If \code{TRUE}, then a moderator variable will be
#'   simulated. Various options also available for highly specific moderator
#'   specification: \code{"fixed", "random", "mixed1", "mixed2", "ar", "binary",
#'   "skewed", "random0", "ordinal"}.
#' @param m2 Numeric. If \code{m2 >= 1}, then this will determine the number of
#'   interaction effects between the moderator and some node in the network. If
#'   a value between 0 and 1 is provided, then this determines the probability
#'   of any given edge being moderated by the moderator.
#' @param m1 Functions similarly to \code{m2}, except that this argument refers
#'   to the number/probability of main effects of the moderator on any given
#'   node.
#' @param m0 Only relevant when \code{m = "ar"}. Determines the autoregressive
#'   coefficient in the estimated models. Defaults to .3
#' @param lag Numeric value, supposed to indicate the number of lags to simulate
#'   models parameters for. Recommended to leave at \code{1}.
#' @param thetaVar Numeric vector containing the variance associated with each
#'   node (excluding the moderator) in the contemporaneous network. If
#'   \code{NULL}, then it is assumed that the variance for each term is 1.
#' @param mu_SD Numeric vector of length 2. The first value determines the
#'   standard deviation of the means associated with the temporal data, and the
#'   second value determines the standard deviations of the means associated
#'   with the between-subjects network.
#' @param init_beta_SD Similar to \code{mu_SD} except that it applies to the
#'   coefficient estimates.
#' @param fixedMuSD Standard deviation of the random values for the means of the
#'   fixed effects.
#' @param shrink_fixed Numeric value to determine the factor by which to shrink
#'   sampled beta coefficients for fixed effects. Value between 0 and 1, where
#'   higher values are recommended.
#' @param propPos The proportion of edges with a positive sign.
#' @param m1SD Standard deviation of the moderator main effect coefficients.
#' @param m2SD Standard deviation of the moderator interaction effect
#'   coefficients.
#' @param m1_range Numeric vector of length 2. The range of values for moderator
#'   main effect coefficients.
#' @param m2_range Numeric vector of length 2. The range of values for moderator
#'   interaction effect coefficients.
#' @param shrink_deviation Numeric value to determine the factor by which to
#'   shrink contemporaneous coefficients. Value between 0 and 1, where higher
#'   values are recommended.
#' @param getM If \code{TRUE}, only the data for the moderator, the moderator
#'   main effects, and interaction effects are returned.
#' @param contemporaneous Options include \code{"wishart", "randomGGM",
#'   "fixed"}. Determines how the contemporaneous network is sampled. The former
#'   two options sample different matrices for each subject, whereas
#'   \code{"fixed"} only samples one matrix and uses it for all subject
#'   contemporaneous networks.
#' @param GGMsparsity Numeric value between 0 and 1. Determines the sparsity of
#'   sampled network matrices.
#' @param mcenter If \code{TRUE} then the moderator variable is mean-centered.
#' @param skew If \code{TRUE} then random values will be generated to represent
#'   the skewness of the node distributions. Alternatively, a numeric vector of
#'   length \code{nNode} can be provided to specify the skewness of each
#'   variable.
#' @param skewErr The skewness parameter for the \code{alpha} argument in the
#'   \code{\link[sn:rmsn]{sn::rmsn}} function.
#' @param ordinal Logical. Determines whether to sample ordinal variables. If a
#'   numeric value is provided, then this will automatically be assigned to the
#'   \code{nLevels} argument.
#' @param nLevels Number of levels for the ordinal variables. Only relevant if
#'   \code{ordinal} is not \code{FALSE}.
#' @param ordWithin If \code{TRUE}, then variables will ordinalized within
#'   subjects, rather than across subjects.
#' @param minOrd The minimum number of unique values allowed for each variable.
#' @param thresholds List of length \code{k}, where each element is a numeric
#'   vector of length \code{(nLevels - 1)} containing the splitpoints for
#'   grouping each variable into ordered categories.
#' @param mseed Numeric value for the seed to be set when
#' @param onlyNets If \code{TRUE} then only the network models are returned,
#'   without the data. Could be used to create random models and then simulate
#'   data by another method.
#' @param modType Determines the type of moderation to employ, such as
#'   \code{"none", "full", "partial"}. See \code{\link{simNet}}
#'
#' @return Simulated mlGVAR or VAR data.
#' @export
#'
#' @seealso \code{\link{mlGVAR}, \link{lmerVAR}, \link{simNet}, \link{plotNet},
#'   \link{net}, \link{netInts}}
#'
#' @examples
#' set.seed(1)
#' x <- mlGVARsim(nTime = 50, nPerson = 10, nNode = 3, m = TRUE)
mlGVARsim <- function(nTime = 50, nPerson = 10, nNode = 3, m = NULL, m2 = .25, m1 = .7,
                      m0 = 1, lag = 1, thetaVar = NULL, mu_SD = NULL, init_beta_SD = NULL,
                      fixedMuSD = 1, shrink_fixed = 0.9, propPos = .5, m1SD = .1, m2SD = .1,
                      m1_range = NULL, m2_range = NULL, shrink_deviation = 0.9, getM = FALSE,
                      contemporaneous = "wishart", GGMsparsity = .5, mcenter = TRUE, skew = FALSE,
                      skewErr = FALSE, ordinal = FALSE, nLevels = 5, ordWithin = TRUE, minOrd = 3,
                      thresholds = NULL, mseed = NULL, onlyNets = FALSE, modType = 'none'){
  modType <- match.arg(tolower(modType), c('none', 'full', 'partial', 'full2', 'partial2', 'zero'))
  if(identical(m, FALSE)){m <- NULL}
  if(minOrd < 2){minOrd <- 2}
  if(is.numeric(ordinal)){nLevels <- ordinal; ordinal <- TRUE}
  if(is.null(thetaVar)){thetaVar <- rep(1, nNode)}
  if(is.null(mu_SD)){mu_SD <- c(1, 1)}
  if(is.null(init_beta_SD)){init_beta_SD <- c(.1, 1)}
  if(is.null(m1_range)){m1_range <- c(.1, .4)}
  if(is.null(m2_range)){m2_range <- c(.1, .3)}
  if(length(nTime) == 1){nTime <- rep(nTime, nPerson)}
  DF_theta <- nNode * 2
  nTemporal <- nNode^2 * lag
  contemporaneous <- match.arg(contemporaneous, c("wishart", "randomGGM", "fixed"))
  Omega_mu <- simCor(nNode, GGMsparsity)
  Omega_Beta <- simCor(nTemporal, GGMsparsity)
  Omega <- rbind(cbind(Omega_mu, matrix(0, nNode, nTemporal)),
                 cbind(matrix(0, nTemporal, nNode), Omega_Beta))
  SD <- runif(nNode + nTemporal, c(rep(mu_SD[1], nNode), rep(init_beta_SD[1], nNode)),
              c(rep(mu_SD[2], nNode), rep(init_beta_SD[2], nNode)))
  Omega <- diag(SD) %*% Omega %*% diag(SD)
  if(contemporaneous == "wishart"){
    Theta_fixed <- simCor(nNode, GGMsparsity)
    Theta_fixed <- diag(sqrt(thetaVar)) %*% Theta_fixed %*% diag(sqrt(thetaVar))
    Theta <- rWishart(nPerson, DF_theta, Theta_fixed/DF_theta)
  } else if(contemporaneous == "randomGGM"){
    Theta <- lapply(1:nPerson, function(x) simCor(nNode, GGMsparsity))
    Theta <- do.call(abind::abind, c(Theta, along = 3))
    Theta_fixed <- apply(Theta, 1:2, mean)
  } else {
    Theta_fixed <- simCor(nNode, GGMsparsity)
    Theta <- lapply(1:nPerson, function(x) Theta_fixed)
    Theta <- do.call(abind::abind, c(Theta, along = 3))
  }
  if(!is.null(m)){
    if(isTRUE(m)){m <- "random"}
    m <- match.arg(tolower(m), c(
      "fixed", "random", "mixed1", "mixed2", "ar", "binary",
      "skewed", "random0", "ordinal"), several.ok = TRUE)
    if(all(c("fixed", "random") %in% m)){m <- m[-which(m == "fixed")]}
    if("ar" %in% m){
      if(m0 >= 1){m0 <- .3}
      mm <- lapply(seq_len(nPerson), function(z){
        as.numeric(arima.sim(n = nTime[z] + 100, model = list(ar = m0)))})
    } else if("binary" %in% m){
      mm <- lapply(seq_len(nPerson), function(z) rbinom(nTime[z] + 100, 1, .5))
    } else if("skewed" %in% m){
      mm <- lapply(seq_len(nPerson), function(z) sn::rsn(nTime[z] + 100, 0, 1, m0))
    } else {
      mm <- lapply(seq_len(nPerson), function(z) rnorm(nTime[z] + 100, 0, m0))
    }
    if(!is.null(mseed)){set.seed(mseed)}
    if(m2 >= 0 & m2 < 1){
      m20 <- diag(0, nNode)
      if(modType != 'zero'){
        wcond <- function(m20, modType){
          all(m20 == 0) | ifelse(grepl('full', modType), any(diag(m20) != 0), FALSE)
        }
        while(wcond(m20, modType)){
          m20 <- matrix(sample(c(0, 1), nNode^2, TRUE, prob = c(1 - m2, m2)), nNode, nNode)
        }
      }
    } else if(m2 >= 1){
      m20 <- numeric(nNode^2)
      m20[sample(1:(nNode^2), ifelse(m2 > nNode^2, nNode^2, round(m2)))] <- 1
      if(modType == 'zero'){
        m20 <- diag(0, nNode)
      } else if(grepl('full', modType)){
        m20 <- matrix(m20, nNode, nNode)
        while(any(diag(m20) != 0)){
          if(m2 >= (nNode^2) - nNode){
            m20 <- matrix(1, nNode, nNode)
            diag(m20) <- 0
          } else {
            m20 <- numeric(nNode^2)
            m20[sample(1:(nNode^2), round(m2))] <- 1
            m20 <- matrix(m20, nNode, nNode)
          }
        }
      }
    }
    if(isTRUE(m1) | is.null(m1)){m1 <- .7}
    if(m1 >= 0 & m1 < 1){
      m10 <- sample(c(0, 1), nNode, TRUE, prob = c(1 - m1, m1))
    } else if(m1 >= 1){
      m10 <- numeric(nNode)
      m10[sample(1:nNode, ifelse(m1 > nNode, nNode, round(m1)))] <- 1
    }
    if(grepl('2', modType) & !all(m20 == 0)){
      mm12 <- apply(m20, 1, function(z) any(z != 0))
      while(ifelse(grepl('full', modType), !all(m10[mm12] == 0), any(m10[mm12] == 0))){
        if(m1 >= 0 & m1 < 1){
          m10 <- sample(0:1, nNode, TRUE, prob = c(1 - m1, m1))
        } else {
          m10 <- numeric(nNode)
          m10[sample(1:nNode, ifelse(m1 > nNode, nNode, round(m1)))] <- 1
        }
      }
    }
    m1 <- m10 * sample(c(-1, 1), nNode, TRUE, prob = c(1 - propPos, propPos))
    m2 <- m20 * sample(c(-1, 1), nNode^2, TRUE, prob = c(1 - propPos, propPos))
    if(any(c("fixed", "mixed1") %in% m)){
      m1 <- m1 * runif(nNode, min(m1_range), max(m1_range))
      mmb1 <- rep(list(m1), nPerson)
    } else if("random0" %in% m){
      mmb1 <- lapply(1:nPerson, function(z){
        m10 <- rnorm(nNode, 0, m1SD)
        while(any(abs(m10) < min(m1_range)) | any(abs(m10) > max(m1_range))){
          m10 <- rnorm(nNode, 0, m1SD)
        }
        m10 <- m1 * m10
        return(m10)
      })
    } else {
      mmb1 <- lapply(1:nPerson, function(z) m1 * runif(nNode, min(m1_range), max(m1_range)))
    }
    if(any(c("fixed", "mixed2") %in% m)){
      m2 <- matrix(m2 * runif(nNode^2, min(m2_range), max(m2_range)), ncol = nNode)
      mmb2 <- rep(list(m2), nPerson)
    } else if("random0" %in% m){
      mmb2 <- lapply(1:nPerson, function(z){
        m20 <- rnorm(nNode^2, 0, m2SD)
        while(any(abs(m20) < min(m2_range)) | any(abs(m20) > max(m2_range))){
          m20 <- rnorm(nNode^2, 0, m2SD)
        }
        m20 <- c(m2) * m20
        return(matrix(m20, ncol = nNode))
      })
    } else {
      mmb2 <- lapply(1:nPerson, function(z){
        matrix(c(m2) * runif(nNode^2, min(m2_range), max(m2_range)), ncol = nNode)})
    }
    if(getM){return(list(m = mm[[1]], mb1 = mmb1[[1]], mb2 = mmb2[[1]]))}
  } else {
    mm <- mmb1 <- mmb2 <- NULL
  }
  mu_fixed <- rnorm(nNode, 0, fixedMuSD)
  beta_fixed <- rnorm(nTemporal, 0)
  beta_fixed[order(abs(beta_fixed))[1:round(nTemporal * GGMsparsity)]] <- 0
  mat <- matrix(0, nNode, nNode * lag)
  diag(mat) <- 1
  t1 <- 0
  beta_fixed[c(mat) == 1] <- runif(sum(c(mat) == 1), 0, 1)
  if(!is.null(m) & !modType %in% c('none', 'zero')){
    beta_fixed2 <- matrix(beta_fixed, nNode, nNode)
    while(ifelse(grepl('full', modType), !all(beta_fixed2[m2 != 0] == 0), any(beta_fixed2[m2 != 0] == 0))){
      beta_fixed <- matrix(rnorm(nTemporal, 0), nNode, nNode)
      allzeros <- round(nTemporal * GGMsparsity)
      if(grepl('full', modType)){
        beta_fixed[m2 != 0] <- 0
        allzeros <- allzeros - sum(beta_fixed == 0)
        bfix <- beta_fixed[beta_fixed != 0]
        beta_fixed[beta_fixed != 0][order(abs(bfix))[1:allzeros]] <- 0
      } else {
        bfix <- beta_fixed[m2 == 0]
        beta_fixed[m2 == 0][order(abs(bfix))[1:allzeros]] <- 0
      }
      mat <- switch(2 - (lag == 1), diag(nNode), cbind(diag(nNode), matrix(0, nNode, nNode * (lag - 1))))
      beta_fixed[c(mat) == 1] <- runif(sum(c(mat) == 1), 0, 1)
      beta_fixed2 <- matrix(beta_fixed, nNode, nNode)
    }
  }
  if(lag > 0){
    repeat{
      t1 <- t1 + 1
      if(all(skew == FALSE)){
        Pars <- mvtnorm::rmvnorm(nPerson, c(mu_fixed, beta_fixed), sigma = Omega)
      } else {
        if(is.numeric(skew)){skew <- rep(skew, nNode + nNode^2)}
        if(is.character(skew)){if(all(skew == 'random')){skew <- runif(nNode + nNode^2, -10000, 10000)}}
        if(isTRUE(skew) | length(skew) != (nNode + nNode^2)){skew <- rep(3, nNode + nNode^2)}
        Omega <- as.matrix(Matrix::forceSymmetric(Omega))
        Pars <- sn::rmsn(nPerson, c(mu_fixed, beta_fixed), Omega, skew)
      }
      Mus <- Pars[, 1:nNode]
      Betas <- array(c(t(Pars[, -(1:nNode)])), c(nNode, nNode * lag, nPerson))
      if(lag > 1){
        under <- cbind(diag(nNode * (lag - 1)), matrix(0, nNode * (lag - 1), nNode))
        ev <- sapply(seq_len(nPerson), function(i){
          mat <- rbind(Betas[, , i], under)
          eigen(mat)$values
        })
      } else {
        ev <- sapply(seq_len(nPerson), function(i){eigen(Betas[, , i])$values})
      }
      allEV <- c(ev)
      if(all(Re(ev)^2 + Im(ev)^2 < 1)){
        if(nPerson == 1){
          Mus <- matrix(mu_fixed, ncol = nNode, nrow = 1)
          Betas <- array(matrix(beta_fixed, nNode, nNode), c(nNode, nNode, 1))
          Theta <- array(Theta_fixed, c(nNode, nNode, 1))
        }
        DataList <- lapply(1:nPerson, function(p){
          parms <- lapply(seq_len(lag), function(l) array(c(Betas[, , p]), c(nNode, nNode, lag))[, , l])
          if(lag > 0){
            if(onlyNets){
              res <- list(parms = parms, means = Mus[p, ], lags = seq_len(lag),
                          Nt = nTime[p], init = Mus[p, ], burnin = 100,
                          residuals = Theta[, , p], m = mm[[p]], mb1 = mmb1[[p]],
                          mb2 = mmb2[[p]], mcenter = mcenter, skewErr = skewErr)
              return(res)
            }
            res <- trevSimulateVAR(parms, means = Mus[p, ], lags = seq_len(lag),
                                   Nt = nTime[p], init = Mus[p, ], burnin = 100,
                                   residuals = Theta[, , p], m = mm[[p]],
                                   mb1 = mmb1[[p]], mb2 = mmb2[[p]],
                                   mcenter = mcenter, skewErr = skewErr)
          } else {
            res <- mvtnorm::rmvnorm(nTime[p], Mus[p, ], Theta[, , p])
          }
          colnames(res) <- paste0("V", 1:nNode)
          if(ordinal & ordWithin){
            for(vv in 1:ncol(res)){
              tick <- ord <- 0
              while(length(unique(ord)) < minOrd){
                thresh <- switch(2 - is.null(thresholds), rnorm(nLevels - 1), thresholds[[vv]])
                ord <- as.numeric(cut(res[, vv], sort(c(-Inf, thresh, Inf))))
                tick <- tick + 1
                if(tick == 10 & !is.null(thresholds)){thresholds <- NULL}
                if(tick == 20){break}
              }
              res[, vv] <- ord
            }
          }
          res$ID <- p
          res
        })
        if(onlyNets){return(DataList[[1]])}
        Data <- do.call(rbind, DataList)
        if(!any(abs(Data[, 1:nNode]) > 100)){break}
      }
      beta_fixed <- beta_fixed * shrink_fixed
      D <- diag(sqrt(diag(Omega)))
      D[-(1:nNode), -(1:nNode)] <- shrink_deviation * D[-(1:nNode), -(1:nNode)]
      Omega <- D %*% cov2cor(Omega) %*% D
    }
  } else {
    Pars <- mvtnorm::rmvnorm(nPerson, mu_fixed, sigma = Omega)
    Mus <- Pars[, 1:nNode]
    Betas <- array(dim = c(0, 0, nPerson))
    DataList <- lapply(1:nPerson, function(p){
      res <- as.data.frame(mvtnorm::rmvnorm(nTime[p], Mus[p, ], Theta[, , p]))
      colnames(res) <- paste0("V", 1:nNode)
      res$ID <- p
      res
    })
    Data <- do.call(rbind, DataList)
  }
  model <- list(mu = trevModelArray(mean = mu_fixed, SD = mu_SD, subject = lapply(1:nrow(Mus), function(i) Mus[i, ])),
                Beta = trevModelArray(mean = array(beta_fixed, c(nNode, nNode, lag)),
                                      SD = array(sqrt(diag(Omega[-(1:nNode), -(1:nNode)])), c(nNode, nNode, lag)),
                                      subject = lapply(1:nPerson, function(p) array(Betas[, , p], c(nNode, nNode, lag)))),
                Omega_mu = trevModelCov(cov = trevModelArray(mean = Omega[1:nNode, 1:nNode])),
                Theta = trevModelCov(cov = trevModelArray(mean = Theta_fixed, subject = lapply(1:nPerson, function(p) Theta[, , p]))),
                Omega = trevModelCov(cov = trevModelArray(mean = Omega)))
  kappa <- corpcor::pseudoinverse(Theta_fixed)
  beta <- matrix(beta_fixed, nNode, nNode)
  # FUNCTION FOR THIS TO WORK
  trevPDC <- function(beta, kappa){
    if(ncol(beta) == nrow(beta) + 1){beta <- beta[, -1, drop = FALSE]}
    sigma <- solve(kappa)
    t(beta/sqrt(diag(sigma) %o% diag(kappa) + beta^2))
  }
  mod2 <- list(fixedKappa = kappa, fixedPCC = pcor2(Theta_fixed), fixedBeta = beta,
               fixedPDC = trevPDC(beta, kappa), between = pcor2(Omega[1:nNode, 1:nNode]))
  if(ordinal & !ordWithin){
    for(vv in 1:(ncol(Data) - 1)){
      tick <- ord <- 0
      while(length(unique(ord)) < minOrd){
        thresh <- switch(2 - is.null(thresholds), rnorm(nLevels - 1), thresholds[[vv]])
        ord <- as.numeric(cut(Data[, vv], sort(c(-Inf, thresh, Inf))))
        tick <- tick + 1
        if(tick == 10 & !is.null(thresholds)){thresholds <- NULL}
        if(tick == 20){break}
      }
      Data[, vv] <- ord
    }
  }
  Results <- list(data = Data, vars = paste0("V", 1:nNode), idvar = "ID", lag = lag, model = model)
  Results <- append(list(data = Data), append(mod2, Results[-1]))
  if(!is.null(m)){
    if('ordinal' %in% m){
      mm <- lapply(mm, function(z){
        ord <- c()
        while(length(unique(ord)) < minOrd){
          ord <- as.numeric(cut(z, sort(c(-Inf, rnorm(nLevels - 1), Inf))))
        }
        return(ord)
      })
    }
    mm <- lapply(mm, function(z) z[-(1:100)])
    dat2 <- data.frame(Data[, 1:nNode], M = unlist(mm), ID = Data[, "ID"])
    mb1 <- if(any(c("fixed", "mixed1") %in% m)){mmb1[[1]]} else {Reduce("+", mmb1)/nPerson}
    mb2 <- if(any(c("fixed", "mixed2") %in% m)){mmb2[[1]]} else {Reduce("+", mmb2)/nPerson}
    Results$data <- dat2 # NEW
    Results$interactions <- list(m = m, subjects = list(mb1 = mmb1, mb2 = mmb2), mb1 = mb1, mb2 = mb2) # NEW
    #Results$mm <- list(data = dat2, m = m, subjects = list(mb1 = mmb1, mb2 = mmb2), mb1 = mb1, mb2 = mb2)
    #if(nPerson == 1){Results$mm$data <- Results$mm$data[, -ncol(Results$mm$data)]}
  }
  if(nPerson == 1){
    Results$data <- Results$data[, -ncol(Results$data)]
    Results$between <- NULL
    #attr(Results, "simMLgvar") <- TRUE
    #class(Results) <- c('list', 'simMLgvar')
    attr(Results, 'GVARsim') <- TRUE # NEW
    class(Results) <- c('list', 'GVARsim') # NEW
  } else {
    #class(Results) <- attr(Results, "mlVARsim") <- "mlVARsim"
    class(Results) <- attr(Results, "mlGVARsim") <- "mlGVARsim" # NEW
  }
  return(Results)
}

##### simCor: simulate correlation matrix
simCor <- function(x, sparsity = .5, maxiter = 100, ...){
  #ow <- getOption("warn"); options(warn = 2)
  try <- 0
  while(try < maxiter){
    out <- tryCatch({
      cov2cor(solve(diag(x) - simPcor(x, sparsity, ...)))},
      error = function(e){TRUE})
    if(!isTRUE(out)){break}
    try <- try + 1
  }
  #options(warn = ow)
  if(try == maxiter){stop("Estimate failed to converge")}
  return(out)
}

##### trevModelArray: straight from mlVAR
trevModelArray <- function(mean, subject, SE, P, lower, upper, SD){
  Results <- list()
  if(missing(subject)){subject <- NULL}
  if(is.null(subject) && missing(mean)){stop("Either 'subject' or 'mean' must be used")}
  if(missing(mean) && !is.null(subject)){
    N <- length(subject)
    mean <- Reduce("+", subject)/N
  }
  dim <- switch(2 - is.null(dim(mean)), length(mean), dim(mean))
  if(is.null(dim)){dim <- length(mean)}
  if(missing(SD)){
    if(!is.null(subject)){
      N <- length(subject)
      SD <- Reduce("+", lapply(subject, function(x) (x - mean)^2))/(N - 1)
    } else {
      SD <- array(NA, dim)
    }
  }
  if(missing(P)){
    if(!missing(SE)){
      P <- 2 * (1 - pnorm(abs(mean/SE)))
    } else {
      P <- array(NA, dim)
    }
  }
  if(missing(SE)){SE <- array(NA, dim)}
  if(missing(lower)){
    if(!missing(SE)){
      lower <- mean - 1.959964 * SE
    } else {
      lower <- array(NA, dim)
    }
  }
  if(missing(upper)){
    if(!missing(SE)){
      upper <- mean + 1.959964 * SE
    } else {
      upper <- array(NA, dim)
    }
  }
  if(is.null(subject)){subject <- NULL}
  Results[["mean"]] <- mean
  Results[["SD"]] <- SD
  Results[["lower"]] <- lower
  Results[["upper"]] <- upper
  Results[["SE"]] <- SE
  Results[["P"]] <- P
  Results[["subject"]] <- subject
  class(Results) <- c("mlVARarray", "list")
  return(Results)
}

##### trevModelCov: straight from mlVAR
trevModelCov <- function(cov, cor, prec, pcor){
  if(missing(cov)){stop("'cov' can not be missing")}
  if(!missing(cor) && !is(cor, "mlVARarray")){stop("'cor' must be missing or an object of class 'mlVARarray'")}
  if(!missing(prec) && !is(prec, "mlVARarray")){stop("'prec' must be missing or an object of class 'mlVARarray'")}
  if(!missing(pcor) && !is(pcor, "mlVARarray")){stop("'pcor' must be missing or an object of class 'mlVARarray'")}
  cov2corNA <- function(x){
    if(any(is.na(x))){
      x[] <- NA
      return(x)
    } else {
      return(cov2cor(x))
    }
  }
  solveNA <- function(x){
    if(any(is.na(x))){
      x[] <- NA
      return(x)
    } else {
      return(corpcor::pseudoinverse(x))
    }
  }
  cor2pcorNA <- function(x){
    if(any(is.na(x))){
      x[] <- NA
      return(x)
    } else {
      return(corpcor::cor2pcor(x))
    }
  }
  if(missing(cor)){
    cor <- trevModelArray(mean = cov2corNA(cov[["mean"]]),
                          subject = lapply(cov[["subject"]], cov2corNA))
  }
  if(missing(prec)){
    prec <- trevModelArray(mean = solveNA(cov[["mean"]]),
                           subject = lapply(cov[["subject"]], solveNA))
  }
  if(missing(pcor)){
    pcor <- trevModelArray(mean = cor2pcorNA(cov[["mean"]]),
                           subject = lapply(cov[["subject"]], cor2pcorNA))
  }
  Results <- list(cov = cov, cor = cor, prec = prec, pcor = pcor)
  class(Results) <- c("mlVarCov", "list")
  return(Results)
}
