#' Fit cross-sectional and idiographic moderated network models
#'
#' The main function that ties everything together for both cross-sectional and
#' idiographic (temporal) network models, moderated or otherwise.
#'
#' For GGMs, nodewise estimation is utilized to fit models to each node, and
#' then aggregate results into the final network. For temporal networks that
#' represent data for a single subject, SUR estimation based on feasible
#' generalized least squares (FGLS) is used. Also incorporates the variable
#' selection functions to integrate model selection and estimation. Nodewise
#' estimation is used for all GGMs, and SUR estimation is used for temporal
#' networks. See \code{systemfit} package for more information on the latter,
#' particularly via the \code{\link[systemfit:systemfit]{systemfit::systemfit}}
#' function.
#'
#' @param data \code{n x k} dataframe or matrix.
#' @param moderators Numeric or character vector indicating which variables (if
#'   any) to use as moderators.
#' @param type Primarily used to supply a variable selection object, such as
#'   those created with \code{\link{varSelect}} or \code{\link{modSelect}}, or
#'   to indicate that a variable selection method should be employed by setting
#'   the value to \code{"varSelect"}. Currently doesn't support setting the
#'   value to \code{"resample"}, although this will be implemented in the
#'   future. Alternatively, this can be used to specify the type of variable for
#'   each node. In this case it should be either a single value --
#'   \code{"gaussian"} or \code{"binomial"} -- or can be a vector of length
#'   \code{k} to specify which of those two types apply to each variable. These
#'   dictate which family to use for the call to
#'   \code{\link[stats:glm]{stats::glm}}. Cannot use binomial models for SUR
#'   networks.
#' @param lags Logical or numeric, to indicate whether to fit a SUR model or
#'   not. Set to \code{TRUE} or 1 for a SUR model fit to temporal data for a
#'   single subject.
#' @param seed Only useful if \code{type = "varSelect"}, and if the
#'   \code{varSeed} argument is not specified in the \code{...}
#' @param folds Can be used to specify the number of folds in cross-validation
#'   when \code{type = "varSelect"} and \code{criterion = "CV"}. Overwritten if
#'   \code{nfolds} argument is provided.
#' @param gamma Only useful if \code{type = "varSelect"} and the criterion is
#'   set to \code{"EBIC"}. This is the hyperparameter for the calculation of
#'   EBIC.
#' @param which.lam Only useful if \code{criterion = "CV"}, or if a variable
#'   selection object based on cross-validation is supplied for \code{type}.
#'   Options include \code{"min"}, which uses the lambda value that minimizes
#'   the objective function, or \code{"1se"} which uses the lambda value at 1
#'   standard error above the value that minimizes the objective function.
#' @param rule Only applies to GGMs (including between-subjects networks) when a
#'   threshold is supplied. The \code{"AND"} rule will only preserve edges when
#'   both corresponding coefficients have p-values below the threshold, while
#'   the \code{"OR"} rule will preserve an edge so long as one of the two
#'   coefficients have a p-value below the supplied threshold.
#' @param threshold Determines whether to employ a p-value threshold on the
#'   model. If \code{TRUE} then this defaults to .05. Not recommended, as
#'   thresholds can be applied post-hoc through the plotting functions, or via
#'   the \code{\link{net}} and \code{\link{netInts}} functions. Recommended to
#'   leave as \code{FALSE}.
#' @param scale Determines whether to standardize all variables or not.
#' @param std Only applies to SUR networks. Logical. Provides input to the
#'   \code{method} argument of the
#'   \code{\link[systemfit:systemfit]{systemfit::systemfit}} function. If
#'   \code{TRUE}, then the \code{method} will be \code{"SUR"}. If \code{FALSE},
#'   then the \code{method} will be \code{"OLS"}. These two methods only differ
#'   when constraints are applied. When a saturated model is fit, both methods
#'   produce the same results.
#' @param center Determines whether to mean-center variables or not.
#' @param covariates Either a numeric value or character string -- this could
#'   also be a vector -- to indicate which variables (if any) should be treated
#'   as covariates in the model.
#' @param verbose Logical. Determines whether to return information about the
#'   progress of the model fitting -- especially when variable selection is
#'   employed -- as well as prints the amount of time it takes to fit the model
#'   to the console.
#' @param exogenous Logical. Indicates whether moderator variables should be
#'   treated as exogenous or not. If they are exogenous, they will not be
#'   modeled as outcomes/nodes in the network. If the number of moderators
#'   reaches \code{k - 1} or \code{k}, then \code{exogenous} will automatically
#'   be \code{FALSE}.
#' @param mval Numeric value to set the moderator variable to when computing
#'   model coefficients. Useful to create conditional networks -- i.e., those
#'   whose values are conditioned on specific values of the moderator. Excellent
#'   when the moderator is a categorical variable, or when it's desired to have
#'   model estimates at +/- 1 SD around the mean of the moderator. These values
#'   must be supplied explicitly. Can only specify a single value for a given
#'   model.
#' @param residMat Character string indicating which type of residual covariance
#'   matrix to compute for SUR models. Options include \code{"res", "dfres",
#'   "sigma"}. \code{"sigma"} uses the residual covariance matrix as computed by
#'   the \code{systemfits} package. \code{"res"} and \code{"dfres"} compute the
#'   matrix based directly on the residual values. \code{"dfres"} is the sample
#'   estimator that uses \code{N - 1} in the denominator, while \code{"res"}
#'   just uses \code{N}. Input for \code{\link{SURnet}} function.
#' @param medges Only relevant when \code{lags = 1} and \code{exogenous =
#'   FALSE}. Determines the linetype of moderated edges (corresponds to the lty
#'   argument of \code{plot()}).
#' @param pcor Logical. Determines whether to operationalize the adjacency
#'   matrix as the partial correlation matrix of the data, or to use nodewise
#'   estimation. Only relevant for unmoderated networks.
#' @param maxiter See argument of \code{SURfit()} function.
#' @param getLL Logical. Determines whether to return log-likelihood statistics
#'   with model results. Recommended to keep \code{TRUE}.
#' @param saveMods Logical. Determines whether to save the \code{fitobj} element
#'   of the output, which contains the nodewise models, or the SUR model output
#'   of \code{\link[systemfit:systemfit]{systemfit::systemfit}}.
#' @param binarize Logical. Determines whether to convert the output to a
#'   binary, unweighted network. Only relevant for GGMs.
#' @param fitCoefs Determines whether to use the \code{\link{getFitCIs}}
#'   function on the output. Not recommended to use. The downside is that this
#'   will overwrite the \code{fitobj} element of the output which contains the
#'   actual models. Better to leave this as \code{FALSE}, and then use the
#'   \code{\link{getFitCIs}} function on the object separately.
#' @param detrend Logical. Determines whether to remove linear trends from time
#'   series variables. Only applies to temporal networks.
#' @param beepno Character string or numeric value to indicate which variable
#'   (if any) encodes the survey number within a single day. Must be used in
#'   conjunction with \code{dayno} argument. Only relevant to temporal data.
#' @param dayno Character string or numeric value to indicate which variable (if
#'   any) encodes the survey number within a single day. Must be used in
#'   conjunction with \code{beepno} argument. Only relevant to temporal data.
#' @param ... Additional arguments.
#'
#' @return A ggm or SUR network
#' @export
#'
#' @examples
#' fit1 <- fitNetwork(ggmDat)
#'
#' \donttest{
#' fit2 <- fitNetwork(ggmDat, 'M', type = 'varSelect', criterion = 'BIC')
#' }
#'
#' fit3 <- fitNetwork(gvarDat, 'M', lags = 1)
fitNetwork <- function(data, moderators = NULL, type = "gaussian", lags = NULL,
                       seed = NULL, folds = 10, gamma = 0.5, which.lam = 'min',
                       rule = "OR", threshold = FALSE, scale = FALSE, std = TRUE,
                       center = TRUE, covariates = NULL, verbose = FALSE, exogenous = TRUE,
                       mval = NULL, residMat = "sigma", medges = 1, pcor = FALSE, maxiter = 100,
                       getLL = TRUE, saveMods = TRUE, binarize = FALSE, fitCoefs = FALSE,
                       detrend = FALSE, beepno = NULL, dayno = NULL, ...){
  # START
  t1 <- Sys.time()

  # Alternative way of specifying that there is a lag
  if(!identical(detrend, FALSE) | (!is.null(beepno) & !is.null(dayno))){lags <- 1}

  # Check if moderators are specified by name -- NEW
  if(isTRUE(is.character(moderators)) & !identical(moderators, 'all')){
    moderators <- which(colnames(data) %in% moderators)
  }

  # Check if covariates are specified by name -- NEW
  if(isTRUE(is.character(covariates))){
    covariates <- which(colnames(data) %in% covariates)
  }

  # Check for missing values
  if(any(is.na(data))){
    ww <- which(apply(data, 1, function(z) any(is.na(z))))
    if(is.null(lags) | identical(lags, 0)){
      data <- data[-ww, ]
      warning(paste0(length(ww), ' cases deleted due to missingness'))
    } else {
      stop(paste0(length(ww), ' rows contain missing values'))
    }
  }

  # Collect additional arguments
  args <- tryCatch({list(...)}, error = function(e){list()})

  # Variable selection -- REMOVE SECOND PART
  if(identical(type, 'varSelect')){
    vargs <- list(data = data, m = moderators, lags = lags, exogenous = exogenous,
                  center = center, scale = scale, gamma = gamma, verbose = verbose)
    vargs$method <- ifelse(!is.null(moderators), 'glinternet',
                           ifelse('method' %in% names(args), args$method, 'glmnet'))
    otherargs <- c('criterion', 'nfolds', 'varSeed', 'useSE', 'nlam')
    if(any(otherargs %in% names(args))){
      if('criterion' %in% names(args)){ # This conditional is new
        if(toupper(args$criterion) == 'CV'){
          if(!'nfolds' %in% names(args)){args$nfolds <- folds}
        }
      }
      vargs <- append(vargs, args[intersect(names(args), otherargs)])
    }
    if(!is.null(seed)){vargs$varSeed <- seed}
    type <- do.call(varSelect, vargs)
  } else if(identical(type, 1)){
    type <- 'g'
    lags <- 1
  }

  # Ensure that if lags == 0, then NULL
  if(!is.null(lags)){if(all(lags == 0)){lags <- NULL}}

  # Obviously related to bootNet; ensures that dataset will be data.frame
  if(!is.null(lags)){
    if(lags != FALSE){lags <- 1}
    if("samp_ind" %in% names(attributes(data))){samp_ind <- attr(data, "samp_ind")}
    data <- data.frame(data)
    if(exists("samp_ind", inherits = FALSE)){attr(data, "samp_ind") <- samp_ind}
  } else {
    data <- data.frame(data)
  }

  # Checking for weird moderator specification
  if(!is.null(moderators)){
    if(all(moderators == 0)){moderators <- NULL}
    if(all(moderators == "all")){moderators <- 1:ncol(data)}
  }

  # Covariates must be specified as a list? Or numeric?
  if(ifelse(!is.null(covariates), ifelse(!is(covariates, 'list'), TRUE, FALSE), FALSE)){
    if(length(covariates) >= ncol(data) - 1){stop("Must have at least 2 outcome variables")}
  }

  # adjust which.lam
  if(grepl("min", which.lam)){which.lam <- "min"} else {which.lam <- "1se"}

  ### SETUP -- create output object
  output <- list()
  output$call <- list(type = type, moderators = moderators, mval = mval, lags = lags,
                      which.lam = which.lam, rule = rule, threshold = threshold,
                      center = center, scale = scale)

  # Add extra arguments
  if(length(args) != 0){output$call <- append(output$call, args)}

  # Removing mval if it's null
  if(is.null(mval)){output$call$mval <- NULL}

  # Make adjustments to TYPE if variable selection was used; otherwise, setup gaussian and binomial
  if(class(type) == "list"){
    output$call$type <- lapply(type, function(z) unlist(z[[ifelse(which.lam == "min", 1, 2)]]))
    if(attr(type, "criterion") != "CV"){output$call$which.lam <- attr(type, "criterion")}
    if(attr(type, "criterion") == "EBIC"){output$call$gamma <- attr(type, "gamma")}
  } else {
    output$call$which.lam <- NULL
    if(length(type) == 1){
      type <- rep("gaussian", ncol(data))
      bb <- unname(which(apply(data, 2, function(z) dim(table(z)) <= 2)))
      if(length(bb) > 0){type[bb] <- "binomial"}
    }
    if("gaussian" %in% type){type[type == "gaussian"] <- "g"}
    if("binomial" %in% type){
      type[type == "binomial"] <- "c"
      if(all(type == "c")){output$call[c("center", "scale")] <- center <- scale <- FALSE}
    }
    output$call$type <- type
  }

  ### CROSS-SECTIONAL
  if(is.null(lags)){

    # Remove lags from output
    output$call$lags <- NULL

    # Ensure THRESHOLD is not a character
    if(is.character(threshold)){output$call$threshold <- threshold <- FALSE}

    # Update output with names of covariates
    if(!is.null(covariates)){
      if(class(covariates) == "list"){covs <- names(covariates)}
      if(class(covariates) %in% c("numeric", "integer")){covs <- colnames(data)[covariates]}
      output$call$covariates <- covs
    }

    # If moderators, collect names and remove 'pcor' option; otherwise check pcor option
    if(!is.null(moderators)){
      output$call$moderators <- ifelse(is.list(moderators), list(names(moderators)),
                                       list(colnames(data)[moderators]))[[1]]
      if("pcor" %in% names(output$call)){output$call$pcor <- NULL}
    } else if("pcor" %in% names(output$call)){
      if(pcor != FALSE){
        if(threshold == FALSE){output$call$pcor <- TRUE}
        output$call$rule <- NULL
      }
    }

    # Resolve having both moderators and covariates specified
    if(!is.null(covariates) & !is.null(moderators)){
      if(!is.list(covariates) & !is.list(moderators)){
        if(max(moderators) > min(covariates)){
          if(is.character(type)){output$call$type <- type <- type[-covariates]}
          covs <- list(data[, covariates])
          data <- data[, -covariates]
          names(covs) <- output$call$covariates
          covariates <- covs
          moderators <- which(colnames(data) %in% output$call$moderators)
        }
      }
    }

    # NODEWISE -- COMPUTE REGRESSIONS
    mods0 <- nodewise(data = data, mods = moderators, varMods = type,
                      lambda = which.lam, center = center, scale = scale,
                      covariates = covariates, exogenous = exogenous)

    # Provide warning if too many parameters
    if(any(sapply(mods0$models, function(z) length(coef(z))) >= nrow(data))){
      warning('Model is overspecified; not enough cases to estimate all parameters')
    }

    # MODNET -- FORM THE FINAL NETWORK MODEL
    out <- modNet(models = mods0, threshold = threshold, rule = rule,
                  mval = mval, pcor = pcor, binarize = binarize)

    # Append the result to the output list
    output <- append(output, out)

    # Then add mods0?
    output$mods0 <- mods0

    # I'll have to look into this...
    if("moderator" %in% names(attributes(out))){
      attr(output, "moderator") <- attr(out, "moderator")
      output$call$exogenous <- attr(mods0$models, "exogenous")
    }

    # And this...
    if("mval" %in% names(attributes(out))){attr(output, "mval") <- attr(out, "mval")}

    # Set an attribute
    attributes(output)$ggm <- TRUE

    # Compute the log-likelihood
    if(getLL){
      output$modLL <- tryCatch({modLL(output, all = TRUE)}, error = function(e){list()})
      if(length(output$modLL) == 0){
        output$modLL <- NULL
      } else {
        dd <- output$data
        output[c('data', 'mods0')] <- NULL
        output$data <- dd
        output$mods0 <- mods0
      }
    }

    # Update the output classes
    class(output) <- c('list', 'ggm')

    # Prune off models; fitCoefs?
    if(fitCoefs | !saveMods){
      output$mods0$models <- NULL
      output$fitobj <- if(fitCoefs){getFitCIs(output)} else {NULL}
    }

    # Print the time
    if(verbose){print(Sys.time() - t1)}

    # Release the kraken
    return(output)

  } else {
    ### TEMPORAL
    if(!is.null(beepno) & !is.null(dayno)){
      beepday <- list(beepno, dayno)
      stopifnot(sum(sapply(c(1, nrow(data)), function(z) all(sapply(beepday, length) == z))) == 1)
      if(all(sapply(beepday, length) == 1)){
        if(is.character(beepno)){beepno <- which(colnames(data) == beepno)}
        if(is.character(dayno)){dayno <- which(colnames(data) == dayno)}
        data0 <- data[, -c(beepno, dayno)]
        beepno <- data[, beepno]
        dayno <- data[, dayno]
        data <- data0
        if(exists("samp_ind", inherits = FALSE)){
          attr(data, "samp_ind") <- samp_ind
        }
      }
    } else if(!is.null(beepno) & is.null(dayno) | !is.null(dayno) & is.null(beepno)){
      stop('Must specify both beepno AND dayno, or neither')
    }
    if(!identical(detrend, FALSE)){
      if(is(detrend, 'list')){
        stopifnot(length(detrend) %in% 1:2)
        if(is.null(names(detrend))){
          if(length(detrend) == 1){
            detrend <- detrend[[1]]
          } else {
            names(detrend)[order(sapply(detrend, length))] <- c('timevar', 'vars')
          }
        }
      }
      timevar <- switch(2 - isTRUE(detrend), NULL, ifelse(
        is(detrend, 'list'), detrend$timevar, ifelse(
          is(detrend, 'numeric'), colnames(data)[detrend], detrend)))
      dvars <- switch(2 - is(detrend, 'list'), detrend$vars, NULL)
      data <- detrender(data = data, timevar = timevar, vars = dvars, verbose = verbose)
      if(exists("samp_ind", inherits = FALSE)){
        attr(data, "samp_ind") <- samp_ind
      }
    }
    if(!is.null(beepno) & !is.null(dayno)){
      #consec <- mgm:::beepday2consec(beepvar = beepno, dayvar = dayno)
      #consec <- mgm:::lagData(data = data, lags = 1, consec = consec)[[3]][-1]
      consec <- makeConsec(beepvar = beepno, dayvar = dayno)
      consec <- lagData(data = data, lags = 1, consec = consec)[[3]][-1]
      output$call$consec <- which(consec)
    } else {
      consec <- NULL
    }
    output$call$rule <- NULL
    if(threshold != FALSE){output$call$pcor <- ifelse(is.logical(pcor), "none", pcor)}
    if(class(type) == "list"){
      if(length(type) == 1){stop("Need more than one outcome to construct network")}
      if(!attr(type, "method") %in% c("regsubsets", "glmnet")){
        exogenous <- attr(type, "exogenous")
        moderators <- attr(type, "moderators")
        if(!is.character(moderators)){
          attr(type, "moderators") <- moderators <- colnames(data)[moderators]
        }
        anymods <- any(sapply(moderators, function(z){
          any(grepl(paste0(z, ":|:", z), unlist(output$call$type)))
        }))
        moderators <- which(colnames(data) %in% moderators)
        if("covs" %in% names(attributes(type))){covariates <- attr(type, "covs")}
      }
    }
    if(!is.null(union(moderators, covariates))){
      mco <- sort(union(moderators, covariates))
      stopifnot(identical(mco, sort(c(moderators, covariates))))
      if(exogenous & length(mco) >= ncol(data) - 1){exogenous <- FALSE}
    }
    if(!is.null(moderators)){
      output$call$moderators <- colnames(data)[moderators]
      output$call$exogenous <- exogenous
    }
    if(exists("anymods", inherits = FALSE)){
      if(!anymods){
        covariates <- union(covariates, moderators)
        moderators <- NULL
      }
    }
    fit <- SURfit(data = data, varMods = type, m = moderators, mod = which.lam,
                  center = center, scale = scale, exogenous = exogenous,
                  covs = covariates, sur = std, maxiter = maxiter, consec = consec)
    dat <- lagMat(data = data, type = type, m = moderators, covariates = covariates,
                  center = center, scale = scale, exogenous = exogenous, consec = consec)
    net <- SURnet(fit = fit, dat = dat, s = residMat, m = moderators, pcor = pcor,
                  threshold = threshold, mval = mval, medges = medges)
    if(!is.null(covariates)){
      output$call$covariates <- net$call$covariates <- colnames(data)[covariates]}
    output$SURnet <- net
    attributes(output$SURnet)$SURnet <- TRUE
    output$SURfit <- fit
    if("mnet" %in% names(net)){
      attr(output, "mnet") <- net$call$moderators
      if("mval" %in% names(net$call)){attr(output, "mval") <- net$call$mval}
    }
    attributes(output)$SURnet <- TRUE
    if(getLL){
      output$SURll <- tryCatch({SURll(output, all = TRUE, s = residMat)}, error = function(e){list()})
      if(length(output$SURll) == 0){output$SURll <- NULL}
    }
    attr(output, 'rank') <- sum(sapply(lapply(output$SURnet$mods, '[[', 'model'), nrow))
    class(output) <- c('list', 'SURnet')
    if(fitCoefs | !saveMods){
      output$SURfit <- if(fitCoefs){getFitCIs(output)} else {NULL}
    }
    if(verbose){print(Sys.time() - t1)}
    return(output)
  }
}

##### nodewise: nodewise regression, with option to simulate/add moderators
nodewise <- function(data, mods = NULL, varMods = NULL, lambda = "min", center = TRUE,
                     scale = FALSE, covariates = NULL, exogenous = TRUE){
  # Good lord...
  data <- dat <- dat0 <- data.frame(data)

  # Save them colnames
  vs <- colnames(data)

  # Variable selection/TYPE
  if(!is.null(varMods)){

    # Variable selection component
    if(class(varMods) == "list"){
      if(all(sapply(varMods, class) != "list")){
        varMods <- lapply(varMods, list)
        lambda <- "min"
      }
      if(is.null(names(varMods))){
        if((is.null(mods) | (!is.null(mods) & class(mods) == "list")) &
           (is.null(covariates) | (!is.null(covariates) & class(covariates) == "list"))){
          names(varMods) <- colnames(data)
        } else {
          if(!is.null(mods) & class(mods) %in% c("numeric", "integer")){mc <- mods} else {mc <- NULL}
          if(!is.null(covariates) & class(covariates) %in% c("numeric", "integer")){mc <- c(mc, covariates)}
          names(varMods) <- colnames(data)[-mc]
        }
      }
      if(!is.null(mods)){if(length(mods) > 1){mods <- NULL}}
      if(length(varMods) == ncol(data)){exogenous <- FALSE}
      if(is.null(unlist(lapply(varMods, '[[', "mod1se")))){lambda <- "min"}
      lambda <- ifelse(lambda == "min", 1, 2)
      type <- unname(sapply(varMods, attr, "family"))
      if(all(type == "c")){center <- FALSE}

      # TYPE component
    } else {
      type <- varMods
      if(is.null(mods)){
        varMods <- NULL
      } else if(length(mods) == 1){
        varMods <- NULL
      } else {
        varMods <- setNames(lapply(1:ncol(data), function(z){
          mains <- vs[-z]
          ints <- apply(combn(mains, 2), 2, paste0, collapse = ":")
          if(!z %in% mods){ints <- ints[apply(combn(mains, 2), 2, function(zz) any(vs[mods] %in% zz))]}
          return(list(mod0 = c(mains, ints)))
        }), vs)
        if(exogenous & length(mods) < (ncol(data) - 1)){
          varMods <- varMods[vs[-mods]]
        } else {
          exogenous <- FALSE
        }
        mods <- NULL
        lambda <- 1
      }
    }
  } else if(!is.null(mods)){
    if(length(mods) > 1){
      varMods <- lapply(1:ncol(data), function(z){
        mains <- vs[-z]
        ints <- apply(combn(mains, 2), 2, paste0, collapse = ":")
        if(!z %in% mods){ints <- ints[apply(combn(mains, 2), 2, function(zz) any(vs[mods] %in% zz))]}
        return(list(mod0 = c(mains, ints)))
      })
      names(varMods) <- vs
      type <- unname(ifelse(apply(data, 2, function(z) dim(table(z))) <= 2, "c", "g"))
      exogenous <- FALSE
      mods <- NULL
      lambda <- 1
    }
  }

  # intMatrix function
  intMatrix <- function(data, mods = NULL, covariates = NULL){
    if(class(mods) == "list"){stopifnot(!is.null(names(mods)))}
    data <- as.data.frame(data)
    vars <- colnames(data)
    eqs <- list()
    if(!is.null(mods)){
      for(i in 1:length(vars)){
        modTerms <- list()
        for(j in 1:length(mods)){modTerms[[j]] <- paste0(names(mods[j]), " + ", paste(paste0(vars[-i], ":", names(mods[j])), collapse = " + "))}
        modTerms <- paste(modTerms, collapse = " + ")
        eqs[[i]] <- paste0(vars[i], " ~ ", paste(vars[-i], collapse = " + "), " + ", modTerms)
      }
    } else {
      for(i in 1:length(vars)){eqs[[i]] <- paste0(vars[i], " ~ ", paste(vars[-i], collapse = " + "))}
    }
    if(!is.null(covariates)){
      for(i in 1:length(vars)){eqs[[i]] <- paste0(eqs[[i]], " + ", paste(names(covariates), collapse = " + "))}
    }
    eqs <- setNames(lapply(eqs, as.formula), vars)
    return(eqs)
  }

  # Covariates
  if(!is.null(covariates)){
    if(class(covariates) %in% c("numeric", "integer")){
      if(length(covariates) > 1){
        covs <- as.list(data[, covariates])
      } else {
        covs <- list(data[, covariates])
        names(covs) <- colnames(data)[covariates]
      }
      data <- dat <- dat0 <- data[, -covariates]
      covariates <- covs
    }
  }

  # Moderators
  if(!is.null(mods)){
    if(class(mods) %in% c("numeric", "integer")){
      mod <- list(data[, mods])
      names(mod) <- colnames(data)[mods]
      data <- data[, -mods]
      if(length(type) > ncol(data)){type <- c(type[-mods], type[mods])}
      mods <- mod
    }
    if(length(mods) != 1){stop("Cannot specify more than on exogenous moderator")}
    dat <- dat0 <- data.frame(data, mods)
  }

  # Centering and standardizing
  if(center != FALSE){
    binary <- unname(which(apply(dat, 2, function(z) dim(table(z)) <= 2)))
    if(length(binary) == ncol(dat) | ifelse(!is.null(mods), length(binary) == ncol(dat) - 1, FALSE)){
      type <- "binomial"
    }
    if(!is.null(mods) & dim(table(mods[[1]])) <= 2 | center != TRUE){
      if(length(binary) > 0){
        dat[, -union(binary, ncol(dat))] <- apply(dat[, -union(binary, ncol(dat))], 2, scale, TRUE, scale)
      } else {
        dat[, -ncol(dat)] <- apply(dat[, -ncol(dat)], 2, scale, TRUE, scale)
      }
    } else {
      if(length(binary) > 0){
        dat[, -binary] <- apply(dat[, -binary], 2, scale, TRUE, scale)
      } else {
        dat <- apply(dat, 2, scale, TRUE, scale)
      }
    }
    if(!is.null(covariates)){
      covariates <- lapply(covariates, function(z){
        ifelse(dim(table(z)) <= 2 | center != TRUE, list(z), list(scale(z, TRUE, scale)))[[1]]
      })
    }
    dat <- dat0 <- data.frame(dat)
  }

  # Return the covariates
  if(!is.null(covariates)){dat <- data.frame(dat, covariates)}

  # Variable selection?
  if(!is.null(varMods)){
    ints <- as.list(paste(names(varMods), "~", lapply(varMods, function(z) paste0(z[[lambda]], collapse = " + "))))
    if(!is.null(covariates) & ifelse("covariates" %in% names(attributes(varMods)), FALSE, TRUE)){
      ints <- lapply(ints, paste0, " + ", paste(names(covariates), collapse = " + "))
    }
    names(ints) <- names(varMods)
  } else {
    ints <- intMatrix(data, mods, covariates)
    if(!is.null(mods) & !exogenous){ints <- append(ints, list(as.formula(paste0(names(mods), " ~ .^2"))))}
  }

  # Check TYPE, FIT MODELS
  if(exists("type", inherits = FALSE)){
    if(length(type) == 1){
      type <- rep(match.arg(type, c("g", "c", "gaussian", "binomial")), length(ints))
    }
    if(any(type %in% c("g", "c"))){
      type <- unname(sapply(type, switch, "g" = "gaussian", "c" = "binomial"))
    }
    m <- suppressWarnings(lapply(1:length(ints), function(z){
      #if(type[z] == "gaussian"){mglm <- lm(ints[[z]], dat)}
      #if(type[z] == "binomial"){mglm <- glm(ints[[z]], data = dat, family = type[z])}
      mglm <- glm(ints[[z]], data = dat, family = type[z])
      attr(mglm, "family") <- type[z]
      return(mglm)
    }))
  } else {
    m <- lapply(ints, function(z) lm(z, dat))
  }

  # Not sure what's going on here...
  if(!is.null(mods) & !exogenous){
    data0 <- data
    data <- dat0
  }

  # Collect coefficients
  mm <- lapply(lapply(m, coef), function(z) z[which(names(z) %in% colnames(data))])

  # Pvalues for coefficients
  ps1 <- lapply(m, function(z){
    z1 <- is.na(coef(z))
    z2 <- t(data.frame(summary(z)$coefficients))
    z0 <- ifelse(any(z1), list(names(z1[!z1])), list(names(coef(z))))[[1]]
    z3 <- z2[4, which(z0 %in% colnames(data)), drop = FALSE]
    rownames(z3) <- NULL
    return(z3[1, ])
  })

  # Pvalues when moderators included
  if(!is.null(mods)){
    mname <- names(mods)
    m2 <- lapply(lapply(m, coef), function(z) z[grep(paste0(":", mname, "$"), names(z))])
    ps2 <- lapply(m, function(z){
      z2 <- t(data.frame(summary(z)$coefficients))
      z3 <- z2[4, grep(paste0(":", mname, "$"), colnames(z2)), drop = FALSE]
      rownames(z3) <- NULL
      return(z3[1, ])
    })
    vs <- ifelse(exogenous, list(colnames(data)), list(colnames(data0)))[[1]]
    mx <- paste0(vs, ":", mname)
    psx <- bx <- suppressWarnings(diag(mx))
    rownames(psx) <- rownames(bx) <- vs
    colnames(psx) <- colnames(bx) <- mx
    diag(psx) <- diag(bx) <- 1
    m3 <- lapply(m, function(z) summary(z)$coefficients)
    names(m3) <- vs
    bm <- do.call(rbind, lapply(m3, function(z) z[rownames(z) == mname, ]))
    if(!exogenous){
      m2 <- m2[-ncol(data)]
      ps2 <- ps2[-ncol(data)]
      m3 <- m3[-ncol(data)]
    }
  }

  # More summary when covariates are included
  if(!is.null(covariates)){
    m4 <- lapply(m, function(z) summary(z)$coefficients)
    names(m4) <- ifelse(!is.null(mods), list(vs), list(colnames(data)))[[1]]
    cnames <- names(covariates)
    bcovs <- list()
    for(i in 1:length(cnames)){
      bcovs[[i]] <- do.call(rbind, lapply(m4, function(z) z[rownames(z) == cnames[i], ]))
    }
    names(bcovs) <- cnames
  }

  #b <- ps <- matrix(NA, nrow = ncol(data), ncol = ncol(data))
  #for(i in 1:ncol(data)){
  #  b[i, match(names(mm[[i]]), colnames(data))] <- mm[[i]]
  #  ps[i, match(names(ps1[[i]]), colnames(data))] <- ps1[[i]]
  #  if(!is.null(mods)){
  #    if(exogenous | (!exogenous & i < ncol(data))){
  #      bx[i, match(names(m2[[i]]), mx)] <- m2[[i]]
  #      psx[i, match(names(ps2[[i]]), mx)] <- ps2[[i]]
  #    }
  #  }
  #}
  #for(i in 1:length(ints)){
  #  b[i, match(names(mm[[i]]), names(ints))] <- mm[[i]]
  #  ps[i, match(names(ps1[[i]]), names(ints))] <- ps1[[i]]
  #  if(!is.null(mods)){
  #    if(exogenous | (!exogenous & i < ncol(data))){
  #      bx[i, match(names(m2[[i]]), mx)] <- m2[[i]]
  #      psx[i, match(names(ps2[[i]]), mx)] <- ps2[[i]]
  #    }
  #  }
  #}

  # Setup results in matrix format
  b <- ps <- matrix(NA, nrow = length(ints), ncol = length(ints))
  for(i in 1:length(ints)){
    mm0 <- mm[[i]][which(names(mm[[i]]) %in% names(ints))]
    ps0 <- ps1[[i]][which(names(ps1[[i]]) %in% names(ints))]
    b[i, match(names(mm0), names(ints))] <- mm0
    ps[i, match(names(ps0), names(ints))] <- ps0
    if(!is.null(mods)){
      if(exogenous | (!exogenous & i < ncol(data))){
        bx[i, match(names(m2[[i]]), mx)] <- m2[[i]]
        psx[i, match(names(ps2[[i]]), mx)] <- ps2[[i]]
      }
    }
  }
  diag(b) <- diag(ps) <- 1
  if(any(is.na(b))){b[is.na(b)] <- 0}
  if(any(is.na(ps))){ps[is.na(ps)] <- 0}
  out <- list(models = setNames(m, names(ints)), B = list(b = b, ps = ps))

  # Set attributes for moderators
  if(is.null(mods)){
    attributes(out$models)$noMods <- TRUE
    if(length(ints) != ncol(data)){attr(out$models, 'exogenous') <- TRUE}
  } else {
    attributes(out$models)$exogenous <- exogenous
    if(nrow(bm) >= 1 & exogenous){out$Bm <- bm}
    out$Bx <- list(bx = bx, px = psx)
  }

  # Not sure if this is about moderators or variable selection, or both
  if(!is.null(varMods)){
    if(!any(grepl(":", unlist(sapply(out$models, function(z) names(coef(z))))))){
      attributes(out$models)$noMods <- TRUE
      out$Bx <- NULL
    }
    attributes(out$models)$varMods <- c("min", "1se")[lambda]
  }

  # Setup covariate output
  if(!is.null(covariates)){
    out$dat <- dat0
    out$covariates <- list(covs = do.call(cbind.data.frame, covariates), Bcovs = bcovs)
  } else {
    out$dat <- dat
  }

  # Gotta figure out this "exists type" thing
  if(exists("type", inherits = FALSE)){attr(out, "type") <- type}

  # RETURN
  return(out)
}

##### modNet: create moderated network from nodewise regression models
modNet <- function(models, data = NULL, threshold = FALSE, rule = "AND", mval = NULL,
                   pcor = FALSE, useCIs = FALSE, nsims = 5000, mlty = 2, binarize = FALSE){
  # Get models and data
  if("models" %in% names(models)){
    mods0 <- models
    models <- models$models
    if(is.null(data)){
      if("noMods" %in% names(attributes(models)) & length(models) == ncol(mods0$dat)){
        data <- mods0$dat
      } else if(isTRUE(attr(models, "exogenous"))){
        #data <- mods0$dat[, -ncol(mods0$dat)]
        data <- mods0$dat[, names(models)]
      } else {
        data <- mods0$dat
      }
    }
  }

  # Setup
  p <- ncol(data)
  n <- nrow(data)
  vs <- colnames(data)

  # Collect coefficients
  mods <- lapply(models, function(z){
    z2 <- matrix(coef(z), ncol = 1)
    rownames(z2) <- names(coef(z))
    return(z2)
  })

  # Removing intercept? And turning matrices into vectors
  mods2 <- lapply(mods, function(z) z[which(rownames(z) %in% vs), ])

  # Get pvalues
  pvals <- lapply(models, function(z){
    z1 <- is.na(coef(z))
    z2 <- t(data.frame(summary(z)$coefficients))
    z0 <- ifelse(any(z1), list(names(z1[!z1])), list(names(coef(z))))[[1]]
    z3 <- z2[4, which(z0 %in% vs), drop = FALSE]
    rownames(z3) <- NULL
    return(z3[1, ])
  })

  # Get standard errors
  ses <- lapply(models, function(z){
    z1 <- is.na(coef(z))
    z2 <- t(data.frame(summary(z)$coefficients))
    z0 <- ifelse(any(z1), list(names(z1[!z1])), list(names(coef(z))))[[1]]
    z3 <- z2[2, which(z0 %in% vs), drop = FALSE]
    rownames(z3) <- NULL
    return(z3[1, ])
  })

  # Turn them all into matrices
  b <- matrix(0, p, p)
  pvals2 <- ses2 <- matrix(1, p, p)
  for(i in 1:p){
    b[i, match(names(mods2[[i]]), vs)] <- mods2[[i]]
    ses2[i, match(names(ses[[i]]), vs)] <- ses[[i]]
    pvals2[i, match(names(pvals[[i]]), vs)] <- pvals[[i]]
  }

  # Compute results for each model separately
  results <- lapply(1:p, function(z){
    notype <- !"type" %in% names(attributes(mods0))
    if(notype | ifelse(!notype, ifelse(attr(mods0, "type")[z] == "gaussian", TRUE, FALSE), TRUE)){
      yhat <- predict(models[[z]])
      deviance <- sum((data[, z] - yhat)^2)
      s <- sqrt(deviance/n)
      LL_model <- sum(dnorm(data[, z], mean = yhat, sd = s, log = TRUE))
      k <- nrow(mods[[z]]) + 1
      aic <- (2 * k) - (2 * LL_model)
      bic <- (log(n) * k) - (2 * LL_model)
    } else {
      deviance <- deviance(models[[z]])
      LL_model <- as.numeric(logLik(models[[z]]))
      aic <- AIC(models[[z]])
      bic <- BIC(models[[z]])
    }
    pees <- as.matrix(summary(models[[z]])$coefficients[, 4], ncol = 1) ###
    return(list(deviance = deviance, LL_model = LL_model, AIC = aic, BIC = bic, model = mods[[z]], pvals = pees)) ###
  })

  # Not sure how this works...
  if(!"noMods" %in% names(attributes(models))){pcor <- FALSE}

  # Formatting base on whether moderators are included? Not sure... maybe related to variable selection
  if("noMods" %in% names(attributes(models))){
    inds <- ints <- mval <- vars1 <- intMats <- NULL
  } else if("varMods" %in% names(attributes(models))){
    nmods <- ifelse(attr(models, "exogenous") == TRUE, length(mods), length(mods) - 1)
    inds0 <- lapply(mods, function(z) rownames(z)[grepl(":", rownames(z))])[1:nmods]
    inds1 <- unlist(lapply(inds0, function(z) gsub(":.*", "", z)))
    inds <- data.frame(outcome = rep(1:nmods, sapply(inds0, length)), interaction = unlist(inds0))
    vars <- lapply(models, vcov)[1:nmods]
    vars1 <- ints <- list()
    for(i in 1:nrow(inds)){
      ints[[i]] <- mods[[inds[i, 1]]][inds[i, 2], 1]
      vars1[[i]] <- c(vars[[inds[i, 1]]][inds1[i], inds1[i]], vars[[inds[i, 1]]][inds[i, 2], inds[i, 2]], vars[[inds[i, 1]]][inds1[i], inds[i, 2]])
      names(vars1[[i]]) <- c(inds1[i], inds[i, 2], "cov")
    }
  } else {
    nmods <- ifelse(attr(models, "exogenous") == TRUE, length(mods), length(mods) - 1)
    inds <- t(combn(1:nmods, 2))
    ints <- vector("list", nrow(inds))
    for(i in 1:nrow(inds)){
      ints[[i]][1] <- mods[[inds[i, 1]]][nmods + inds[i, 2], ]
      ints[[i]][2] <- mods[[inds[i, 2]]][nmods + inds[i, 1] + 1, ]
    }
    vars <- lapply(models, vcov)[1:nmods]
    vars1 <- list()
    for(i in 1:nmods){
      vars1[[i]] <- vector("list", nmods - 1)
      for(j in 1:(nmods - 1)){
        vars1[[i]][[j]] <- c(vars[[i]][j + 1, j + 1], vars[[i]][j + nmods + 1, j + nmods + 1], vars[[i]][j + 1, j + nmods + 1])
        names(vars1[[i]][[j]]) <- c(colnames(vars[[i]])[j + 1], colnames(vars[[i]])[j + nmods + 1], "cov")
      }
      names(vars1[[i]]) <- colnames(vars[[i]])[2:nmods]
    }
    names(vars1) <- colnames(data)[1:nmods]
  }

  # Function for converting betas to GGM
  b2ggm <- function(b, rule = "AND", pcor = FALSE, threshold = FALSE, n = NULL){
    rule <- match.arg(tolower(rule), c("and", "or"))
    if(rule == "and" & pcor == FALSE){
      bb <- cbind(b[upper.tri(b)], t(b)[upper.tri(t(b))])
      notBoth <- !apply(bb, 1, function(z) (z[1] == 0 & z[2] == 0) | (z[1] != 0 & z[2] != 0))
      if(any(notBoth)){
        bb[notBoth, ] <- 0
        b[upper.tri(b)] <- bb[, 1]
        b <- t(b)
        b[upper.tri(b)] <- bb[, 2]
        b <- t(b)
      }
    }
    if(pcor != FALSE){
      bb <- sign(b) * sqrt(b * t(b))
      if(pcor == 'cor'){
        diag(bb) <- 1
        bb <- corpcor::pcor2cor(bb)
      } else if(grepl('cor_auto', pcor)){
        bb <- qgraph::cor_auto(data, npn.SKEPTIC = TRUE)
        if(pcor == 'cor_auto2'){bb <- corpcor::cor2pcor(bb)}
      }
      diag(bb) <- 0
      if(threshold != FALSE){
        pcor <- ifelse(isTRUE(pcor) | grepl('cor', pcor), 'none', pcor)
        if(is.character(threshold)){pcor <- threshold}
        if(!is.numeric(threshold)){threshold <- .05}
        dimnames(bb) <- rep(list(paste0("X", 1:ncol(bb))), 2)
        sigMat <- ifelse(psych::corr.p(bb, n, adjust = pcor)[[4]] <= threshold, 1, 0)
        sigMat0 <- matrix(0, ncol(bb), ncol(bb))
        sigMat0[upper.tri(sigMat0)] <- sigMat[upper.tri(sigMat)]
        sigMat0 <- as.matrix(Matrix::forceSymmetric(sigMat0))
        bb <- bb * sigMat0
        bb <- unname(bb)
      }
      return(bb)
    } else {
      return((b + t(b))/2)
    }
  }

  # If there are moderators in the model, AND you want a conditional network
  if(!is.null(mval)){
    margSE <- function(x, vars){sqrt(vars[1] + ((x^2) * vars[2]) + (2 * x * vars[3]))}
    if("varMods" %in% names(attributes(models))){
      inds1.1 <- match(inds1, vs)
      for(i in 1:length(ints)){
        b[inds[i, 1], inds1.1[i]] <- b[inds[i, 1], inds1.1[i]] + (mval * ints[[i]])
        ses[[inds[i, 1]]][inds1[i]] <- margSE(mval, vars1[[i]])
      }
      for(i in 1:p){ses2[i, match(names(ses[[i]]), vs)] <- ses[[i]]}
    } else {
      for(i in 1:length(ints)){
        b[inds[i,1], inds[i,2]] <- b[inds[i,1], inds[i,2]] + (mval * ints[[i]][1])
        b[inds[i,2], inds[i,1]] <- b[inds[i,2], inds[i,1]] + (mval * ints[[i]][2])
      }
      inds2 <- rbind(inds, cbind(inds[, 2], inds[, 1]))
      inds2 <- inds2[order(inds2[, 1]), ]
      inds3 <- cbind(inds2[, 1], rep(c(1:(p - 1)), p))
      for(i in 1:nrow(inds3)){ses[[inds3[i, 1]]][inds3[i, 2]] <- margSE(mval, vars1[[inds3[i, 1]]][[inds3[i, 2]]])}
      for(i in 1:p){ses2[i, match(names(ses[[i]]), vs)] <- ses[[i]]}
    }
    bmval <- b
    dimnames(bmval) <- rep(list(colnames(data)), 2)
    dfs <- matrix(sapply(models, function(z) z$df.residual), ncol = p, nrow = p)
    pvals2 <- (2 * pt(abs(b/ses2), df = dfs, lower.tail = FALSE))
    if(any(is.na(pvals2))){pvals2[is.na(pvals2)] <- 1}
  }

  # Thresholding!
  if(threshold != FALSE & pcor == FALSE){
    if(threshold == TRUE){threshold <- .05}
    b <- b * ifelse(pvals2 <= threshold, 1, 0)
  }

  # PRODUCE THE NETWORK
  bb <- b2ggm(b, rule = rule, pcor = pcor, threshold = threshold, n = n)

  # getEdgeColors function -- CHECK
  getEdgeColors <- function(adjMat){
    obj <- sign(as.vector(adjMat))
    colMat <- rep(NA, length(obj))
    if(any(obj == 1)){colMat[obj == 1] <- "darkgreen"}
    if(any(obj == 0)){colMat[obj == 0] <- "darkgrey"}
    if(any(obj == -1)){colMat[obj == -1] <- "red"}
    colMat <- matrix(colMat, ncol = ncol(adjMat), nrow = nrow(adjMat))
    colnames(colMat) <- colnames(adjMat)
    rownames(colMat) <- rownames(adjMat)
    colMat
  }

  # Making the edge colors; more uncertainty around noMods
  if(!"noMods" %in% names(attributes(models))){
    if(useCIs & nsims > 0 & attr(models, "exogenous") == TRUE){
      cis <- margCIs(mods = mods0, alpha = ifelse(threshold == FALSE, .05, threshold), nsims = nsims)
      modEdgesNW <- getInts(x = cis, allInts = TRUE)
    } else {
      modEdgesNW <- ifelse(mods0$Bx$px <= ifelse(threshold == FALSE, .05, threshold), 1, 0)
      if(any(mods0$Bx$px == 0)){modEdgesNW <- modEdgesNW * ifelse(mods0$Bx$px == 0, 0, 1)}
      colnames(modEdgesNW) <- rownames(modEdgesNW)
    }
    modEdges <- t(modEdgesNW) * modEdgesNW
    modEdgesNW <- (modEdgesNW * (mlty - 1)) + 1
    modEdges <- (modEdges * (mlty - 1)) + 1
    rule <- match.arg(tolower(rule), c("and", "or"))
    if(rule == "or"){
      modEdges <- modEdgesNW + t(modEdgesNW)
      modEdges[modEdges == 2] <- 1
      modEdges[modEdges > 2] <- mlty
    }
    intMat1 <- mods0$Bx$bx
    intPs <- mods0$Bx$px
    if(any(intPs == 0)){intPs[intPs == 0] <- 1}
    if(threshold != FALSE){intMat1 <- intMat1 * ifelse(intPs <= threshold, 1, 0)}
    intMat2 <- b2ggm(intMat1, rule = rule, pcor = FALSE)
    diag(intMat1) <- diag(intMat2) <- 0
    intMats <- list(avgInts = t(intMat2), nwInts = list(adj2NW = t(intMat1), pvals2NW = t(intPs)))
  } else {
    modEdges <- modEdgesNW <- matrix(1, p, p)
  }

  # Format results
  diag(modEdges) <- diag(modEdgesNW) <- 0
  dimnames(b) <- dimnames(bb) <- dimnames(pvals2) <- rep(list(colnames(data)), 2)
  names(mods) <- names(models) <- names(results) <- colnames(data)
  out <- list(adjMat = bb, edgeColors = getEdgeColors(bb), modEdges = t(modEdges),
              nodewise = list(adjNW = t(b), edgeColsNW = getEdgeColors(t(b)),
                              pvalsNW = t(pvals2), modEdgesNW = t(modEdgesNW)),
              interactions = list(intMats = intMats, inds = inds, ints = ints, coefvars = vars1))

  # Making binary edges grey
  if(binarize){
    out$adjMat[out$adjMat != 0] <- 1
    out$edgeColors[out$edgeColors != 'darkgrey'] <- 'darkgrey'
  }

  # More noMods uncertainty...
  if("noMods" %in% names(attributes(models))){
    out[c("interactions", "modEdges")] <- NULL
    out[["nodewise"]]["modEdgesNW"] <- NULL
  } else if(attr(models, "exogenous") == FALSE){
    out[["modEdges"]] <- rbind(cbind(out[["modEdges"]], 1), 1)
    out[["nodewise"]][["modEdgesNW"]] <- rbind(cbind(out[["nodewise"]][["modEdgesNW"]], 1), 1)
    colnames(out$modEdges)[p] <- rownames(out$modEdges)[p] <- colnames(mods0$dat)[p]
    colnames(out$nodewise$modEdgesNW)[p] <- rownames(out$nodewise$modEdgesNW)[p] <- colnames(mods0$dat)[p]
    attributes(out)$moderator <- colnames(mods0$dat)[ncol(mods0$dat)]
  } else {
    if(useCIs){out$interactions$cis <- cis}
    attributes(out)$moderator <- colnames(mods0$dat)[ncol(mods0$dat)]
    if("Bm" %in% names(mods0)){
      mp <- ncol(bb)
      madj <- medges <- matrix(0, mp + 1, mp + 1)
      md <- matrix(FALSE, mp + 1, mp + 1)
      madj[1:mp, 1:mp] <- bb
      dimnames(madj) <- dimnames(medges) <- dimnames(md) <- rep(list(colnames(mods0$dat)), 2)
      if(nrow(mods0$Bm) != p){mpp <- which(colnames(data) %in% rownames(mods0$Bm))} else {mpp <- 1:mp}
      if(threshold != FALSE){
        madj[(mp + 1), mpp] <- ifelse(mods0$Bm[, 4] <= threshold, mods0$Bm[, 1], 0)
      } else {
        madj[(mp + 1), mpp] <- mods0$Bm[, 1]
      }
      medges[1:mp, 1:mp] <- t(modEdges)
      medges[(mp + 1), 1:mp] <- 1
      md[(mp + 1), 1:mp] <- TRUE
      ints <- out$interactions
      out$interactions <- NULL
      out$mnet <- list(adjMat = madj, edgeColors = getEdgeColors(madj), modEdges = medges, d = md)
      out$interactions <- ints
    }
  }

  # Format output
  out$mods <- results
  out$fitobj <- models
  out$data <- data
  if(!is.null(mval)){
    attributes(out)$mval <- mval
    if(threshold != FALSE){out$nodewise[[paste0("adjMval", mval)]] <- t(bmval)}
  }

  # Set ggm attribute
  attributes(out)$ggm <- TRUE

  # RETURN
  return(out)
}

#' Provides model coefficients with confidence intervals
#'
#' Requires that either \code{fitobj} or \code{SURfit} is included in the object
#' from \code{\link{fitNetwork}}. Returns a list of nodewise model coefficients,
#' including confidence intervals computed from the estimated standard errors.
#'
#' The \code{select} column in the output indicates whether the variable would
#' be selected given the supplied alpha level.
#'
#' @param fit Output from \code{\link{fitNetwork}}, or either the
#'   \code{fixedNets} or \code{betweenNet} element of the output from
#'   \code{\link{mlGVAR}}
#' @param allNames Character vector containing all the predictor names. Do not
#'   change, as these are automatically detected.
#' @param alpha Type 1 error rate. The complement of the confidence level.
#'
#' @return List of tables containing model coefficients along with confidence
#'   intervals
#' @export
#'
#' @seealso \code{\link{fitNetwork}, \link{plotCoefs}}
#'
#' @examples
#' x <- fitNetwork(ggmDat)
#' getFitCIs(x)
getFitCIs <- function(fit, allNames = NULL, alpha = .05){
  if("SURnet" %in% names(fit)){
    if(!'SURfit' %in% names(fit)){stop('Requires SURfit')}
    if(is(fit$SURfit, 'fitCoefs')){return(fit$SURfit)}
    fit <- append(fit$SURnet, list(fitobj = fit$SURfit$eq))
  }
  if(is.null(allNames)){
    allNames <- lapply(fit$mods, function(z) rownames(z$model)[-1])
  } else if(all(c("mods", "call") %in% names(allNames))){
    allNames <- lapply(allNames$mods, function(z) rownames(z$model)[-1])
  }
  if('fitobj' %in% names(fit)){
    fitobj <- fit$fitobj
  } else {
    fitobj <- fit
  }
  if(is(fitobj, 'fitCoefs')){return(fitobj)}
  fitCoefs <- suppressWarnings(suppressMessages(
    lapply(seq_along(fitobj), function(z){
      coefs1 <- summary(fitobj[[z]])$coefficients[-1, , drop = FALSE]
      if(nrow(coefs1) == 0){return(NA)}
      coefs2 <- confint(fitobj[[z]], level = 1 - alpha)[-1, , drop = FALSE]
      coefs <- cbind.data.frame(coefs1, coefs2)
      colnames(coefs) <- c("b", "se", "t", "P", "lower", "upper")
      coefs <- coefs[match(allNames[[z]], rownames(coefs)), ]
      coefs <- data.frame(predictor = factor(allNames[[z]]), lower = coefs$lower,
                          b = coefs$b, upper = coefs$upper, Pvalue = coefs$P,
                          select = ifelse(is.na(coefs$b), FALSE, TRUE))
      rownames(coefs) <- 1:nrow(coefs)
      if(any(is.na(coefs))){coefs[is.na(coefs)] <- NA}
      return(coefs)
    })))
  names(fitCoefs) <- names(fitobj)
  class(fitCoefs) <- c('list', 'fitCoefs')
  return(fitCoefs)
}
