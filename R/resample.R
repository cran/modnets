#' Bootstrapping or multi-sample splits for variable selection
#'
#' Multiple resampling procedures for selecting variables for a final network
#' model. There are three resampling methods that can be parameterized in a
#' variety of different ways. The ultimate goal is to fit models across iterated
#' resamples with variable selection procedures built in so as to home in on the
#' best predictors to include within a given model. The methods available
#' include: bootstrapped resampling, multi-sample splitting, and stability
#' selection.
#'
#' Sampling methods can be specified via the \code{sampMethod} argument.
#' \describe{ \item{Bootstrapped resampling}{Standard bootstrapped resampling,
#' wherein a bootstrapped sample of size \code{n} is drawn with replacement at
#' each iteration. Then, a variable selection procedure is applied to the
#' sample, and the selected model is fit to obtain the parameter values.
#' P-values and confidence intervals for the parameter distributions are then
#' estimated.} \item{Multi-sample splitting}{Involves taking two disjoint
#' samples from the original data -- a training sample and a test sample. At
#' each iteration the variable selection procedure is applied to the training
#' sample, and then the resultant model is fit on the test sample. Parameters
#' are then aggregated based on the coefficients in the models fit to the test
#' samples.} \item{Stability selection}{Stability selection begins the same as
#' multi-sample splitting, in that two disjoint samples are drawn from the data
#' at each iteration. However, the variable selection procedure is then applied
#' to each of the two subsamples at each iteration. The objective is to compute
#' the proportion of times that each predictor was selected in each subsample
#' across iterations, as well as the proportion of times that it was
#' simultaneously selected in both disjoint samples. At the end of the
#' resampling, the final model is selected by setting a frequency threshold
#' between 0 and 1, indicating the minimum proportion of samples that a variable
#' would have to have been selected to be retained in the final model.} }
#'
#' For the bootstrapping and multi-sample split methods, p-values are aggregated
#' for each parameter using a method developed by Meinshausen, Meier, & Buhlmann
#' (2009) that employs error control based on the false-discovery rate. The same
#' procedure is employed for creating adjusted confidence intervals.
#'
#' A key distinguishing feature of the bootstrapping procedure implemented in
#' this function versus the \code{\link{bootNet}} function is that the latter is
#' designed to estimate the parameter distributions of a single model, whereas
#' the version here is aimed at using the bootstrapped resamples to select a
#' final model. In a practical sense, this boils down to using the bootstrapping
#' method in the \code{\link{resample}} function to perform variable selection
#' at each iteration of the resampling, rather than taking a single constrained
#' model and applying it equally at all iterations.
#'
#' @param data \code{n x k} dataframe. Cannot supply a matrix as input.
#' @param m Character vector or numeric vector indicating the moderator(s), if
#'   any. Can also specify \code{"all"} to make every variable serve as a
#'   moderator, or \code{0} to indicate that there are no moderators. If the
#'   length of \code{m} is \code{k - 1} or longer, then it will not be possible
#'   to have the moderators as exogenous variables. Thus, \code{exogenous} will
#'   automatically become \code{FALSE}.
#' @param niter Number of iterations for the resampling procedure.
#' @param sampMethod Character string indicating which type of procedure to use.
#'   \code{"bootstrap"} is a standard bootstrapping procedure. \code{"split"} is
#'   the multi-sample split procedure where the data are split into disjoint
#'   training and test sets, the variables to be modeled are selected based on
#'   the training set, and then the final model is fit to the test set.
#'   \code{"stability"} is stability selection, where models are fit to each of
#'   two disjoint subsamples of the data, and it is calculated how frequently
#'   each variable is selected in each subset, as well how frequently they are
#'   simultaneously selected in both subsets at each iteration.
#' @param criterion The criterion for the variable selection procedure. Options
#'   include: \code{"cv", "aic", "bic", "ebic", "cp", "rss", "adjr2", "rsq",
#'   "r2"}. \code{"CV"} refers to cross-validation, the information criteria are
#'   \code{"AIC", "BIC", "EBIC"}, and \code{"Cp"}, which refers to Mallow's Cp.
#'   \code{"RSS"} is the residual sum of squares, \code{"adjR2"} is adjusted
#'   R-squared, and \code{"Rsq"} or \code{"R2"} is R-squared. Capitalization is
#'   ignored. For methods based on the LASSO, only \code{"CV", "AIC", "BIC",
#'   "EBIC"} are available. For methods based on subset selection, only
#'   \code{"Cp", "BIC", "RSS", "adjR2", "R2"} are available.
#' @param method Character string to indicate which method to use for variable
#'   selection. Options include \code{"lasso"} and \code{"glmnet"}, both of
#'   which use the LASSO via the \code{glmnet} package (either with
#'   \code{\link[glmnet:glmnet]{glmnet::glmnet}} or
#'   \code{\link[glmnet:cv.glmnet]{glmnet::cv.glmnet}}, depending upon the
#'   criterion). \code{"subset", "backward", "forward", "seqrep"}, all call
#'   different types of subset selection using the
#'   \code{\link[leaps:regsubsets]{leaps::regsubsets}} function. Finally
#'   \code{"glinternet"} is used for applying the hierarchical lasso, and is the
#'   only method available for moderated network estimation (either with
#'   \code{\link[glinternet:glinternet]{glinternet::glinternet}} or
#'   \code{\link[glinternet:glinternet.cv]{glinternet::glinternet.cv}},
#'   depending upon the criterion). If one or more moderators are specified,
#'   then \code{method} will automatically default to \code{"glinternet"}.
#' @param rule Only applies to GGMs (including between-subjects networks) when a
#'   threshold is supplied. The \code{"AND"} rule will only preserve edges when
#'   both corresponding coefficients have p-values below the threshold, while
#'   the \code{"OR"} rule will preserve an edge so long as one of the two
#'   coefficients have a p-value below the supplied threshold.
#' @param gamma Numeric value of the hyperparameter for the \code{"EBIC"}
#'   criterion. Only relevant if \code{criterion = "EBIC"}. Recommended to use a
#'   value between 0 and .5, where larger values impose a larger penalty on the
#'   criterion.
#' @param nfolds Only relevant if \code{criterion = "CV"}. Determines the number
#'   of folds to use in cross-validation.
#' @param nlam if \code{method = "glinternet"}, determines the number of lambda
#'   values to evaluate in the selection path.
#' @param which.lam Character string. Only applies if \code{criterion = "CV"}.
#'   Options include \code{"min"}, which uses the lambda value that minimizes
#'   the objective function, or \code{"1se"} which uses the lambda value at 1
#'   standard error above the value that minimizes the objective function.
#' @param threshold Logical or numeric. If \code{TRUE}, then a default value of
#'   .05 will be set. Indicates whether a threshold should be placed on the
#'   models at each iteration of the sampling. A significant choice by the
#'   researcher.
#' @param bonf Logical. Determines whether to apply a bonferroni adjustment on
#'   the distribution of p-values for each coefficient.
#' @param alpha Type 1 error rate. Defaults to .05.
#' @param exogenous Logical. Indicates whether moderator variables should be
#'   treated as exogenous or not. If they are exogenous, they will not be
#'   modeled as outcomes/nodes in the network. If the number of moderators
#'   reaches \code{k - 1} or \code{k}, then \code{exogenous} will automatically
#'   be \code{FALSE}.
#' @param split If \code{sampMethod == "split"} or \code{sampMethod =
#'   "stability"} then this is a value between 0 and 1 that indicates the
#'   proportion of the sample to be used for the training set. When
#'   \code{sampMethod = "stability"} there isn't an important distinction
#'   between the labels "training" and "test", although this value will still
#'   cause the two samples to be taken of complementary size.
#' @param center Logical. Determines whether to mean-center the variables.
#' @param scale Logical. Determines whether to standardize the variables.
#' @param varSeed Numeric value providing a seed to be set at the beginning of
#'   the selection procedure. Recommended for reproducible results. Importantly,
#'   this seed will be used for the variable selection models at each iteration
#'   of the resampler. Caution this means that while each model is run with a
#'   different sample, it will always have the same seed.
#' @param seed Can be a single value, to set a seed before drawing random seeds
#'   of length \code{niter} to be used across iterations. Alternatively, one can
#'   supply a vector of seeds of length \code{niter}. It is recommended to use
#'   this argument for reproducibility over the \code{varSeed} argument.
#' @param verbose Logical. Determines whether information about the modeling
#'   progress should be displayed in the console.
#' @param lags Numeric or logical. Can only be 0, 1 or \code{TRUE} or
#'   \code{FALSE}. \code{NULL} is interpreted as \code{FALSE}. Indicates whether
#'   to fit a time-lagged network or a GGM.
#' @param binary Numeric vector indicating which columns of the data contain
#'   binary variables.
#' @param type Determines whether to use gaussian models \code{"g"} or binomial
#'   models \code{"c"}. Can also just use \code{"gaussian"} or
#'   \code{"binomial"}. Moreover, a vector of length \code{k} can be provided
#'   such that a value is given to every variable. Ultimately this is not
#'   necessary, though, as such values are automatically detected.
#' @param saveMods Logical. Indicates whether to save the models fit to the
#'   samples at each iteration or not.
#' @param saveData Logical. Determines whether to save the data from each
#'   subsample across iterations or not.
#' @param saveVars Logical. Determines whether to save the variable selection
#'   models at each iteration.
#' @param fitit Logical. Determines whether to fit the final selected model on
#'   the original sample. If \code{FALSE}, then this can still be done with
#'   \code{\link{fitNetwork}} and \code{\link{modSelect}}.
#' @param nCores Numeric value indicating the number of CPU cores to use for the
#'   resampling. If \code{TRUE}, then the
#'   \code{\link[parallel:detectCores]{parallel::detectCores}} function will be
#'   used to maximize the number of cores available.
#' @param cluster Character vector indicating which type of parallelization to
#'   use, if \code{nCores > 1}. Options include \code{"mclapply"} and
#'   \code{"SOCK"}.
#' @param block Logical or numeric. If specified, then this indicates that
#'   \code{lags != 0} or \code{lags != NULL}. If numeric, then this indicates
#'   that block bootstrapping will be used, and the value specifies the block
#'   size. If \code{TRUE} then an appropriate block size will be estimated
#'   automatically.
#' @param beepno Character string or numeric value to indicate which variable
#'   (if any) encodes the survey number within a single day. Must be used in
#'   conjunction with \code{dayno} argument.
#' @param dayno Character string or numeric value to indiciate which variable
#'   (if any) encodes the survey number within a single day. Must be used in
#'   conjunction with \code{beepno} argument.
#' @param ... Additional arguments.
#'
#' @return \code{resample} output
#' @export
#' @references Meinshausen, N., Meier, L., & Buhlmann, P. (2009). P-values for
#'   high-dimensional regression. Journal of the American Statistical
#'   Association. 104, 1671-1681.
#'
#'   Meinshausen, N., & Buhlmann, P. (2010). Stability selection. Journal of the
#'   Royal Statistical Society: Series B (Statistical Methodology). 72, 417-423
#'
#' @seealso \code{\link{plot.resample}, \link{modSelect}, \link{fitNetwork},
#'   \link{bootNet}, \link{mlGVAR}, \link{plotNet}, \link{plotCoefs},
#'   \link{plotBoot}, \link{plotPvals}, \link{plotStability}, \link{net},
#'   \link{netInts}, \link[glinternet:glinternet]{glinternet::glinternet},
#'   \link[glinternet:glinternet.cv]{glinternet::glinternet.cv},
#'   \link[glmnet:glmnet]{glmnet::glmnet},
#'   \link[glmnet:cv.glmnet]{glmnet::cv.glmnet},
#'   \link[leaps:regsubsets]{leaps::regsubsets}}
#'
#' @examples
#' \donttest{
#' fit1 <- resample(ggmDat, m = 'M', niter = 10)
#'
#' net(fit1)
#' netInts(fit1)
#'
#' plot(fit1)
#' plot(fit1, what = 'coefs')
#' plot(fit1, what = 'bootstrap', multi = TRUE)
#' plot(fit1, what = 'pvals', outcome = 2, predictor = 4)
#'
#' fit2 <- resample(gvarDat, m = 'M', niter = 10, lags = 1, sampMethod = 'stability')
#'
#' plot(fit2, what = 'stability', outcome = 3)
#' }
resample <- function(data, m = NULL, niter = 10, sampMethod = "bootstrap", criterion = "AIC",
                     method = "glmnet", rule = "OR", gamma = .5, nfolds = 10,
                     nlam = 50, which.lam = "min", threshold = FALSE, bonf = FALSE,
                     alpha = .05, exogenous = TRUE, split = .5, center = TRUE,
                     scale = FALSE, varSeed = NULL, seed = NULL, verbose = TRUE,
                     lags = NULL, binary = NULL, type = 'g', saveMods = TRUE,
                     saveData = FALSE, saveVars = FALSE, fitit = TRUE,
                     nCores = 1, cluster = 'mclapply', block = FALSE,
                     beepno = NULL, dayno = NULL, ...){
  # The 'split' argument, for sampMethod %in% c('split', 'stability')
  # Is the proportion of the data to be included in the training set
  t1 <- Sys.time() # RESAMPLE START
  if(is.null(lags) & !identical(block, FALSE)){
    message('Assuming lags = 1 since block resampling was requested')
    lags <- 1
  }
  if(!is.null(lags)){
    lags <- switch(2 - identical(as.numeric(lags), 0), NULL, 1)
    if(any(!sapply(c(beepno, dayno), is.null))){
      stopifnot(!is.null(beepno) & !is.null(dayno))
      data <- getConsec(data = data, beepno = beepno, dayno = dayno)
    }
  }
  if(any(is.na(data))){
    ww <- which(apply(data, 1, function(z) any(is.na(z))))
    stop(paste0(length(ww), ' rows contain missing values'))
  }
  consec <- switch(2 - (!is.null(lags) & 'samp_ind' %in% names(attributes(data))),
                   attr(data, 'samp_ind'), NULL)
  N <- ifelse(is.null(consec), nrow(data) - ifelse(is.null(lags), 0, lags), length(consec))
  atts <- attributes(data)
  data <- data.frame(data)
  attributes(data) <- atts
  if(!is.null(m)){if(all(m == 0)){m <- NULL}}
  method <- ifelse(!is.null(m), 'glinternet', ifelse(
    !method %in% c('glmnet', 'subset'), 'glmnet', method))
  if(isTRUE(block)){block <- floor(3.15 * N^(1/3))}
  criterion <- toupper(match.arg(tolower(criterion), c(
    "cv", "aic", "bic", "ebic", "cp", "rss", "adjr2", "rsq", "r2")))
  if(!criterion %in% c("AIC", "BIC", "EBIC", "CV")){method <- "subset"}
  method <- match.arg(method, c("glinternet", "hiernet", "subset", "glmnet"))
  if(is.character(m)){
    m <- ifelse(all(m == "all"), list(1:ncol(data)), list(which(colnames(data) %in% m)))[[1]]
  }
  if(!is.null(binary)){
    type <- rep("g", ifelse(is.null(m) | !exogenous, ncol(data), ncol(data) - 1))
    type[intersect(seq_along(type), binary)] <- "c"
  }
  if(method == "glinternet" & is.null(m)){m <- 1:ncol(data)}
  if(!is.null(m)){if(length(m) > 1){exogenous <- FALSE}}
  if(!is.character(sampMethod)){sampMethod <- c("split", "bootstrap", "stability")[sampMethod]}
  sampMethod <- match.arg(tolower(sampMethod), c("split", "bootstrap", "stability"))
  preout <- list(niter = niter, criterion = criterion, sampMethod = sampMethod,
                 method = method, moderators = m, rule = rule, alpha = alpha,
                 bonf = bonf, gamma = gamma, nfolds = nfolds, which.lam = which.lam,
                 split = split, center = center, scale = scale, exogenous = exogenous,
                 varSeed = varSeed, type = type, lags = lags, block = block)
  args <- tryCatch({list(...)}, error = function(e){list()}) ##### RESAMPLE ARGS
  if(length(args) != 0){preout <- append(preout, args)}
  if(is.null(lags)){preout$lags <- NULL} else {preout$rule <- NULL}
  if(criterion != "CV"){preout[c("nfolds", "which.lam")] <- NULL; which.lam <- "min"}
  if(nlam != 50 & method != "subset"){preout$nlam <- nlam}
  if(sampMethod == "bootstrap"){preout$split <- NULL}
  if(criterion != "EBIC"){preout$gamma <- NULL}
  if(is.null(varSeed)){preout$varSeed <- NULL}
  lam <- ifelse(grepl("min", which.lam), 1, 2)
  if(length(seed) <= 1){
    if(length(seed) == 1){set.seed(seed)}
    seeds <- sample(1:10000000, niter, replace = FALSE)
  } else {
    seeds <- seed
  }
  sampler2 <- function(N, size = N, block = FALSE, method = 'bootstrap', consec = NULL){
    allinds <- switch(2 - is.null(consec), seq_len(N), consec)
    if(identical(block, FALSE)){
      sample(allinds, size, replace = (method == 'bootstrap'))
    } else if(method == 'bootstrap'){
      stopifnot(block < N & block > 1)
      nblox <- 1
      while(N > block * nblox){nblox <- nblox + 1}
      possible <- seq_len(N - block)
      starts <- replicate(nblox, sample(possible, 1))
      sampinds <- c(sapply(starts, function(z) allinds[z:(z + block - 1)]))
      if(length(sampinds) > N){sampinds <- sampinds[1:N]}
      return(sampinds)
    } else {
      stop('Multi-split block resampling not yet implemented')
    }
  }
  sampInd <- samps <- vars <- vars1 <- fits <- train <- test <- list()
  if((exogenous | is.null(m)) & sampMethod != 'bootstrap'){
    mno <- as.numeric(!is.null(m) & !identical(m, 0))
    lno <- as.numeric(!is.null(lags) & !identical(lags, 0))
    minsize <- sampleSize(p = ncol(data) - mno, m = mno,
                          lags = lno, print = FALSE)
    minn <- ifelse(split > .5, N - floor(N * split), floor(N * split))
    if(minn < minsize & split != .5){
      split <- .5
      minn <- ifelse(split > .5, N - floor(N * split), floor(N * split))
      if(minn > minsize){warning('Split size set to .5 to produce large enough subsamples')}
    }
    if(minn < minsize){stop('Sample-split size is smaller than required')}
  }
  if(nCores > 1 | isTRUE(nCores)){
    if(isTRUE(nCores)){nCores <- parallel::detectCores()}
    if(grepl('Windows', sessionInfo()$running)){cluster <- 'SOCK'}
    if(sampMethod != 'bootstrap'){
      if(split <= 0 | split >= 1){stop('split size must be between 0 and 1')}
      n <- floor(N * split)
    } else {
      n <- N
    }
    for(i in 1:niter){ #### SAMPLING 1
      set.seed(seeds[i])
      sampInd[[i]] <- sampler2(N = N, size = n, block = block, method = sampMethod, consec = consec)
      if(sampMethod != 'bootstrap'){
        if(is.null(lags)){
          train[[i]] <- data[sampInd[[i]], ]
          test[[i]] <- data[-sampInd[[i]], ]
        } else {
          train[[i]] <- test[[i]] <- data
          attr(train[[i]], 'samp_ind') <- sampInd[[i]]
          attr(test[[i]], 'samp_ind') <- (1:N)[-sampInd[[i]]]
        }
        samps[[i]] <- vector('list', length = 2)
        samps[[i]][[1]] <- train[[i]]
        samps[[i]][[2]] <- test[[i]]
        names(samps[[i]]) <- c('train', 'test')
      } else {
        if(is.null(lags)){
          samps[[i]] <- data[sampInd[[i]], ]
        } else {
          samps[[i]] <- data
          attr(samps[[i]], 'samp_ind') <- sampInd[[i]]
        }
      }
    } ##### SAMPLING 1 END
    staticArgs <- list(m = m, criterion = criterion, center = center, method = method,
                       nfolds = nfolds, gamma = gamma, lags = lags, nlam = nlam,
                       type = type, varSeed = varSeed, exogenous = exogenous,
                       scale = scale, verbose = FALSE, saveMods = saveMods,
                       rule = rule, which.lam = which.lam)
    sampFun <- function(samp, seed, sampMethod, args){
      set.seed(seed)
      criterion <- args$criterion
      saveMods <- args$saveMods
      vars <- vars1 <- fits <- list()
      if(sampMethod != 'bootstrap'){
        vargs1 <- append(list(data = samp$train), args[intersect(names(args), formalArgs('varSelect'))])
        vars <- do.call(varSelect, vargs1)
        if(!saveMods){
          for(j in seq_len(length(vars))){
            if(criterion != 'CV'){vars[[j]]$fitobj$fit$fitted <- NA}
            if(criterion == 'CV'){
              vars[[j]]$fitobj$fitCV$glinternetFit$fitted <- NA
              vars[[j]]$fitobj$fit0 <- vars[[j]]$fitobj$fit1se <- NA
            }
          }
        }
        if(sampMethod == 'stability'){
          vargs2 <- replace(vargs1, 'data', list(data = samp$test))
          vars1 <- do.call(varSelect, vargs2)
          if(!saveMods){
            for(j in seq_len(length(vars1))){
              if(criterion != 'CV'){vars1[[j]]$fitobj$fit$fitted <- NA}
              if(criterion == 'CV'){
                vars1[[j]]$fitobj$fitCV$glinternetFit$fitted <- NA
                vars1[[j]]$fitobj$fit0 <- vars1[[j]]$fitobj$fit1se <- NA
              }
            }
          }
        } else {
          fitargs <- append(list(data = samp$test, moderators = args$m, type = vars, threshold = FALSE),
                            args[setdiff(names(args), c('saveMods', 'm', 'type'))])
          fitargs <- fitargs[intersect(names(fitargs), formalArgs('fitNetwork'))]
          fits <- do.call(fitNetwork, fitargs)
          if(!saveMods){fits$mods0 <- NULL}
        }
      } else {
        vargs1 <- append(list(data = samp), args[intersect(names(args), formalArgs('varSelect'))])
        vars <- do.call(varSelect, vargs1)
        if(!saveMods){
          for(j in seq_len(length(vars))){
            if(criterion != 'CV'){vars[[j]]$fitobj$fit$fitted <- NA}
            if(criterion == 'CV'){
              vars[[j]]$fitobj$fitCV$glinternetFit$fitted <- NA
              vars[[j]]$fitobj$fit0 <- vars[[j]]$fitobj$fit1se <- NA
            }
          }
        }
        fitargs <- append(list(data = samp, moderators = args$m, type = vars, threshold = FALSE),
                          args[setdiff(names(args), c('saveMods', 'm', 'type'))])
        fitargs <- fitargs[intersect(names(fitargs), formalArgs('fitNetwork'))]
        fits <- do.call(fitNetwork, fitargs)
        if(!saveMods){fits$mods0 <- NULL}
      }
      out <- list(vars = vars, vars1 = vars1, fits = fits)
      return(out)
    }
    obj0 <- c('sampFun', 'varSelect', 'fitNetwork', 'Matrix', 'net', 'netInts')
    obj1 <- switch(2 - is.null(lags), c('nodewise', 'modNet', 'modLL'),
                   c('lagMat', 'SURfit', 'SURnet', 'SURll', 'surEqs', 'getCoefs', 'systemfit'))
    obj2 <- switch(2 - (method == 'glinternet'), c('fitHierLASSO', ifelse(criterion == 'CV', 'glinternet.cv', 'glinternet')),
                   ifelse(method == 'glmnet', ifelse(criterion == 'CV', 'cv.glmnet', 'glmnet'), 'regsubsets'))
    if(method == 'glmnet'){obj2 <- c(obj2, 'lassoSelect')}
    objects <- c(obj0, obj1, obj2)
    if(tolower(cluster) != 'mclapply'){
      cluster <- match.arg(toupper(cluster), c('SOCK', 'FORK'))
      cl <- parallel::makeCluster(nCores, type = cluster)
      if(cluster == 'SOCK'){parallel::clusterExport(cl, objects, envir = environment())}
    } else {
      cl <- nCores
    }
    if(verbose){
      pbapply::pboptions(type = 'timer', char = '-')
      cl_out <- pbapply::pblapply(1:niter, function(i){
        sampFun(samp = samps[[i]], seed = seeds[i],
                sampMethod = sampMethod, args = staticArgs)
      }, cl = cl)
    } else if(tolower(cluster) == 'mclapply'){
      cl_out <- parallel::mclapply(1:niter, function(i){
        sampFun(samp = samps[[i]], seed = seeds[i],
                sampMethod = sampMethod, args = staticArgs)
      }, mc.cores = nCores)
    } else {
      cl_out <- parallel::parLapply(cl, 1:niter, function(i){
        sampFun(samp = samps[[i]], seed = seeds[i],
                sampMethod = sampMethod, args = staticArgs)
      })
    }
    if(tolower(cluster) != 'mclapply'){parallel::stopCluster(cl)}
    vars <- lapply(cl_out, '[[', 1)
    vars1 <- lapply(cl_out, '[[', 2)
    fits <- lapply(cl_out, '[[', 3)
    rm(cl_out, staticArgs, objects, obj1, obj2, obj0, sampFun, cl)
  } else {
    if(sampMethod != "bootstrap"){
      if(split <= 0 | split >= 1){stop("split size must be between 0 and 1")}
      train <- list(); test <- list()
      n <- floor(N * split)
      for(i in 1:niter){ ##### SAMPLING 2
        set.seed(seeds[i]) # RESAMPLE SPLIT
        sampInd[[i]] <- sampler2(N = N, size = n, block = block, method = sampMethod, consec = consec)
        if(is.null(lags)){
          train[[i]] <- data[sampInd[[i]], ]
          test[[i]] <- data[-sampInd[[i]], ]
        } else {
          train[[i]] <- test[[i]] <- data
          attr(train[[i]], "samp_ind") <- sampInd[[i]]
          attr(test[[i]], "samp_ind") <- (1:N)[-sampInd[[i]]]
        }
        samps[[i]] <- vector("list", length = 2)
        samps[[i]][[1]] <- train[[i]]
        samps[[i]][[2]] <- test[[i]]
        names(samps[[i]]) <- c("train", "test")
        if(verbose){ ##### SAMPLING 2 END
          if(i == 1){cat("\n")}
          t2 <- Sys.time()
          cat("************ Sample Split: ", i, "/", niter, " ************\n", sep = "")
        }
        set.seed(seeds[i])
        vars[[i]] <- varSelect(data = train[[i]], m = m, criterion = criterion, center = center,
                               method = method, nfolds = nfolds, gamma = gamma, lags = lags,
                               nlam = nlam, type = type, varSeed = varSeed, exogenous = exogenous,
                               scale = scale, verbose = ifelse(sampMethod == "stability", ifelse(
                                 verbose, "pbar", FALSE), verbose))
        if(!saveMods){
          for(j in seq_len(length(vars[[i]]))){
            if(criterion != "CV"){vars[[i]][[j]]$fitobj$fit$fitted <- NA}
            if(criterion == "CV"){
              vars[[i]][[j]]$fitobj$fitCV$glinternetFit$fitted <- NA
              vars[[i]][[j]]$fitobj$fit0 <- vars[[i]][[j]]$fitobj$fit1se <- NA
            }
          }
        } # RESAMPLE STABILITY
        if(sampMethod == "stability"){
          vars1[[i]] <- varSelect(data = test[[i]], m = m, criterion = criterion, center = center,
                                  method = method, nfolds = nfolds, gamma = gamma, lags = lags,
                                  nlam = nlam, type = type, varSeed = varSeed, exogenous = exogenous,
                                  scale = scale, verbose = ifelse(verbose, "pbar", FALSE))
          if(!saveMods){
            for(j in seq_len(length(vars1[[i]]))){
              if(criterion != "CV"){vars1[[i]][[j]]$fitobj$fit$fitted <- NA}
              if(criterion == "CV"){
                vars1[[i]][[j]]$fitobj$fitCV$glinternetFit$fitted <- NA
                vars1[[i]][[j]]$fitobj$fit0 <- vars1[[i]][[j]]$fitobj$fit1se <- NA
              }
            }
          }
          if(verbose & !method %in% c("subset", "glmnet")){
            t3 <- Sys.time() - t2
            cat(paste0("Time: ", round(t3, 2), " ", attr(t3, "units")), "\n\n")
          } # RESAMPLE SPLIT FIT
        } else {
          fits[[i]] <- fitNetwork(data = test[[i]], moderators = m, type = vars[[i]],
                                  threshold = FALSE, which.lam = which.lam, lags = lags,
                                  gamma = gamma, center = center, rule = rule, scale = scale,
                                  exogenous = exogenous, verbose = FALSE, ...)
          if(!saveMods){fits[[i]]$mods0 <- NULL}
        }
      }
    }
    if(sampMethod == "bootstrap"){
      for(i in 1:niter){ ##### SAMPLING 3
        set.seed(seeds[i]) # RESAMPLE BOOT
        sampInd[[i]] <- sampler2(N = N, block = block, consec = consec) # Currently only for sampMethod == 'bootstrap'
        if(is.null(lags)){
          samps[[i]] <- data[sampInd[[i]], ]
        } else {
          samps[[i]] <- data
          attr(samps[[i]], "samp_ind") <- sampInd[[i]]
        }
        if(verbose == TRUE){ ##### SAMPLING 3 END
          if(i == 1){cat("\n")}
          cat("************* Bootstrap: ", i, "/", niter, " **************\n", sep = "")
        }
        set.seed(seeds[i])
        vars[[i]] <- varSelect(data = samps[[i]], m = m, criterion = criterion,
                               center = center, method = method, gamma = gamma,
                               lags = lags, type = type, varSeed = varSeed,
                               scale = scale, exogenous = exogenous,
                               verbose = verbose)
        if(!saveMods){
          for(j in seq_len(length(vars[[i]]))){
            if(criterion != "CV"){vars[[i]][[j]]$fitobj$fit$fitted <- NA}
            if(criterion == "CV"){
              vars[[i]][[j]]$fitobj$fitCV$glinternetFit$fitted <- NA
              vars[[i]][[j]]$fitobj$fit0 <- vars[[i]][[j]]$fitobj$fit1se <- NA
            }
          }
        }
        fits[[i]] <- fitNetwork(data = samps[[i]], moderators = m, type = vars[[i]],
                                threshold = FALSE, which.lam = which.lam, lags = lags,
                                gamma = gamma, center = center, rule = rule, scale = scale,
                                exogenous = exogenous, verbose = FALSE, ...)
        if(!saveMods){fits[[i]]$mods0 <- NULL}
      }
    }
  }
  fit0 <- fitNetwork(data = data, moderators = m, type = type, rule = rule, scale = scale,
                     threshold = FALSE, gamma = gamma, center = center, lags = lags,
                     exogenous = exogenous, verbose = FALSE, ...)
  if(!is.null(lags)){
    if(sampMethod != "stability"){
      fits <- lapply(fits, function(fi){
        append(fi$SURnet[setdiff(names(fi$SURnet), 'data')], append(
          switch(2 - ('SURll' %in% names(fi)), list(SURll = fi$SURll), list()),
          list(data = fi$SURnet$data, fitobj = fi$SURfit$eq)))
      })
    }
    fit0 <- append(fit0$SURnet[setdiff(names(fit0$SURnet), 'data')], append(
      switch(2 - ('SURll' %in% names(fit0)), list(SURll = fit0$SURll), list()),
      list(data = fit0$SURnet$data, fitobj = fit0$SURfit$eq)))
  }
  p <- length(fit0$mods)
  vs <- names(fit0$mods)
  allNames <- lapply(fit0$mods, function(z) rownames(z$model)[-1])
  ### RESAMPLE VARMODS
  varMods0 <- lapply(vars, function(z) lapply(z, function(zz) zz[[lam]]))
  if(sampMethod == "stability"){
    if(!is.numeric(threshold)){threshold <- .6}
    varMods1 <- lapply(vars1, function(z) lapply(z, '[[', lam))
    simultvars <- lapply(seq_len(niter), function(z){
      lapply(seq_len(p), function(zz){
        varSamp1 <- !is.na(match(allNames[[zz]], varMods0[[z]][[zz]]))
        varSamp2 <- !is.na(match(allNames[[zz]], varMods1[[z]][[zz]]))
        return(cbind(varSamp1, varSamp2, varSamp1 & varSamp2))
      })
    })
    varFreqs <- suppressMessages(suppressWarnings(lapply(seq_len(p), function(z){
      simultSel <- t(sapply(seq_len(niter), function(zz) simultvars[[zz]][[z]][, 3]))
      s1sel <- t(sapply(seq_len(niter), function(zz) simultvars[[zz]][[z]][, 1]))
      s2sel <- t(sapply(seq_len(niter), function(zz) simultvars[[zz]][[z]][, 2]))
      colnames(simultSel) <- colnames(s1sel) <- colnames(s2sel) <- allNames[[z]]
      freqs <- colMeans(simultSel)
      s1freqs <- colMeans(s1sel)
      s2freqs <- colMeans(s2sel)
      coefs1 <- summary(fit0$fitobj[[z]])$coefficients[-1, , drop = FALSE]
      coefs2 <- confint(fit0$fitobj[[z]], level = 1 - alpha)[-1, , drop = FALSE]
      outFreqs <- data.frame(predictor = factor(allNames[[z]]), select = freqs >= threshold,
                             lower = coefs2[, 1], b = coefs1[, 1], upper = coefs2[, 2],
                             Pvalue = coefs1[, 4], freq = freqs, split1 = s1freqs, split2 = s2freqs)
      rownames(outFreqs) <- 1:nrow(outFreqs)
      return(outFreqs)
    })))
    preout$nlam <- nlam
    attributes(varFreqs) <- list(
      names = vs, type = type, criterion = criterion, rule = rule, gamma = gamma,
      center = center, scale = scale, exogenous = exogenous, moderators = m,
      threshold = FALSE, lags = lags, method = unique(sapply(vars, attr, "method")))
    if(is.null(lags)){attr(varFreqs, "lags") <- NULL}
    split1 <- lapply(seq_len(niter), function(z){list(vars = varMods0[[z]], varMods = vars[[z]])})
    split2 <- lapply(seq_len(niter), function(z){list(vars = varMods1[[z]], varMods = vars1[[z]])})
    names(split1) <- names(split2) <- paste0("iter", 1:niter)
    out <- list(call = preout, freqs = varFreqs, split1 = split1, split2 = split2)
    if(!(method == "subset" & criterion == "CV")){
      out$stability <- suppressWarnings(tryCatch({
        append(stability(out), list(seeds = seeds))}, error = function(e){
          list(split1 = split1, split2 = split2, seeds = seeds)}))
    } else {
      out$seeds <- seeds
    }
    if(!saveVars & criterion != "CV"){out <- out[-c(3:4)]}
    if(fitit){
      tryfit <- tryCatch({modSelect(obj = out, data = data, fit = TRUE, saveMods = saveMods)},
                         error = function(e){TRUE})
      if(!isTRUE(tryfit)){out$fit0 <- tryfit}
    }
    out$results <- structure(cbind.data.frame(
      outcome = rep(gsub('[.]y$', '', vs), sapply(out$freqs, nrow)),
      do.call(rbind, out$freqs)), row.names = 1:sum(sapply(out$freqs, nrow)))
    if(saveData){out$data <- data}
    if(verbose){cat("\n"); print(Sys.time() - t1)}
    attr(out, 'resample') <- TRUE
    class(out) <- c('list', 'resample')
    return(out)
  } else {
    mod0Freqs <- lapply(1:p, function(z){
      data.frame(table(unlist(lapply(varMods0, function(zz) zz[[z]]))))})
    mod0Freqs2 <- lapply(1:p, function(z) mod0Freqs[[z]][, 2])
    for(i in 1:p){
      names(mod0Freqs2[[i]]) <- as.character(mod0Freqs[[i]][, 1])
      mod0Freqs2[[i]] <- mod0Freqs2[[i]][match(allNames[[i]], names(mod0Freqs2[[i]]))]
      names(mod0Freqs2[[i]]) <- allNames[[i]]
      if(any(is.na(mod0Freqs2[[i]]))){mod0Freqs2[[i]][is.na(mod0Freqs2[[i]])] <- 0}
      mod0Freqs2[[i]] <- mod0Freqs2[[i]]/niter
    }
    mod0Freqs <- lapply(mod0Freqs2, as.matrix, ncol = 1)
    names(mod0Freqs) <- names(fit0$mods)
    mod0sizes <- sapply(varMods0, function(z) sapply(z, length))
    for(i in 1:niter){
      if(any(mod0sizes[, i] != max(mod0sizes[, i]))){
        ms0 <- which(mod0sizes[, i] != max(mod0sizes[, i]))
        for(j in seq_along(ms0)){
          varMods0[[i]][[ms0[j]]] <- c(
            varMods0[[i]][[ms0[j]]],
            rep("", max(mod0sizes[, i]) - length(varMods0[[i]][[ms0[j]]])))
        }
      }
    }
    mod0coefs <- finalCoefs0 <- list()
    if(nCores > 1 & FALSE){
      summaryFun <- function(i, fits, p, alpha, allNames, bonf, mod0sizes){
        mod0coefs <- suppressWarnings(suppressMessages(
          lapply(fits[[i]]$fitobj, function(z){
            coefs1 <- t(summary(z)$coefficients[-1, , drop = FALSE])
            if(length(coefs1) != 0){
              coefs2 <- t(confint(z, level = 1 - alpha)[-1, , drop = FALSE])
            } else {
              vsx <- allNames[[as.character(z$terms[[2]])]]
              coefs1 <- matrix(NA, nrow = 4, ncol = length(vsx))
              coefs2 <- matrix(NA, nrow = 2, ncol = length(vsx))
              colnames(coefs1) <- colnames(coefs2) <- vsx
            }
            return(rbind(coefs1, coefs2))
          })))
        for(j in 1:p){
          rownames(mod0coefs[[j]]) <- c("b", "se", "t", "P", "lower", "upper")
          mod0coefs[[j]] <- mod0coefs[[j]][, match(allNames[[j]], colnames(mod0coefs[[j]]))]
          colnames(mod0coefs[[j]]) <- allNames[[j]]
          if(bonf){mod0coefs[[j]]["P", ] <- pmin(mod0coefs[[j]]["P", ] * mod0sizes[j, i], 1)}
        }
        return(mod0coefs)
      }
      cl <- parallel::makeCluster(nCores, type = 'SOCK')
      parallel::clusterExport(cl, 'summaryFun', envir = environment())
      mod0coefs <- parallel::parLapply(cl, 1:niter, summaryFun, fits = fits, p = p,
                                       alpha = alpha, allNames = allNames,
                                       bonf = bonf, mod0sizes = mod0sizes)
      parallel::stopCluster(cl)
    } else {
      if(verbose){
        npb <- seq(43, 53, by = 2)[sum(niter >= c(10, 100, 1000, 10000, 100000)) + 1]
        cat("\n"); cat(paste0(rep("#", npb), collapse = ""), "\n")
        cat(paste0("Estimating ", (1 - alpha) * 100, "% CIs\n"))
        pb <- txtProgressBar(min = 0, max = niter, style = 1, char = "-", width = npb)
      }
      for(i in 1:niter){
        mod0coefs[[i]] <- suppressWarnings(suppressMessages(
          lapply(fits[[i]]$fitobj, function(z){
            coefs1 <- t(summary(z)$coefficients[-1, , drop = FALSE])
            if(length(coefs1) != 0){
              coefs2 <- t(confint(z, level = 1 - alpha)[-1, , drop = FALSE])
            } else {
              vsx <- allNames[[as.character(z$terms[[2]])]]
              coefs1 <- matrix(NA, nrow = 4, ncol = length(vsx))
              coefs2 <- matrix(NA, nrow = 2, ncol = length(vsx))
              colnames(coefs1) <- colnames(coefs2) <- vsx
            }
            return(rbind(coefs1, coefs2))
          })))
        for(j in 1:p){
          rownames(mod0coefs[[i]][[j]]) <- c("b", "se", "t", "P", "lower", "upper")
          mod0coefs[[i]][[j]] <- mod0coefs[[i]][[j]][, match(allNames[[j]], colnames(mod0coefs[[i]][[j]]))]
          colnames(mod0coefs[[i]][[j]]) <- allNames[[j]]
          if(bonf){mod0coefs[[i]][[j]]["P", ] <- pmin(mod0coefs[[i]][[j]]["P", ] * mod0sizes[j, i], 1)}
        }
        if(verbose){setTxtProgressBar(pb, i); if(i == niter){close(pb)}}
      }
    }
    n <- ifelse(sampMethod == "bootstrap", 0, n)
    for(j in 1:p){
      finalCoefs0[[j]] <- vector("list", length = 6)
      finalCoefs0[[j]][[1]] <- t(sapply(mod0coefs, function(z) z[[j]]["b", ]))
      finalCoefs0[[j]][[2]] <- t(sapply(mod0coefs, function(z) z[[j]]["se", ]))
      if(any(is.na(finalCoefs0[[j]][[2]]))){finalCoefs0[[j]][[2]][is.na(finalCoefs0[[j]][[2]])] <- Inf}
      finalCoefs0[[j]][[3]] <- t(sapply(mod0coefs, function(z) z[[j]]["P", ]))
      if(any(is.na(finalCoefs0[[j]][[3]]))){finalCoefs0[[j]][[3]][is.na(finalCoefs0[[j]][[3]])] <- 1}
      finalCoefs0[[j]][[4]] <- t(sapply(mod0coefs, function(z) z[[j]]["lower", ]))
      if(any(is.na(finalCoefs0[[j]][[4]]))){finalCoefs0[[j]][[4]][is.na(finalCoefs0[[j]][[4]])] <- -Inf}
      finalCoefs0[[j]][[5]] <- t(sapply(mod0coefs, function(z) z[[j]]["upper", ]))
      if(any(is.na(finalCoefs0[[j]][[5]]))){finalCoefs0[[j]][[5]][is.na(finalCoefs0[[j]][[5]])] <- Inf}
      finalCoefs0[[j]][[6]] <- unname((N - n - mod0sizes - 1)[j, ])
      names(finalCoefs0[[j]]) <- c("b", "se", "P", "lower", "upper", "df.res")
    }
    names(finalCoefs0) <- names(fit0$mods)
    gammas <- seq(ceiling(alpha * niter)/niter, 1 - 1/niter, by = 1/niter)
    penalty <- ifelse(length(gammas) > 1, (1 - log(min(gammas))), 1)
    adjCIs0 <- list()
    for(i in 1:p){
      k <- ncol(finalCoefs0[[i]]$P)
      pvals0 <- gamma0 <- c()
      for(j in 1:k){
        quantGam0 <- quantile(finalCoefs0[[i]][["P"]][, j], gammas)/gammas
        pvals0[j] <- pmin((min(quantGam0) * penalty), 1)
        gamma0[j] <- gammas[which.min(quantGam0)]
      }
      sInd <- rep(1:k, each = niter)
      adjCIs0[[i]] <- suppressMessages(
        t(mapply(adjustedCI, lci = split(finalCoefs0[[i]][["lower"]], sInd),
                 rci = split(finalCoefs0[[i]][["upper"]], sInd),
                 centers = split(finalCoefs0[[i]][["b"]], sInd),
                 ses = split(finalCoefs0[[i]][["se"]], sInd),
                 df.res = list(df.res = finalCoefs0[[i]][["df.res"]]),
                 gamma.min = min(gammas), ci.level = (1 - alpha), var = 1:k))
      )
      colnames(adjCIs0[[i]]) <- c("lower", "upper")
      rownames(adjCIs0[[i]]) <- rownames(mod0Freqs[[i]])
      err1 <- unname(apply(adjCIs0[[i]], 1, function(z) all(z == 0)))
      if(any(err1)){
        message(paste0("Errors in ", sum(err1), " adjusted CI", ifelse(sum(err1) > 1, "s", ""), "\n"))
        err2 <- which(err1)
        for(e in 1:length(err2)){
          centers <- finalCoefs0[[i]]$b[, err2[e]]
          centers <- centers[!is.infinite(centers)]
          margerrs <- sd(centers, na.rm = TRUE) * qt(c(alpha/2, 1 - alpha/2), length(centers) - 1)
          adjCIs0[[i]][err2[e], ] <- mean(centers, na.rm = TRUE) + margerrs
        }
      }
      adjCIs0[[i]] <- data.frame(adjCIs0[[i]])
      if(!all(is.na(adjCIs0[[i]]))){
        if(any(adjCIs0[[i]][!is.na(adjCIs0[[i]]$lower), "lower"] == -Inf)){
          adjCIs0[[i]][which(adjCIs0[[i]]$lower == -Inf), ] <- NA
        }
      }
      selectp <- ifelse(is.na(pvals0), FALSE, ifelse(pvals0 <= alpha, TRUE, FALSE))
      select_ci <- ifelse(
        is.na(adjCIs0[[i]]$lower), FALSE, ifelse(
          adjCIs0[[i]]$lower <= 0 & adjCIs0[[i]]$upper >= 0, FALSE, TRUE))
      adjCIs0[[i]] <- data.frame(predictor = factor(rownames(adjCIs0[[i]])),
                                 select_ci = select_ci, lower = adjCIs0[[i]]$lower,
                                 b = rowMeans(adjCIs0[[i]]), upper = adjCIs0[[i]]$upper,
                                 Pvalue = pvals0, select = selectp, gamma = gamma0,
                                 freq = mod0Freqs[[i]][, 1])
      rownames(adjCIs0[[i]]) <- 1:nrow(adjCIs0[[i]])
      if(any(err1)){attr(adjCIs0[[i]], "err") <- as.character(adjCIs0[[i]]$predictor)[err2]}
    }
    attributes(adjCIs0) <- list(
      names = vs, type = type, criterion = criterion, rule = rule, gamma = gamma,
      center = center, scale = scale, exogenous = exogenous, moderators = m,
      threshold = FALSE, lags = lags, method = unique(sapply(vars, attr, "method")))
    if(is.null(lags)){attr(adjCIs0, "lags") <- NULL}
    if("err" %in% unlist(lapply(adjCIs0, function(z) names(attributes(z))))){
      attr(adjCIs0, "err") <- names(which(sapply(lapply(adjCIs0, function(z){
        names(attributes(z))}), function(zz) "err" %in% zz)))
    }
    for(i in 1:p){
      finalCoefs0[[i]]$se[is.infinite(finalCoefs0[[i]]$se)] <- NA
      finalCoefs0[[i]]$lower[is.infinite(finalCoefs0[[i]]$lower)] <- NA
      finalCoefs0[[i]]$upper[is.infinite(finalCoefs0[[i]]$upper)] <- NA
    }
  }
  sci <- sapply(adjCIs0, function(z) !all(z$select_ci == z$select))
  if(any(sci)){
    message(paste0('CIs different than p-values for:\n', paste(
      names(sci)[sci], collapse = '\n')))
  }
  varMods0 <- lapply(varMods0, function(z) do.call(cbind.data.frame, z))
  boots <- lapply(1:niter, function(z){
    p1 <- list(vars = varMods0[[z]])
    if(saveMods){
      fits[[z]]$fitobj <- NULL
      if(is.null(lags)){fits[[z]]$mods0$models <- NULL}
      if(!is.null(lags)){attr(fits[[z]], 'SURnet') <- TRUE}
      if(!saveData){fits[[z]]$data <- NULL}
      p1$fit <- fits[[z]]
    }
    if(saveVars){p1$varMods <- vars[[z]]}
    return(append(p1, list(samp_inds = sampInd[[z]])))
  })
  names(boots) <- paste0("iter", 1:niter)
  out <- list(call = preout, adjCIs = adjCIs0)
  out$samples <- list(coefs = finalCoefs0, iters = boots, seeds = seeds)
  if(fitit){
    tryfit <- tryCatch({modSelect(obj = out, data = data, fit = TRUE, saveMods = saveMods)},
                       error = function(e){TRUE})
    if(!isTRUE(tryfit)){out$fit0 <- tryfit}
  }
  out$results <- structure(cbind.data.frame(
    outcome = rep(gsub('[.]y$', '', vs), sapply(out$adjCIs, nrow)),
    do.call(rbind, out$adjCIs)), row.names = 1:sum(sapply(out$adjCIs, nrow)))
  out$data <- data
  if(verbose){
    if(nCores == 1){cat("\n")}
    print(Sys.time() - t1)
  }
  attr(out, "time") <- Sys.time() - t1
  attr(out, 'resample') <- TRUE
  class(out) <- c('list', 'resample')
  out
}

#' Select a model based on output from \code{\link{resample}}
#'
#' Creates the necessary input for fitNetwork when selecting variables based on
#' the \code{\link{resample}} function. The purpose of making this function
#' available to the user is to that different decisions can be made about how
#' exactly to use the \code{\link{resample}} output to select a model, as
#' sometimes there is more than one option for choosing a final model.
#'
#' @param obj \code{\link{resample}} output
#' @param data The dataframe used to create the \code{resample} object.
#'   Necessary if \code{ascall = TRUE} or \code{fit = TRUE}.
#' @param fit Logical. Determines whether to fit the selected model to the data
#'   or just return the model specifications. Must supply a dataset in the
#'   \code{data} argument as well.
#' @param select Character string, referring to which variable of the output
#'   should be used as the basis for selecting variables. If the resampling
#'   method was either \code{"bootstrap"} or \code{"split"}, then setting
#'   \code{select = "select"} will select variables based on the aggregated
#'   p-values being below a pre-specified threshold. Setting \code{select =
#'   "select_ci"}, however, will use the adjusted confidence intervals rather
#'   than p-values to select variables. Alternatively, if \code{select = "freq"}
#'   then the \code{thresh} argument can be used to indicate the minimum
#'   selection frequency across iterations. In this case, variables are selected
#'   based on how frequently they were selected in the resampling procedure.
#'   This also works if \code{select} is simply set a numeric value (this value
#'   will serve as the value for \code{thresh}).
#'
#'   When the resampling method was \code{"stability"}, the default option of
#'   \code{select = "select"} chooses variables based on the original threshold
#'   provided to the \code{\link{resample}} function, and relies on the
#'   simultaneous selection proportion (the \code{"freq"} column in the
#'   \code{"results"} element). Alternatively, if \code{select} is a numeric
#'   value, or a value for \code{thresh} is provided, that new frequency
#'   selection threshold will determine the choice of variables. Alternatively,
#'   one can specify \code{select = "split1"} or \code{select = "split2"} to
#'   base the threshold on the selection frequency in one of the two splits
#'   rather than on the simultaneous selection frequency which is likely to be
#'   the most conservative.
#'
#'   For all types of \code{resample} objects, when \code{select = "Pvalue"}
#'   then \code{thresh} can be set to a numeric value in order to select
#'   variables based on aggregated p-values. For the \code{"bootstrapping"} and
#'   \code{"split"} methods this allows one to override the original threshold
#'   (set as part of \code{\link{resample}}) if desired.
#' @param thresh Numeric value. If \code{select = "Pvalue"}, then this value
#'   will be the p-value threshold. Otherwise, this value will determine the
#'   minimum frequency selection threshold.
#' @param ascall Logical. Determines whether to return a list with arguments
#'   necessary for fitting the model with \code{do.call} to
#'   \code{\link{fitNetwork}}. Only possible if a dataset is supplied.
#' @param type Should just leave as-is. Automatically taken from the
#'   \code{resample} object.
#' @param ... Additional arguments.
#'
#' @return A call ready for \code{\link{fitNetwork}}, a fitted network model, or
#'   a list of selected variables for each node along with relevant attributes.
#'   Essentially, the output is either the selected model itself or a list of
#'   the necessary parameters to fit it.
#' @export
#'
#' @seealso \code{\link{resample}}
#'
#' @examples
#' \donttest{
#' res1 <- resample(ggmDat, m = 'M', niter = 10)
#' mods1 <- modSelect(res1)
#' fit1 <- fitNetwork(ggmDat, morderators = 'M', type = mods1)
#'
#' res2 <- resample(ggmDat, m = 'M', sampMethod = 'stability')
#' fit2 <- modSelect(res2, data = ggmDat, fit = TRUE, thresh = .7)
#' }
modSelect <- function(obj, data = NULL, fit = FALSE, select = "select",
                      thresh = NULL, ascall = TRUE, type = "gaussian", ...){
  if(is.null(data)){ascall <- FALSE}
  cis <- c("adjCIs", "freqs")[which(c("adjCIs", "freqs") %in% names(obj))]
  if(length(cis) != 0){obj <- obj[[cis]]}
  allNames <- lapply(lapply(obj, '[[', "predictor"), as.character)
  vs <- names(allNames) <- names(obj)
  if(is.numeric(select)){
    for(i in seq_along(obj)){obj[[i]]$select <- obj[[i]]$freq >= select}
    select <- "select"
  } else if(!is.null(thresh)){
    if(grepl('select|ci', select)){select <- 'freq'}
    select <- match.arg(select, c("freq", "split1", "split2", "Pvalue"))
    for(i in seq_along(obj)){
      if(select != "Pvalue"){
        obj[[i]]$select <- obj[[i]][[select]] >= thresh
      } else {
        obj[[i]]$select <- obj[[i]][[select]] <= thresh
      }
    }
    select <- "select"
  } else if(select == "ci"){
    select <- "select_ci"
  }
  if(select == 'select_ci'){
    sel <- sapply(obj, function(z) all(z$select_ci == z$select))
    if(all(sel)){message('CIs lead to same decisions as p-values')}
  }
  varMods <- lapply(seq_along(obj), function(z){
    list(mod0 = as.character(obj[[z]][which(obj[[z]][[select]]), "predictor"]))
  })
  if(any(sapply(lapply(varMods, '[[', "mod0"), length) == 0)){
    v0 <- which(sapply(lapply(varMods, '[[', "mod0"), length) == 0)
    for(jj in seq_along(v0)){varMods[[v0[jj]]]$mod0 <- "1"}
  }
  attributes(varMods) <- attributes(obj)
  type <- ifelse("type" %in% names(attributes(obj)),
                 list(attr(obj, "type")), list(type))[[1]]
  for(i in seq_along(obj)){
    if(any(grepl(":", varMods[[i]]$mod0))){
      ints <- varMods[[i]]$mod0[grepl(":", varMods[[i]]$mod0)]
      mains <- unique(unlist(strsplit(ints, ":")))
      mains <- mains[order(match(mains, vs))]
      varMods[[i]]$mod0 <- union(varMods[[i]]$mod0, union(mains, ints))
      varMods[[i]]$mod0 <- varMods[[i]]$mod0[match(allNames[[i]], varMods[[i]]$mod0)]
      varMods[[i]]$mod0 <- varMods[[i]]$mod0[!is.na(varMods[[i]]$mod0)]
    }
    attr(varMods[[i]], "family") <- ifelse(
      length(type) == 1, ifelse(type %in% c("g", "gaussian"), "g", "c"),
      ifelse(type[i] %in% c("g", "gaussian"), "g", "c")
    )
  }
  if(!is.null(data)){
    if(!is(data, "data.frame")){
      data <- data.frame(data, check.names = FALSE)
    }
    args <- list(...)
    atts <- attributes(obj)
    atts[c("names", "criterion", "method")] <- NULL
    atts$data <- data
    atts$type <- varMods
    if("lags" %in% names(attributes(atts$type))){
      if(attr(atts$type, "lags") == 1 & "moderators" %in% names(attributes(atts$type))){
        m <- attr(atts$type, "moderators")
        attr(atts$type, "moderators") <- colnames(data)[m]
        attr(varMods, "moderators") <- colnames(data)[m]
      }
    }
    if("err" %in% names(atts)){atts$err <- NULL}
    dupArgs <- intersect(names(args), names(atts))
    newArgs <- setdiff(names(args), names(atts))
    attr(varMods, "call") <- atts <- append(
      replace(atts, dupArgs, args[dupArgs]), args[newArgs])
    if(fit){return(do.call(fitNetwork, atts))}
  }
  if(ascall){return(attr(varMods, "call"))} else {return(varMods)}
}

#' Plot method for output of resample function
#'
#' Allows one to plot results from the \code{\link{resample}} function based on
#' a few different options.
#'
#' For the \code{what} argument, the options correspond with calls to the
#' following functions: \itemize{ \item{\code{"network": \link{plotNet}}}
#' \item{\code{"bootstrap": \link{plotBoot}}} \item{\code{"coefs":
#' \link{plotCoefs}}} \item{\code{"stability": \link{plotStability}}}
#' \item{\code{"pvals": \link{plotPvals}}} }
#'
#' \code{"bootstrap"} and \code{"pvals"} only available for bootstrapped and
#' multi-sample split resampling. \code{"stability"} only available for
#' stability selection.
#'
#' @param x Output from the \code{\link{resample}} function.
#' @param what Can be one of three options for all \code{\link{resample}}
#'   outputs: \code{what = "network"} will plot the final network model selected
#'   from resampling. \code{what = "bootstrap"} will run \code{\link{bootNet}}
#'   based on the final model to create bootstrapped estimates of confidence
#'   bands around each edge estimate. \code{what = "coefs"} will plot the
#'   confidence intervals based on the model parameters in the final network.
#'   Additionally, if the object was fit with \code{sampMethod = "stability"}, a
#'   stability plot can be created with the \code{"stability"} option.
#'   Otherwise, if \code{sampMethod = "bootstrap"} or \code{sampMethod =
#'   "split"}, a plot of the empirical distribution function of p-values can be
#'   displayed with the \code{"pvals"} option.
#' @param ... Additional arguments.
#'
#' @return Plots various aspects of output from the \code{\link{resample}}
#'   function.
#' @export
#'
#' @examples
#' \donttest{
#' fit1 <- resample(ggmDat, m = 'M', niter = 10)
#'
#' net(fit1)
#' netInts(fit1)
#'
#' plot(fit1)
#' plot(fit1, what = 'coefs')
#' plot(fit1, what = 'bootstrap', multi = TRUE)
#' plot(fit1, what = 'pvals', outcome = 2, predictor = 4)
#'
#' fit2 <- resample(gvarDat, m = 'M', niter = 10, lags = 1, sampMethod = 'stability')
#'
#' plot(fit2, what = 'stability', outcome = 3)
#' }
plot.resample <- function(x, what = 'network', ...){
  args <- tryCatch({list(...)}, error = function(e){list()})
  if(isTRUE(what) | is(what, 'numeric')){
    args$threshold <- ifelse(!'threshold' %in% names(args), what, args$threshold)
    what <- 'network'
  }
  what <- match.arg(tolower(what), c('network', 'bootstrap', 'coefs', 'stability', 'pvals'))
  if(what == 'stability'){stopifnot('stability' %in% names(x))}
  if(what == 'pvals'){stopifnot(x$call$sampMethod != 'stability')}
  FUN <- switch(what, network = plotNet, bootstrap = plotBoot, coefs = plotCoefs,
                stability = plotStability, pvals = plotPvals)
  if(what == 'bootstrap' & !is.null(x$call$moderators)){
    if('fit0' %in% names(x)){
      if(!any(grepl(':', unlist(x$fit0$call$type)))){
        x$call <- replace(x$call, 'moderators', list(moderators = NULL))
      }
    }
  }
  args$XX <- x
  names(args)[which(names(args) == 'XX')] <- ifelse(what == 'coefs', 'fit', 'x')
  #names(args)[which(names(args) == 'XX')] <- switch(what, network = 'x', bootstrap = 'x', coefs = 'fit')
  do.call(FUN, args)
}

##### adjustedCI: uses the union-bound approach for obtaining adjCIs. Adapted from hdi package
adjustedCI <- function(lci, rci, centers, ses, df.res, gamma.min, ci.level, var,
                       multi.corr = FALSE, verbose = FALSE, s0 = list(s0 = NA)){
  findInsideGamma <- function(low, high, ci.info, verbose){
    range.length <- 10
    test.range <- seq(low, high, length.out = range.length)
    cover <- mapply(doesItCoverGamma, beta.j = test.range, ci.info = list(ci.info = ci.info))
    while(!any(cover)){
      range.length <- 10 * range.length
      test.range <- seq(low, high, length.out = range.length)
      cover <- mapply(doesItCoverGamma, beta.j = test.range, ci.info = list(ci.info = ci.info))
      if(range.length > 10^3){
        message("FOUND NO INSIDE POINT")
        message("number of splits")
        message(length(ci.info$centers))
        message("centers")
        message(ci.info$centers)
        message("ses")
        message(ci.info$ses)
        message("df.res")
        message(ci.info$df.res)
        return(NULL)
        #stop("couldn't find an inside point between low and high. The confidence interval doesn't exist!")
      }
    }
    if(verbose){cat("Found an inside point at granularity of", range.length, "\n")}
    min(test.range[cover])
  }
  doesItCoverGamma <- function(beta.j, ci.info){
    if(missing(ci.info)){stop("ci.info is missing to the function doesItCoverGamma")}
    centers <- ci.info$centers
    ci.lengths <- ci.info$ci.lengths
    no.inf.ci <- ci.info$no.inf.ci
    ses <- ci.info$ses
    df.res <- ci.info$df.res
    gamma.min <- ci.info$gamma.min
    multi.corr <- ci.info$multi.corr
    s0 <- ci.info$s0
    alpha <- 1 - ci.info$ci.level
    pval.rank <- rank(-abs(beta.j - centers)/(ci.lengths/2))
    nsplit <- length(pval.rank) + no.inf.ci
    gamma.b <- pval.rank/nsplit
    if(multi.corr){
      if(any(is.na(s0))){stop("need s0 information to be able to create multiple testing corrected pvalues")}
      level <- (1 - alpha * gamma.b/(1 - log(gamma.min) * s0))
    } else {
      level <- (1 - alpha * gamma.b/(1 - log(gamma.min)))
    }
    a <- (1 - level)/2
    a <- 1 - a
    if(all(gamma.b <= gamma.min)){
      TRUE
    } else {
      fac <- qt(a, df.res)
      nlci <- centers - fac * ses
      nrci <- centers + fac * ses
      all(nlci[gamma.b > gamma.min] <= beta.j) && all(nrci[gamma.b > gamma.min] >= beta.j)
    }
  }
  bisectionBounds <- function(shouldcover, shouldnotcover, ci.info, verbose){
    reset.shouldnotcover <- FALSE
    if(doesItCoverGamma(beta.j = shouldnotcover, ci.info = ci.info)){
      reset.shouldnotcover <- TRUE
      if(verbose){cat("finding a new shouldnotcover bound\n")}
      while(doesItCoverGamma(beta.j = shouldnotcover, ci.info = ci.info)){
        old <- shouldnotcover
        shouldnotcover <- shouldnotcover + (shouldnotcover - shouldcover)
        shouldcover <- old
      }
      if(verbose){
        cat("new\n")
        cat("shouldnotcover", shouldnotcover, "\n")
        cat("shouldcover", shouldcover, "\n")
      }
    }
    if(!doesItCoverGamma(beta.j = shouldcover, ci.info = ci.info)){
      if(reset.shouldnotcover){stop("Problem: we first reset shouldnotcover and are now resetting shouldcover, this is not supposed to happen")}
      if(verbose){cat("finding a new shouldcover bound\n")}
      while(!doesItCoverGamma(beta.j = shouldcover, ci.info = ci.info)){
        old <- shouldcover
        shouldcover <- shouldcover + (shouldcover - shouldnotcover)
        shouldnotcover <- old
      }
      if(verbose){
        cat("new\n")
        cat("shouldnotcover", shouldnotcover, "\n")
        cat("shouldcover", shouldcover, "\n")
      }
    }
    return(list(shouldcover = shouldcover, shouldnotcover = shouldnotcover))
  }
  bisectionCoverage <- function(outer, inner, ci.info, verbose, eps.bound = 10^(-7)){
    checkBisBounds(shouldcover = inner, shouldnotcover = outer, ci.info = ci.info, verbose = verbose)
    eps <- 1
    while(eps > eps.bound){
      middle <- (outer + inner)/2
      if(doesItCoverGamma(beta.j = middle, ci.info = ci.info)){
        inner <- middle
      } else {
        outer <- middle
      }
      eps <- abs(inner - outer)
    }
    solution <- (inner + outer)/2
    if(verbose){cat("finished bisection...eps is", eps, "\n")}
    return(solution)
  }
  checkBisBounds <- function(shouldcover, shouldnotcover, ci.info, verbose){
    if(doesItCoverGamma(beta.j = shouldnotcover, ci.info = ci.info)){
      stop("shouldnotcover bound is covered! we need to decrease it even more! (PLZ implement)")
    } else if(verbose){cat("shouldnotcover bound is not covered, this is good")}
    if(doesItCoverGamma(beta.j = shouldcover, ci.info = ci.info)){
      if(verbose){cat("shouldcover is covered!, It is a good covered bound")}
    } else {
      stop("shouldcover is a bad covered bound, it is not covered!")
    }
  }
  inf.ci <- is.infinite(lci) | is.infinite(rci)
  no.inf.ci <- sum(inf.ci)
  if(verbose){cat("number of Inf ci:", no.inf.ci, "\n")}
  if((no.inf.ci == length(lci)) || (no.inf.ci >= (1 - gamma.min) * length(lci))){return(c(-Inf, Inf))}
  lci <- lci[!inf.ci]
  rci <- rci[!inf.ci]
  centers <- centers[!inf.ci]
  df.res <- df.res[!inf.ci]
  ses <- ses[!inf.ci]
  s0 <- s0[!inf.ci]
  ci.lengths <- rci - lci
  ci.info <- list(lci = lci, rci = rci, centers = centers, ci.lengths = ci.lengths,
                  no.inf.ci = no.inf.ci, ses = ses, s0 = s0, df.res = df.res,
                  gamma.min = gamma.min, multi.corr = multi.corr, ci.level = ci.level)
  inner <- findInsideGamma(low = min(centers), high = max(centers), ci.info = ci.info, verbose = verbose)
  if(is.null(inner)){
    alpha <- 1 - ci.level
    beta.j <- seq(min(centers), max(centers), length = length(centers))
    pval.rank <- rank(-abs(beta.j - centers)/(ci.lengths/2))
    nsplit <- length(pval.rank) + no.inf.ci
    gamma.b <- pval.rank/nsplit
    level <- (1 - alpha * gamma.b/(1 - log(gamma.min)))
    a <- 1 - ((1 - level)/2)
    fac <- qt(a, df.res)
    nlci <- centers - fac * ses
    nrci <- centers + fac * ses
    return(c(0, 0))
  }
  outer <- min(lci)
  new.bounds <- bisectionBounds(shouldcover = inner, shouldnotcover = outer, ci.info = ci.info, verbose = verbose)
  inner <- new.bounds$shouldcover
  outer <- new.bounds$shouldnotcover
  l.bound <- bisectionCoverage(outer = outer, inner = inner, ci.info = ci.info, verbose = verbose)
  if(verbose){cat("lower bound ci aggregated is", l.bound, "\n")}
  outer <- max(rci)
  new.bounds <- bisectionBounds(shouldcover = inner, shouldnotcover = outer, ci.info = ci.info, verbose = verbose)
  inner <- new.bounds$shouldcover
  outer <- new.bounds$shouldnotcover
  u.bound <- bisectionCoverage(inner = inner, outer = outer, ci.info = ci.info, verbose = verbose)
  if(verbose){cat("upper bound ci aggregated is", u.bound, "\n")}
  return(c(l.bound, u.bound))
}

##### stability: get empirical selection probabilities from resampled data
stability <- function(obj){
  objcoefs <- obj$freqs
  p <- length(objcoefs)
  vs <- names(objcoefs)
  allNames <- lapply(lapply(objcoefs, '[[', "predictor"), as.character)
  k <- sapply(allNames, length)
  niter <- obj$call$niter
  nlam <- obj$call$nlam - 1
  splits <- crits <- list()
  for(i in 1:2){
    s0 <- paste0("split", i)
    splits[[i]] <- lapply(1:p, function(z){
      pIts <- lapply(1:niter, function(it){
        lapply(1:k[z], function(j){
          sapply(obj[[s0]][[it]]$varMods[[z]]$allCoefs, '[[', j)
        })
      })
      pOut <- lapply(1:k[z], function(zz){
        pCoefs <- do.call(rbind, lapply(pIts, '[[', zz))
        return(colMeans(abs(sign(pCoefs))))
      })
      pOut2 <- lapply(1:k[z], function(zz){
        pCoefs2 <- do.call(rbind, lapply(pIts, '[[', zz))
        return(abs(sign(pCoefs2)))
      })
      freqOut <- do.call(cbind.data.frame, pOut)
      names(freqOut) <- allNames[[z]]
      return(list(freqs = freqOut, signs = pOut2))
    })
    crits[[i]] <- lapply(1:p, function(z){
      sapply(1:niter, function(zz){
        obj[[s0]][[zz]]$varMods[[z]]$fitobj[[3]]
      })
    })
    names(crits[[i]]) <- names(splits[[i]]) <- vs
  }
  s1 <- sapply(lapply(splits[[1]], '[[', 'freqs'), nrow)
  s2 <- sapply(lapply(splits[[2]], '[[', 'freqs'), nrow)
  if(!all(s1 == s2)){ # WARNING: THIS REMOVES SOME ELEMENTS OF THE PATHS WHEN THEY DO NOT MATCH
    swhich <- which(s1 != s2)
    smax <- sapply(swhich, function(z) which.max(c(s1[z], s2[z])))
    smin <- sapply(swhich, function(z) which.min(c(s1[z], s2[z])))
    for(i in seq_along(swhich)){
      n <- nrow(splits[[smin[i]]][[swhich[i]]]$freqs)
      nn <- ncol(splits[[smin[i]]][[swhich[i]]]$freqs)
      splits[[smax[i]]][[swhich[i]]]$freqs <- splits[[smax[i]]][[swhich[i]]]$freqs[1:n, ]
      for(j in 1:nn){
        splits[[smax[i]]][[swhich[i]]]$signs[[j]] <- splits[[smax[i]]][[swhich[i]]]$signs[[j]][, 1:n]
      }
    }
  }
  simulVars <- lapply(1:p, function(z){
    simul1 <- sapply(1:k[z], function(zz){
      colMeans(splits[[1]][[z]]$signs[[zz]] * splits[[2]][[z]]$signs[[zz]])
    })
    simul1 <- data.frame(simul1)
    colnames(simul1) <- allNames[[z]]
    return(simul1)
  })
  names(simulVars) <- vs
  for(i in 1:2){for(j in 1:p){splits[[i]][[j]] <- splits[[i]][[j]]$freqs}}
  output <- list(split1 = append(splits[[1]], list(crits = crits[[1]])),
                 split2 = append(splits[[2]], list(crits = crits[[2]])),
                 simult = simulVars)
  return(output)
}
