#' Bootstrapping network estimation for moderated networks
#'
#' Follows closely to the methods of bootstrapping found in the \code{bootnet}
#' package. An essential goal behind this function is to expand the methods in
#' \code{bootnet} to encompass moderated networks.
#'
#' Can be used to perform bootstrapped network estimation, as well as perform a
#' case-drop bootstrap. Details on these two methods can be found in the help
#' page for the \code{bootnet::bootnet} function.
#'
#' The defining feature of \code{\link{bootNet}} that differentiates it from the
#' \code{\link{resample}} function when \code{sampMethod = "bootstrap"} is that
#' the *same model is fit at every iteration* in \code{\link{bootNet}}. The only
#' time that models may differ across iterations is if a \code{threshold} is
#' specified. When \code{threshold = FALSE}, then the saturated model is fit to
#' each bootstrapped sample. Alternatively, bootstrapping can be performed with
#' respect to a specific constrained model. In this case, the constrained model
#' (variable selection model; output of \code{\link{varSelect}} or
#' \code{\link{resample}}) can be supplied to the \code{type} argument, and thus
#' this function provides a way to estimate the posterior distributions of the
#' nodes based on a constrained model.
#'
#' In addition to expanding \code{bootnet} to handle moderated networks, there
#' are also some additional features such as the capacity to perform the block
#' bootstrap for temporal networks via the \code{block} argument. The block
#' bootstrap is \strong{highly} recommended for resampling temporal networks.
#'
#' Another feature of this function is that it can be used on outputs from the
#' \code{\link{resample}} function. This can be used as a way to evaluate the
#' iterations of \code{\link{resample}} beyond just using it for variable
#' selection.
#'
#' @section Warning:
#'
#'   Importantly, if output from the \code{\link{resample}} function is used as
#'   input for the \code{\link{bootNet}} function, and the user wishes to use
#'   the model selected by the \code{\link{resample}} function as the comparison
#'   to the bootstrapped results, you must add the \code{fit0} argument to this
#'   function. Use the fitted object in the \code{\link{resample}} output as the
#'   input for the undocumented \code{fit0} argument for the
#'   \code{\link{bootNet}} function.
#'
#' @param data Dataframe or matrix.
#' @param m Numeric or character string. Indicates which variable should be
#'   treated as a moderator (if any).
#' @param nboots Number of bootstrapped samples.
#' @param lags Numeric or logical, to indicate whether or not a temporal network
#'   is being estimated. Maximum of 1 lag -- meaningful values are either 1 or
#'   \code{TRUE}.
#' @param caseDrop Logical. Determines whether to do a caseDrop bootstrap
#'   procedure or not.
#' @param rule Only applies to GGMs (including between-subjects networks) when a
#'   threshold is supplied. The \code{"AND"} rule will only preserve edges when
#'   both corresponding coefficients have p-values below the threshold, while
#'   the \code{"OR"} rule will preserve an edge so long as one of the two
#'   coefficients have a p-value below the supplied threshold.
#' @param ci Numeric, between 0 and 1. The level of the confidence intervals
#'   estimated. Defaults at .95
#' @param caseMin Numeric. The minimum proportion of the sample that should be
#'   taken when \code{caseDrop = TRUE}. Provide a value between 0 and 1. The
#'   value indicates the smallest proportion of the total sample size to test in
#'   the case-dropping procedure,
#' @param caseMax Numeric. The maximum proportion of the sample that should be
#'   taken when \code{caseDrop = TRUE}. Provide a value between 0 and 1. The
#'   value indicates the largest proportion of the total sample size to test in
#'   the case-dropping procedure,
#' @param caseN Numeric. The number of samples to draw at each sample size
#'   tested when \code{caseDrop = TRUE}.
#' @param threshold Logical or numeric. If \code{TRUE}, then a default value of
#'   .05 will be set. Indicates whether a threshold should be placed on the
#'   bootstrapped samples. A significant choice by the researcher. Only applies
#'   when a variable selection procedure is applied, or whether a
#'   \code{\link{resample}} object is used as input.
#' @param fits A list of all fitted models, if available. Not likely to be used.
#' @param type See \code{type} argument in \code{\link{fitNetwork}} function.
#'   This is where a variable selection model can be provided. This will fit the
#'   same selected model across all iterations of the bootstrapping procedure.
#' @param saveMods Logical. Determines whether or not to return all of the
#'   fitted models -- that is, all the models fit to each bootstrapped sample.
#'   Defaults to \code{TRUE}, but if \code{FALSE} then models will not be
#'   returned which can save memory.
#' @param verbose Logical. Determines whether a progress bar should be shown, as
#'   well as whether messages should be shown.
#' @param fitCoefs Logical, refers to the argument in the
#'   \code{\link{fitNetwork}} function. Most likely this should always be
#'   \code{FALSE}.
#' @param size Numeric. Size of sample to use for bootstrapping. Not
#'   recommended.
#' @param nCores If a logical or numeric value is provided, then the
#'   bootstrapping procedure will be parallelized across multiple CPUs. If
#'   numeric, this will specify the number of cores to use for the procedure. If
#'   \code{TRUE}, then the
#'   \code{\link[parallel:detectCores]{parallel::detectCores}} function of the
#'   \code{parallel} package will be run to maximize the number of cores
#'   available. Defaults to 1, which does not run any parallelization functions.
#' @param cluster Character string to indicate which type of parallelization
#'   function to use, if \code{nCores > 1}. Options are \code{"mclapply"} or
#'   \code{"SOCK"}.
#' @param block Logical or numeric. If specified, then this indicates that
#'   \code{lags != 0} or \code{lags != NULL}. If numeric, then this indicates
#'   that block bootstrapping will be used, and the value specifies the block
#'   size. If \code{TRUE} then an appropriate block size will be estimated
#'   automatically.
#' @param maxiter The maximum number of iterations for the algorithm to go
#'   through before stopping. In some circumstances, iterated versions of the
#'   model based on subsamples of the data may not be possible to fit. In these
#'   cases, \code{maxiter} specifies the number of attempts that are made with
#'   different versions of the sample before stopping the algorithm.
#' @param directedDiag logical
#' @param beepno Character string or numeric value to indicate which variable
#'   (if any) encodes the survey number within a single day. Must be used in
#'   conjunction with \code{dayno} argument.
#' @param dayno Character string or numeric value to indicate which variable (if
#'   any) encodes the survey number within a single day. Must be used in
#'   conjunction with \code{beepno} argument.
#' @param ... Additional arguments.
#'
#' @return A \code{bootNet} object
#' @export
#'
#' @seealso \code{\link{summary.bootNet}, \link{fitNetwork}, \link{varSelect},
#'   \link{resample}, \link{plotBoot}, \link{plotNet}, \link{net},
#'   \link{netInts}}
#'
#' @examples
#' \donttest{
#' boot1 <- bootNet(ggmDat, 'M')
#' summary(boot1)
#'
#' boot2 <- bootNet(gvarDat, 'M', lags = 1)
#'
#' mod1 <- varSelect(gvarDat, 'M', lags = 1)
#' boot3 <- bootNet(gvarDat, 'M', lags = 1, type = mod1, caseDrop = TRUE)
#' summary(boot3)
#' }
bootNet <- function(data, m = NULL, nboots = 10, lags = NULL, caseDrop = FALSE, rule = 'OR',
                    ci = .95, caseMin = .05, caseMax = .75, caseN = 10, threshold = FALSE,
                    fits = NULL, type = 'g', saveMods = TRUE, verbose = TRUE, fitCoefs = FALSE,
                    size = NULL, nCores = 1, cluster = 'mclapply', block = FALSE, maxiter = 10,
                    directedDiag = FALSE, beepno = NULL, dayno = NULL, ...){ # Need to add beepno and dayno
  if(identical(m, 0)){m <- NULL}
  args <- tryCatch({list(...)}, error = function(e){list()})
  call <- as.list(match.call())
  if(is.numeric(threshold)){threshold <- ifelse(threshold == 0, 'fit0', ifelse(threshold == 1, 'fits', threshold))}
  sampThresh <- switch(2 - (threshold == 'fit0'), TRUE, ifelse(threshold == 'fits', FALSE, threshold))
  threshold <- ifelse(threshold == 'fits', TRUE, ifelse(threshold == 'fit0', FALSE, threshold))
  if(any(sapply(list(data, fits), function(z) isTRUE(attr(z, 'resample'))))){
    dat <- switch(2 - isTRUE(attr(data, 'resample')), data, replace(fits, 'data', list(data = data)))
    stopifnot('data' %in% names(dat))
    caseDrop <- FALSE
    inds <- do.call(cbind.data.frame, lapply(dat$samples$iters, '[[', 'samp_inds'))
    fits <- lapply(dat$samples$iters, '[[', 'fit')
    attr(fits, 'resample') <- TRUE
    lags <- ifelse('lags' %in% names(dat$call), dat$call$lags, FALSE)
    nboots <- dat$call$niter
    m <- call$m <- dat$call$moderators
    attributes(fits)[c('inds', 'lags', 'm')] <- list(inds, lags, m)
    if('rule' %in% names(dat$call)){rule <- dat$call$rule}
    args0 <- dat$call[intersect(names(dat$call), formalArgs(fitNetwork))]
    args0$threshold <- sampThresh
    data <- args0$data <- dat$data
    if(!'fit0' %in% names(args)){
      fit0 <- do.call(fitNetwork, append(args0, list(saveMods = FALSE)))
    } else {
      fit0 <- switch(2 - isTRUE(args$fit0), modSelect(dat, data, TRUE), args$fit0)
      args$fit0 <- NULL
    }
  }
  ci <- (1 - ci)/2
  lags <- switch(2 - !is.null(lags), ifelse(all(lags != 0), 1, 0), 0)
  consec <- switch(2 - (lags & 'samp_ind' %in% names(attributes(data))), attr(data, 'samp_ind'), NULL)
  n <- n0 <- ifelse(is.null(consec), ifelse(is(data, 'list'), nrow(data[[1]]), nrow(data)), length(consec))
  if(is.null(fits)){
    if(lags & any(!sapply(c(beepno, dayno), is.null))){
      stopifnot(!is.null(beepno) & !is.null(dayno))
      data <- getConsec(data = data, beepno = beepno, dayno = dayno)
      consec <- switch(2 - (lags & 'samp_ind' %in% names(attributes(data))), attr(data, 'samp_ind'), NULL)
      n <- n0 <- ifelse(is.null(consec), nrow(data), length(consec))
    }
    sampInd <- function(size = NULL, n, lags = NULL, consec = NULL,
                        caseDrop = FALSE, block = FALSE, ind_args = NULL){
      if(!is.null(ind_args)){
        call <- as.list(match.call())[-1]
        call <- call[setdiff(names(call), 'ind_args')]
        if('size' %in% names(call)){call$size <- size}
        if(length(call) > 0){ind_args <- replace(ind_args, names(call), call)}
        out <- do.call(sampInd, ind_args)
        return(out)
      }
      if(!caseDrop){
        n0 <- n
        n <- ifelse(is.null(size), n - lags, size)
        allinds <- switch(2 - is.null(consec), seq_len(n0 - lags), consec)
        if(!identical(block, FALSE) & lags){
          N <- length(allinds); nblox <- 1
          if(isTRUE(block)){block <- floor(3.15 * N^(1/3))}
          while(N > block * nblox){nblox <- nblox + 1}
          possible <- seq_len(N - block)
          starts <- replicate(nblox, sample(possible, 1))
          out <- c(sapply(starts, function(z) allinds[z:(z + block - 1)]))
          if(length(out) > N){out <- out[1:N]}
        } else {
          out <- sample(allinds, n, replace = TRUE)
        }
      } else {
        possible <- switch(2 - is.null(consec), seq_len(n - lags), consec)
        if(lags & block){
          start <- sample(seq_len(n - size), 1)
          section <- start:(start + size - 1)
          out <- possible[section]
        } else {
          out <- sample(possible, size, replace = FALSE)
        }
      }
      return(out)
    }
    indArgs <- list(size = size, n = n, lags = lags, consec = consec, caseDrop = caseDrop, block = block)
    if(!caseDrop){
      inds <- replicate(nboots, list(sampInd(ind_args = indArgs)))
    } else {
      mm <- as.numeric(!is.null(m))
      p <- ncol(data) - mm
      fixMax <- ifelse(lags, round((1 - caseMax) * n) < ((p + 1) * (mm + 2)),
                       round((1 - caseMax) * n) < sampleSize(p = p, m = mm, print = FALSE))
      if(fixMax){
        caseMax <- 1 - ifelse(lags, ((p + 1) * (mm + 2)), sampleSize(p = p, m = mm, print = FALSE))/n
        message(paste0('caseMax too large -- setting to max value: ', round(caseMax, 2)))
      }
      subCases <- round((1 - seq(caseMin, caseMax, length = caseN)) * n)
      subNs <- sample(subCases, nboots, replace = TRUE)
      inds <- lapply(subNs, sampInd, ind_args = indArgs)
    }
    if(verbose & identical(nCores, 1)){pb <- txtProgressBar(min = 0, max = nboots + 1, style = 3)}
    args0 <- list(data = data, moderators = m, type = type, lags = lags, rule = rule,
                  threshold = sampThresh, verbose = FALSE, saveMods = FALSE,
                  fitCoefs = fitCoefs, maxiter = maxiter)
    args <- args[setdiff(names(args), names(args0))]
    args0 <- append(args0, args[intersect(names(args), formalArgs(fitNetwork))])
    fit0 <- do.call(fitNetwork, args0)
    if('fit0' %in% names(args)){fit0 <- args$fit0}
    if(verbose & identical(nCores, 1)){setTxtProgressBar(pb, 1)}
    args0$threshold <- threshold
    if(nCores > 1 | isTRUE(nCores)){
      if(isTRUE(nCores)){nCores <- parallel::detectCores()}
      if(grepl('Windows', sessionInfo()$running)){cluster <- 'SOCK'}
      if(tolower(cluster) != 'mclapply'){
        cluster <- match.arg(toupper(cluster), c('SOCK', 'FORK'))
        cl <- parallel::makeCluster(nCores, type = cluster)
        if(cluster == 'SOCK'){
          obj1 <- switch(2 - !as.logical(lags), c('nodewise', 'modNet', 'modLL'),
                         c('lagMat', 'SURfit', 'SURnet', 'SURll', 'surEqs', 'getCoefs', 'systemfit'))
          obj1 <- c(obj1, 'lags', 'inds', 'data', 'args0', 'm', 'caseDrop', 'indArgs', 'sampInd', 'maxiter')
          parallel::clusterExport(cl, c('fitNetwork', 'Matrix', 'net', 'netInts', obj1), envir = environment())
        }
      } else {
        cl <- nCores
      }
      if(verbose){
        pbapply::pboptions(type = 'timer', char = '-')
        fits <- suppressWarnings(pbapply::pblapply(seq_len(nboots), function(z){
          newargs <- args0
          inds0 <- inds[[z]]
          fit <- tt <- 0
          while(identical(fit, 0)){
            newargs$data <- switch(lags + 1, data[inds0, ], structure(data, samp_ind = inds0))
            fit <- tryCatch({do.call(fitNetwork, newargs)}, error = function(e){0})
            if(identical(fit, 0)){
              tt <- tt + 1
              if(tt == maxiter){break}
              inds0 <- sampInd(size = length(inds0), ind_args = indArgs)
            } else {
              attr(fit, 'inds') <- inds0
            }
          }
          if(identical(fit, 0)){fit <- structure(list(), class = 'try-error')}
          return(fit)
        }, cl = cl))
      } else if(tolower(cluster) == 'mclapply'){
        fits <- suppressWarnings(parallel::mclapply(seq_len(nboots), function(z){
          newargs <- args0
          inds0 <- inds[[z]]
          fit <- tt <- 0
          while(identical(fit, 0)){
            newargs$data <- switch(lags + 1, data[inds0, ], structure(data, samp_ind = inds0))
            fit <- tryCatch({do.call(fitNetwork, newargs)}, error = function(e){0})
            if(identical(fit, 0)){
              tt <- tt + 1
              if(tt == maxiter){break}
              inds0 <- sampInd(size = length(inds0), ind_args = indArgs)
            } else {
              attr(fit, 'inds') <- inds0
            }
          }
          if(identical(fit, 0)){fit <- structure(list(), class = 'try-error')}
          return(fit)
        }, mc.cores = nCores))
      } else {
        fits <- suppressWarnings(parallel::parLapply(cl, seq_len(nboots), function(z){
          newargs <- args0
          inds0 <- inds[[z]]
          fit <- tt <- 0
          while(identical(fit, 0)){
            newargs$data <- switch(lags + 1, data[inds0, ], structure(data, samp_ind = inds0))
            fit <- tryCatch({do.call(fitNetwork, newargs)}, error = function(e){0})
            if(identical(fit, 0)){
              tt <- tt + 1
              if(tt == maxiter){break}
              inds0 <- sampInd(size = length(inds0), ind_args = indArgs)
            } else {
              attr(fit, 'inds') <- inds0
            }
          }
          if(identical(fit, 0)){fit <- structure(list(), class = 'try-error')}
          return(fit)
        }))
      }
      if(tolower(cluster) != 'mclapply'){parallel::stopCluster(cl)}
      rm(cl)
    } else {
      fits <- suppressWarnings(lapply(seq_len(nboots), function(z){
        newargs <- args0
        inds0 <- inds[[z]]
        fit <- tt <- 0
        while(identical(fit, 0)){
          newargs$data <- switch(lags + 1, data[inds0, ], structure(data, samp_ind = inds0))
          fit <- tryCatch({do.call(fitNetwork, newargs)}, error = function(e){0})
          if(identical(fit, 0)){
            tt <- tt + 1
            if(tt == maxiter){break}
            inds0 <- sampInd(size = length(inds0), ind_args = indArgs)
          } else {
            attr(fit, 'inds') <- inds0
          }
        }
        if(identical(fit, 0)){fit <- structure(list(), class = 'try-error')}
        if(verbose){setTxtProgressBar(pb, z + 1)}
        return(fit)
      }))
    }
    if(any(sapply(fits, class) == 'try-error')){
      err <- which(sapply(fits, class) == 'try-error')
      inds <- inds[-err]
      nboots <- nboots - length(err)
      fits <- fits[-err]
      if(caseDrop){subNs <- subNs[-err]}
      if(verbose){message(paste0('\n', length(err), ' iterations failed'))}
    }
    inds <- lapply(fits, attr, 'inds')
    fits <- structure(fits, class = c('list', 'bootFits'), inds = inds, lags = lags, m = m)
    net2 <- net3 <- FALSE
  } else {
    if(is(fits, 'bootNet')){
      if(!'fit0' %in% names(args)){fit0 <- fits$fit0}
      if('bootFits' %in% names(fits)){fits <- fits$bootFits} else {stop('No fits')}
    }
    if('lags' %in% names(attributes(fits))){lags <- attr(fits, 'lags')}
    if('inds' %in% names(attributes(fits))){inds <- attr(fits, 'inds')}
    if('m' %in% names(attributes(fits))){m <- attr(fits, 'm')}
    net2 <- ifelse('net2' %in% names(args), args$net2, FALSE)
    net3 <- ifelse('net3' %in% names(args), args$net3, FALSE)
    nboots <- length(fits)
    if(length(args) > 0){
      if('fit0' %in% names(args)){fit0 <- args$fit0}
      if('inds' %in% names(args)){inds <- args$inds}
      if('subNs' %in% names(args)){subNs <- args$subNs}
    }
  }
  if(directedDiag){
    strFUN <- switch(2 - net3, function(x){rowSums(x)}, function(x){colSums(x)})
  } else {
    strFUN <- switch(2 - net3, function(x){diag(x) <- 0; rowSums(x)}, function(x){diag(x) <- 0; colSums(x)})
  }
  nodes <- switch(lags + 1, names(fit0$mods), gsub('[.]y$', '', names(fit0$SURnet$mods)))
  v <- length(nodes)
  node01 <- node1 <- rep(nodes[-v], (v - 1):1)
  node02 <- node2 <- matrix(nodes, v, v)[lower.tri(matrix(nodes, v, v))]
  if(lags & net2){
    node1 <- rep(nodes, each = v)
    node2 <- rep(nodes, v)
  } else if(lags){
    m <- NULL
  }
  edges <- paste0(node1, ifelse(lags & net2, '-->', '--'), node2)
  e <- length(edges)
  edge1 <- rep(edges[-e], (e - 1):1)
  edge2 <- matrix(edges, e, e)[lower.tri(matrix(edges, e, e))]
  adj0 <- switch(factorial(lags + net2),
                 net(fit0, 'PCC', sampThresh, rule)[lower.tri(net(fit0, 'PCC', sampThresh, rule))],
                 c(net(fit0, threshold = sampThresh, rule = rule)))
  which.net <- ifelse(lags & net2, 'beta', 'PCC')
  str0 <- strFUN(abs(net(fit0, which.net, sampThresh, rule)))
  ei0 <- strFUN(net(fit0, which.net, sampThresh, rule))
  adj1 <- do.call(cbind, lapply(fits, function(z){
    switch(factorial(lags + net2),
           net(z, 'PCC', threshold, rule)[lower.tri(net(z, 'PCC', threshold, rule))],
           c(net(z, threshold = threshold, rule = rule)))
  }))
  colnames(adj1) <- names(inds) <- paste0('boot', 1:nboots)
  if(!caseDrop){
    adj2 <- as.data.frame(t(apply(adj1, 1, function(z) quantile(z, c(ci, 1 - ci)))))
    colnames(adj2) <- c('boot_lower', 'boot_upper')
    adjmeans <- data.frame(boot_mean = rowMeans(adj1), boot_sd = apply(adj1, 1, sd))
    adj2 <- data.frame(edge = edges, node1 = node1, node2 = node2, adjmeans, adj2,
                       sample = adj0, sample_lower = adj0 + (adjmeans$boot_sd * qnorm(ci)),
                       sample_upper = adj0 + (adjmeans$boot_sd * qnorm(1 - ci)))
  } else {
    uniqueNs <- sort(unique(subNs))
    uN <- length(uniqueNs)
    if(any(adj0 == 0)){
      groupEdges <- lapply(split(data.frame(t(adj1))[, -which(adj0 == 0)], subNs),
                           function(z){apply(z, 1, function(zz) cor(zz, adj0[adj0 != 0]))})
    } else {
      groupEdges <- lapply(split(data.frame(t(adj1)), subNs),
                           function(z){apply(z, 1, function(zz) cor(zz, adj0))})
    }
    adjmeans <- lapply(groupEdges, function(z) c(mean = mean(z), sd = sd(z)))
    qs <- c(.01, .025, .05, .25, .5, .75, .95, .975, .99)
    adjqs <- data.frame(t(sapply(groupEdges, quantile, probs = qs)))
    colnames(adjqs) <- paste0('q', c(1, 2.5, 5, 25, 50, 75, 95, 97.5, 99))
    adjmeans <- data.frame(do.call(rbind, adjmeans))
    adj2 <- data.frame(N = rep(n, uN), drop = round(1 - uniqueNs/n, 2),
                       subN = uniqueNs, n = sapply(groupEdges, length),
                       adjmeans, adjqs)
    adj2 <- adj2[order(adj2$subN, decreasing = FALSE), ]
    rownames(adj2) <- 1:nrow(adj2)
  }
  rownames(adj1) <- edges
  if(!is.null(m)){
    mname <- ifelse(length(fit0$call$moderators) == 1, fit0$call$moderators, 'MOD')
    inodes <- paste0(nodes, ':', mname); iedges <- paste0(edges, '|', mname)
    inode1 <- paste0(node1, ':', mname); inode2 <- paste0(node2, ':', mname)
    iedge1 <- paste0(edge1, '|', mname); iedge2 <- paste0(edge2, '|', mname)
    if(isTRUE(attr(fit0, 'ggm'))){
      ints0 <- netInts(fit0, threshold = threshold, rule = rule, avg = TRUE)
      istr0 <- colSums(abs(ints0))
      iei0 <- colSums(ints0)
      ints0 <- ints0[lower.tri(ints0)]
      ints1 <- do.call(cbind, lapply(1:nboots, function(z){
        i1 <- netInts(fits[[z]], threshold = threshold, rule = rule, avg = TRUE)
        if(length(i1) == 0){i1 <- diag(0, v)}
        return(i1[lower.tri(i1)])
      }))
      #ints0 <- fit0$interactions$intMats$avgInts[lower.tri(fit0$interactions$intMats$avgInts)]
      #istr0 <- colSums(abs(fit0$interactions$intMats$avgInts))
      #iei0 <- colSums(fit0$interactions$intMats$avgInts)
      #ints1 <- do.call(cbind, lapply(1:nboots, function(z){
      #  fits[[z]]$interactions$intMats$avgInts[lower.tri(fits[[z]]$interactions$intMats$avgInts)]
      #}))
    } else {
      ints0 <- c(netInts(fit0, threshold = sampThresh))
      istr0 <- strFUN(abs(netInts(fit0, threshold = sampThresh)))
      iei0 <- strFUN(netInts(fit0, threshold = sampThresh))
      ints1 <- do.call(cbind, lapply(1:nboots, function(z){
        z1 <- netInts(fits[[z]], threshold = threshold)
        if(length(z1) == 0){z1 <- diag(0, v)}
        return(c(z1))
      }))
    }
    colnames(ints1) <- paste0('boot', 1:nboots)
    if(!caseDrop){
      ints2 <- as.data.frame(t(apply(ints1, 1, function(z) quantile(z, c(ci, 1 - ci)))))
      colnames(ints2) <- c('boot_lower', 'boot_upper')
      intsmeans <- data.frame(boot_mean = rowMeans(ints1), boot_sd = apply(ints1, 1, sd))
      ints2 <- data.frame(edge = iedges, node1 = inode1, node2 = inode2,
                          intsmeans, ints2, sample = ints0,
                          sample_lower = ints0 + (intsmeans$boot_sd * qnorm(ci)),
                          sample_upper = ints0 + (intsmeans$boot_sd * qnorm(1 - ci)))
    } else {
      if(any(ints0 == 0)){
        groupInts <- lapply(split(data.frame(t(ints1))[, -which(ints0 == 0)], subNs),
                            function(z){apply(z, 1, function(zz) cor(zz, ints0[ints0 != 0]))})
      } else {
        groupInts <- lapply(split(data.frame(t(ints1)), subNs),
                            function(z){apply(z, 1, function(zz) cor(zz, ints0))})
      }
      intsmeans <- lapply(groupInts, function(z) c(mean = mean(z), sd = sd(z)))
      intsqs <- data.frame(t(sapply(groupInts, quantile, probs = qs)))
      colnames(intsqs) <- paste0('q', c(1, 2.5, 5, 25, 50, 75, 95, 97.5, 99))
      intsmeans <- data.frame(do.call(rbind, intsmeans))
      ints2 <- data.frame(N = rep(n, uN), drop = round(1 - uniqueNs/n, 2),
                          subN = uniqueNs, n = sapply(groupInts, length),
                          intsmeans, intsqs)
      ints2 <- ints2[order(ints2$subN, decreasing = FALSE), ]
      rownames(ints2) <- 1:nrow(ints2)
    }
    rownames(ints1) <- iedges
  }
  strengths <- intStr <- EIs <- intEIs <- list()
  if(lags & net2){
    fnets <- lapply(fits, net, which.net, threshold, rule)
    strengths <- data.frame(do.call(rbind, lapply(fnets, function(z) strFUN(abs(z)))))
    EIs <- data.frame(do.call(rbind, lapply(fnets, function(z) strFUN(z))))
    if(!is.null(m)){
      fnets2 <- lapply(fits, netInts, threshold = threshold)
      intStr <- data.frame(do.call(rbind, lapply(fnets2, function(z) strFUN(abs(z)))))
      intEIs <- data.frame(do.call(rbind, lapply(fnets2, function(z) strFUN(z))))
    }
    node1 <- node01
    node2 <- node02
  } else {
    for(i in 1:length(nodes)){
      strengths[[i]] <- adj1[which(node1 == nodes[i] | node2 == nodes[i]), ]
      EIs[[i]] <- colSums(strengths[[i]])
      strengths[[i]] <- colSums(abs(strengths[[i]]))
      if(!is.null(m)){
        intStr[[i]] <- ints1[which(node1 == nodes[i] | node2 == nodes[i]), ]
        intEIs[[i]] <- colSums(intStr[[i]])
        intStr[[i]] <- colSums(abs(intStr[[i]]))
      }
    }
    strengths <- do.call(cbind.data.frame, strengths)
    EIs <- do.call(cbind.data.frame, EIs)
  }
  colnames(strengths) <- colnames(EIs) <- nodes
  if(!is.null(m)){
    if(!(lags & net2)){
      intStr <- do.call(cbind.data.frame, intStr)
      intEIs <- do.call(cbind.data.frame, intEIs)
    }
    colnames(intStr) <- colnames(intEIs) <- inodes
  }
  if(!caseDrop){
    estDiffs <- function(ests, ci = .025){
      comps <- t(combn(ncol(ests), 2))
      diffOut <- as.data.frame(t(sapply(1:nrow(comps), function(z){
        ints <- quantile(ests[, comps[z, 1]] - ests[, comps[z, 2]], c(ci, 1 - ci))
        return(c(ints, ifelse(ints[1] <= 0 & ints[2] >= 0, FALSE, TRUE)))
      })))
      colnames(diffOut) <- c('lower', 'upper', 'sig')
      diffOut$sig <- as.logical(diffOut$sig)
      return(diffOut)
    }
    sdiffs <- data.frame(node1 = node1, node2 = node2, estDiffs(strengths, ci))
    eidiffs <- data.frame(node1 = node1, node2 = node2, estDiffs(EIs, ci))
    smeans <- data.frame(node = nodes, boot_mean = colMeans(strengths), boot_sd = apply(strengths, 2, sd))
    squants <- as.data.frame(t(apply(strengths, 2, function(z) quantile(z, c(ci, 1 - ci)))))
    eimeans <- data.frame(node = nodes, boot_mean = colMeans(EIs), boot_sd = apply(EIs, 2, sd))
    eiquants <- as.data.frame(t(apply(EIs, 2, function(z) quantile(z, c(ci, 1 - ci)))))
    colnames(squants) <- colnames(eiquants) <- c('boot_lower', 'boot_upper')
    smeans <- data.frame(smeans, squants, sample = str0, sample_lower = str0 + (smeans$boot_sd * qnorm(ci)),
                         sample_upper = str0 + (smeans$boot_sd * qnorm(1 - ci)))
    eimeans <- data.frame(eimeans, eiquants, sample = ei0, sample_lower = ei0 + (eimeans$boot_sd * qnorm(ci)),
                          sample_upper = ei0 + (eimeans$boot_sd * qnorm(1 - ci)))
    rownames(smeans) <- rownames(eimeans) <- 1:v
    edge_diffs <- data.frame(edge1 = edge1, edge2 = edge2, estDiffs(t(adj1), ci))
  } else {
    pairStr <- pairEI <- multiStr <- multiEI <- vector('list', length = uN)
    psmeans <- pemeans <- msmeans <- memeans <- vector('list', length = uN)
    for(i in 1:uN){
      colN <- which(subNs == uniqueNs[i])
      pairStr[[i]] <- apply(t(strengths)[, colN, drop = FALSE], 2, function(z) cor(z, str0))
      pairEI[[i]] <- apply(t(EIs)[, colN, drop = FALSE], 2, function(z) cor(z, ei0))
      psmeans[[i]] <- c(mean = mean(pairStr[[i]]), sd = sd(pairStr[[i]]))
      pemeans[[i]] <- c(mean = mean(pairEI[[i]]), sd = sd(pairEI[[i]]))
      if(i == uN){
        psmeans <- data.frame(do.call(rbind, psmeans))
        pemeans <- data.frame(do.call(rbind, pemeans))
        psqs <- data.frame(t(sapply(pairStr, quantile, probs = qs)))
        peqs <- data.frame(t(sapply(pairEI, quantile, probs = qs)))
        colnames(psqs) <- colnames(peqs) <- paste0('q', c(1, 2.5, 5, 25, 50, 75, 95, 97.5, 99))
        pairStrengths <- data.frame(N = rep(n, uN), drop = round(1 - uniqueNs/n, 2),
                                    subN = uniqueNs, n = sapply(pairStr, length),
                                    psmeans, psqs)
        pairExInf <- data.frame(N = rep(n, uN), drop = round(1 - uniqueNs/n, 2),
                                subN = uniqueNs, n = sapply(pairEI, length),
                                pemeans, peqs)
        pairStrengths <- pairStrengths[order(pairStrengths$subN, decreasing = FALSE), ]
        pairExInf <- pairExInf[order(pairExInf$subN, decreasing = FALSE), ]
        rownames(pairExInf) <- rownames(pairStrengths) <- 1:nrow(pairStrengths)
      }
      if(!is.null(m)){
        multiStr[[i]] <- apply(t(intStr)[, colN, drop = FALSE], 2, function(z) cor(z, istr0))
        multiEI[[i]] <- apply(t(intEIs)[, colN, drop = FALSE], 2, function(z) cor(z, iei0))
        msmeans[[i]] <- c(mean = mean(multiStr[[i]]), sd = sd(multiStr[[i]]))
        memeans[[i]] <- c(mean = mean(multiEI[[i]]), sd = sd(multiEI[[i]]))
        if(i == uN){
          msmeans <- data.frame(do.call(rbind, msmeans))
          memeans <- data.frame(do.call(rbind, memeans))
          msqs <- data.frame(t(sapply(multiStr, quantile, probs = qs)))
          meqs <- data.frame(t(sapply(multiEI, quantile, probs = qs)))
          colnames(msqs) <- colnames(meqs) <- paste0('q', c(1, 2.5, 5, 25, 50, 75, 95, 97.5, 99))
          intStrengths <- data.frame(N = rep(n, uN), drop = round(1 - uniqueNs/n, 2),
                                     subN = uniqueNs, n = sapply(multiStr, length),
                                     msmeans, msqs)
          intExInf <- data.frame(N = rep(n, uN), drop = round(1 - uniqueNs/n, 2),
                                 subN = uniqueNs, n = sapply(multiEI, length),
                                 memeans, meqs)
          intStrengths <- intStrengths[order(intStrengths$subN, decreasing = FALSE), ]
          intExInf <- intExInf[order(intExInf$subN, decreasing = FALSE), ]
          rownames(intExInf) <- rownames(intStrengths) <- 1:nrow(intStrengths)
        }
      }
    }
  }
  if(!is.null(m) & !caseDrop){
    inode1 <- paste0(node1, ':', mname); inode2 <- paste0(node2, ':', mname)
    iedge1 <- paste0(edge1, '|', mname); iedge2 <- paste0(edge2, '|', mname)
    int_squants <- as.data.frame(t(apply(intStr, 2, function(z) quantile(z, c(ci, 1 - ci)))))
    int_eiquants <- as.data.frame(t(apply(intEIs, 2, function(z) quantile(z, c(ci, 1 - ci)))))
    colnames(int_squants) <- colnames(int_eiquants) <- c('boot_lower', 'boot_upper')
    int_edge_diffs <- data.frame(edge1 = iedge1, edge2 = iedge2, estDiffs(t(ints1), ci))
    int_sdiffs <- data.frame(node1 = inode1, node2 = inode2, estDiffs(intStr, ci))
    int_eidiffs <- data.frame(node1 = inode1, node2 = inode2, estDiffs(intEIs, ci))
    int_smeans <- data.frame(node = inodes, boot_mean = colMeans(intStr), boot_sd = apply(intStr, 2, sd))
    int_smeans <- data.frame(int_smeans, int_squants, sample = istr0,
                             sample_lower = istr0 + (int_smeans$boot_sd * qnorm(ci)),
                             sample_upper = istr0 + (int_smeans$boot_sd * qnorm(1 - ci)))
    int_eimeans <- data.frame(node = inodes, boot_mean = colMeans(intEIs), boot_sd = apply(intEIs, 2, sd))
    int_eimeans <- data.frame(int_eimeans, int_eiquants, sample = iei0,
                              sample_lower = iei0 + (int_eimeans$boot_sd * qnorm(ci)),
                              sample_upper = iei0 + (int_eimeans$boot_sd * qnorm(1 - ci)))
    rownames(int_smeans) <- rownames(int_eimeans) <- 1:v
    out <- list(pairwise = list(
      edges = adj2, strength = smeans, EI = eimeans, diffs = list(
        edge_diffs = edge_diffs, str_diffs = sdiffs, EI_diffs = eidiffs)),
      interactions = list(edges = ints2, strength = int_smeans, EI = int_eimeans, diffs = list(
        int_edge_diffs = int_edge_diffs, int_str_diffs = int_sdiffs, int_EI_diffs = int_eidiffs)),
      boots = list(boot_edges = adj1, boot_strengths = strengths, boot_EIs = EIs, boot_int_edges = ints1,
                   boot_int_strengths = intStr, boot_int_EI = intEIs, boot_inds = inds))
  } else if(is.null(m) & !caseDrop){
    out <- list(edges = adj2, strength = smeans, EI = eimeans,
                diffs = list(edge_diffs = edge_diffs, str_diffs = sdiffs, EI_diffs = eidiffs),
                boots = list(boot_edges = adj1, boot_strengths = strengths,
                             boot_EIs = EIs, boot_inds = inds))
  } else {
    pairwise <- list(edges = adj2, strength = pairStrengths, EI = pairExInf)
    boots <- list(boot_edges = adj1, boot_strengths = strengths, boot_EIs = EIs)
    out <- append(pairwise, list(boots = boots))
    if(!is.null(m)){
      interactions <- list(edges = ints2, strength = intStrengths, EI = intExInf)
      boots2 <- list(boot_int_edges = ints1, boot_int_strengths = intStr, boot_int_EI = intEIs)
      out <- list(pairwise = pairwise, interactions = interactions, boots = append(boots, boots2))
    }
    out$boots$boot_inds <- inds
    out$subN <- subNs
  }
  if(lags & !net2){
    if(isTRUE(attr(fits, 'resample'))){
      attr(fits, 'resample') <- NULL
      if(!'nboots' %in% names(call)){call$nboots <- nboots}
      caseDrop <- FALSE
    }
    args1 <- list(data = data, fit0 = fit0, inds = inds, net2 = TRUE,
                  fits = fits, lags = lags, saveMods = FALSE)
    args2 <- append(args1, args[setdiff(names(args), names(args1))])
    args2 <- append(args2, call[-1][setdiff(names(call[-1]), names(args2))])
    if(caseDrop){args2$subNs <- subNs} else {args2$caseDrop <- FALSE}
    out2 <- tryCatch({do.call(match.fun(call[[1]]), args2)}, error = function(e){list()})
    out3 <- tryCatch({do.call(match.fun(call[[1]]), append(args2, list(net3 = TRUE)))}, error = function(e){list()})
    attributes(out)[c('ci', 'caseDrop', 'class')] <- attributes(out2)[c('ci', 'caseDrop', 'class')]
    attr(out, 'net') <- 'contemporaneous'; attr(out2, 'net') <- 'temporal'
    makeTemporal <- function(out2, out3, caseDrop, subNs = NULL){
      makeTemp <- function(out2, out3, caseDrop, subNs = NULL){
        temp1 <- append(list(edges = out2$edges), list(
          strength = list(outStrength = out2$strength, inStrength = out3$strength),
          EI = list(outEI = out2$EI, inEI = out3$EI)
        ))
        boot1 <- append(list(boot_edges = out2$boots$boot_edges), list(
          boot_strengths = list(boot_outStrengths = out2$boots$boot_strengths, boot_inStrengths = out3$boots$boot_strengths),
          boot_EIs = list(boot_outEIs = out2$boots$boot_EIs, boot_inEIs = out3$boots$boot_EIs),
          boot_inds = out2$boots$boot_inds
        ))
        if(!caseDrop){
          diffs1 <- append(list(edge_diffs = out2$diffs$edge_diffs), list(
            str_diffs = list(outStr_diffs = out2$diffs$str_diffs, inStr_diffs = out3$diffs$str_diffs),
            EI_diffs = list(outEI_diffs = out2$diffs$EI_diffs, inEI_diffs = out3$diffs$EI_diffs)
          ))
          out2 <- append(temp1, list(diffs = diffs1, boots = boot1))
        } else {
          out2 <- append(temp1, list(boots = boot1, subN = subNs))
        }
        return(out2)
      }
      if(!'interactions' %in% names(out2)){
        out <- makeTemp(out2, out3, caseDrop, subNs)
      } else {
        bootNames <- names(out2$boots)
        pairs2 <- append(out2$pairwise, list(boots = out2$boots[!grepl('int', names(out2$boots))]))
        pairs3 <- append(out3$pairwise, list(boots = out3$boots[!grepl('int', names(out3$boots))]))
        ints2 <- append(out2$interactions, list(boots = out2$boots[grepl('int|inds', names(out2$boots))]))
        ints3 <- append(out3$interactions, list(boots = out3$boots[grepl('int|inds', names(out3$boots))]))
        names(ints2$boots) <- names(ints3$boots) <- names(pairs2$boots)
        names(ints2$diffs) <- names(ints3$diffs) <- names(pairs2$diffs)
        pairwise <- makeTemp(pairs2, pairs3, caseDrop, subNs)
        interactions <- makeTemp(ints2, ints3, caseDrop, subNs)
        if(caseDrop){pairwise$subN <- interactions$subN <- NULL}
        out <- list(pairwise = pairwise, interactions = interactions)
      }
      return(out)
    }
    atts1 <- attributes(out2)[-1]
    out2 <- makeTemporal(out2, out3, caseDrop, subNs)
    attributes(out2)[names(atts1)] <- atts1
    out <- list(temporal = out2, contemporaneous = out, fit0 = fit0)
    if(caseDrop){out$subN <- subNs}
  }
  names(fits) <- paste0('boot', 1:nboots)
  if(!lags){out$fit0 <- fit0}
  if(saveMods){out$bootFits <- fits}
  attr(out, 'n') <- n
  attr(out, 'ci') <- 1 - ci * 2
  attr(out, 'caseDrop') <- caseDrop
  class(out) <- c('list', 'bootNet')
  return(out)
}



#' Descriptive statistics for \code{\link{bootNet}} objects
#'
#' Currently only works for GGMs, including the between-subjects network
#' returned in the \code{\link{mlGVAR}} output.
#'
#' Outputs correlation-stability (CS) coefficients for the case-drop bootstrap.
#'
#' @param object \code{\link{bootNet}} output
#' @param centrality Logical. Determines whether or not strength centrality and
#'   expected influence should be computed for output.
#' @param cor Numeric value to indicate the correlation stability value to be
#'   computed.
#' @param ci Numeric. Confidence interval level for CS coefficient.
#' @param first Logical. Whether or not to count the first instance that a CS
#'   statistic dips below the requisite threshold. Often times the value of this
#'   will not affect the results. When it does, if \code{first = TRUE} then the
#'   calculation will be more conservative.
#' @param verbose Logical. Whether to write out the full statement of the CS
#'   coefficient in output. Set to \code{FALSE} if you want the details about
#'   the CS coefficient saved as attributes on the output.
#' @param ... Additional arguments.
#'
#' @return A table of descriptives for \code{\link{bootNet}} objects, or
#'   correlation-stability coefficients for the case-drop bootstrap.
#' @export
#' @name bootNetDescriptives
#'
#' @seealso \code{\link{bootNet}}
#'
#' @examples
#' \donttest{
#' boot1 <- bootNet(ggmDat, 'M')
#' summary(boot1)
#'
#' boot2 <- bootNet(gvarDat, 'M', lags = 1)
#'
#' mod1 <- varSelect(gvarDat, 'M', lags = 1)
#' boot3 <- bootNet(gvarDat, 'M', lags = 1, type = mod1, caseDrop = TRUE)
#' summary(boot3)
#' }
summary.bootNet <- function(object, centrality = TRUE, ...){
  if(isTRUE(attr(object, 'caseDrop'))){return(cscoef(object, ...))}
  inds1 <- c('edges', 'strength', 'EI')
  inds2 <- paste0('boot_', inds1)
  inds2[2:3] <- paste0(inds2[2:3], 's')
  inds <- c('pairwise', 'interactions')
  out <- vector('list', 2)
  for(i in 1:ifelse('interactions' %in% names(object), 2, 1)){
    if('interactions' %in% names(object)){
      samps <- lapply(inds1, function(z) object[[inds[i]]][[z]]$sample)
      if(i == 2){inds2 <- setdiff(names(object$boots), c(inds2, 'boot_inds'))}
    } else {
      samps <- lapply(inds1, function(z) object[[z]]$sample)
    }
    cors <- lapply(1:3, function(i){
      boots <- object$boots[[inds2[i]]]
      if(!grepl('edge', inds1[i])){boots <- t(boots)}
      sapply(seq_len(ncol(boots)), function(z) cor0(samps[[i]], boots[, z]))
    })
    mae <- lapply(1:3, function(i){
      boots <- object$boots[[inds2[i]]]
      if(!grepl('edge', inds1[i])){boots <- t(boots)}
      sapply(seq_len(ncol(boots)), function(z) mean(abs(samps[[i]] - boots[, z])))
    })
    boots <- object$boots$boot_edges
    sens <- sapply(seq_len(ncol(boots)), function(z){
      sum(boots[, z] != 0 & samps[[1]] != 0)/sum(samps[[1]] != 0)
    })
    spec <- sapply(seq_len(ncol(boots)), function(z){
      sum(boots[, z] == 0 & samps[[1]] == 0)/sum(samps[[1]] == 0)
    })
    prec <- sapply(seq_len(ncol(boots)), function(z){
      sum(boots[, z] != 0 & samps[[1]] != 0)/sum(boots[, z] != 0)
    })
    accu <- sapply(seq_len(ncol(boots)), function(z){
      tp <- sum(boots[, z] != 0 & samps[[1]] != 0)
      tn <- sum(boots[, z] == 0 & samps[[1]] == 0)
      (tp + tn)/length(samps[[1]])
    })
    if(any(is.na(sens))){sens[is.na(sens)] <- 0}
    if(any(is.na(spec))){spec[is.na(spec)] <- 0}
    if(any(is.na(prec))){prec[is.na(prec)] <- 0}
    if(any(is.na(accu))){accu[is.na(accu)] <- 0}
    out[[i]] <- cbind.data.frame(cor = cors[[1]], mae = mae[[1]], sensitivity = sens,
                                 specificity = spec, precision = prec, accuracy = accu,
                                 N = attr(object, 'n'), p = length(samps[[2]]),
                                 sparsity = sum(samps[[1]] == 0)/length(samps[[1]]),
                                 index = 1, iter = 1:ncol(boots), type = capitalize(inds[i]),
                                 network = 'Between')
    if(centrality){out[[i]] <- cbind.data.frame(out[[i]], strength = cors[[2]], EI = cors[[3]])}
    class(out[[i]]) <- c('data.frame', 'mnetPowerSim')
  }
  out <- do.call(rbind, out)
  return(out)
}

#' @rdname bootNetDescriptives
#' @export
cscoef <- function(object, cor = .7, ci = .95, first = TRUE, verbose = TRUE){
  stopifnot(isTRUE(attr(object, 'caseDrop')))
  ci0 <- ci * 100
  ci <- paste0("q", c((1 - ci)/2, 1 - (1 - ci)/2) * 100)[1]
  lags <- any(c('temporal', 'contemporaneous') %in% names(object))
  inds <- c('edges', switch(2 - !lags, c('strength', 'EI'), c('outStrength', 'inStrength', 'outEI', 'inEI')))
  if(lags){
    obj0 <- object
    object <- object$temporal
  }
  if('interactions' %in% names(object)){
    inds2 <- c('pairwise', 'interactions')
    out <- do.call(rbind, setNames(lapply(inds2, function(z){
      obj2 <- object[[z]]
      setNames(sapply(inds, function(i){
        if(grepl('^out|^in', i)){
          ii <- tolower(gsub('^out|^in', '', i))
          if(ii == 'ei'){ii <- 'EI'}
          dat <- obj2[[ii]][[i]]
          dat <- dat[order(dat$subN), ]
        } else {
          dat <- obj2[[i]][order(obj2[[i]]$subN), ]
        }
        if(!first){
          css <- dat[dat[, ci] > cor, 'drop']
          css <- ifelse(length(css) == 0, 0, max(css))
        } else {
          css <- which(!dat[, ci] > cor)
          css <- dat[ifelse(length(css) == 0, 1, max(css) + 1), 'drop']
        }
        return(css)
      }), inds)
    }), inds2))
  } else {
    out <- setNames(sapply(inds, function(i){
      if(grepl('^out|^in', i)){
        ii <- tolower(gsub('^out|^in', '', i))
        if(ii == 'ei'){ii <- 'EI'}
        dat <- object[[ii]][[i]]
        dat <- dat[order(dat$subN), ]
      } else {
        dat <- object[[i]][order(object[[i]]$subN), ]
      }
      if(!first){
        css <- dat[dat[, ci] > cor, 'drop']
        css <- ifelse(length(css) == 0, 0, max(css))
      } else {
        css <- which(!dat[, ci] > cor)
        css <- dat[ifelse(length(css) == 0, 1, max(css) + 1), 'drop']
      }
      return(css)
    }), inds)
  }
  if(lags){
    call <- replace(as.list(match.call())[-1], c('object', 'verbose'),
                    list(object = obj0$contemporaneous, verbose = FALSE))
    out <- list(temporal = out, contemporaneous = do.call(cscoef, call))
  }
  if(verbose){
    cat(paste0('CS = ', cor, '(', ci0, '%)'), '\n')
    print(out)
  } else {
    attr(out, c('cor', 'ci', 'first')) <- c(cor, ci0/100, first)
    return(out)
  }
}
