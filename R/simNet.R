#' Simulate network structure and data
#'
#' Used for generating moderated and unmoderated adjacency matrices, along with
#' data based on those model structures.
#'
#' If no moderator is specified then data can be generated directly from a
#' partial correlation matrix by setting \code{gibbs = FALSE}, which produces
#' fast simulation results. Alternatively, a Gibbs sampler is used to generate
#' data, which is the default option. For moderated networks, Gibbs sampling is
#' the only method available.
#'
#' @section Warning:
#'
#'   Importantly, the Gibbs sampler can easily diverge given certain model
#'   parameters. Generating network data based on moderator variables can
#'   produce data that quickly take on large values due to the presence of
#'   multiplicative terms. If the simulation fails, first simply try re-running
#'   the function with a different seed; this will often be sufficient to solve
#'   the problem when default parameters are specified. Additionally, one can
#'   increase the value of \code{div}, in case the sampler only diverges
#'   slightly or simply produced an anomalous value. This raises the threshold
#'   of tolerated values before the sampler stops. If supplying user-generated
#'   model matrices (for the \code{b1} and/or \code{b2} arguments) and the
#'   function continues to fail, you will likely need to change the parameter
#'   values in those matrices, as it may not be possible to simulate data under
#'   the given values. If simulating the model matrices inside the function (as
#'   is the default) and the function continues to fail, try adjusting the
#'   following parameters: \enumerate{ \item{Try reducing the value of \code{m2}
#'   to specify fewer interactions}. \item{Try reducing a range with a smaller
#'   maximum for \code{m2_range}, to adjust the range of interaction
#'   coefficients}. \item{Try adjusting the corresponding main effect parameters
#'   for the moderator, \code{m1} and \code{m1_range}}. \item{Try setting
#'   \code{modType = "full"} to reduce the number of main effect parameters}.
#'   \item{Try setting a low value(s) for \code{fixedPar}, in order to provide
#'   parameter values that are known to be lower} }
#'
#'   An alternative approach could be to use the internal function
#'   \code{simNet2}, which is a wrapper designed to re-run \code{simNet} when it
#'   fails and automatically adjust simulation parameters such as \code{div} to
#'   thoroughly test a given parameterization scheme. This function can be
#'   accessed via \code{modnets:::simNet2}. There is not documentation for this
#'   function, so it is recommended to look at the source code if one wishes to
#'   use it This wrapper is also used inside the \code{mnetPowerSim} function.
#'
#' @param N Numeric value. Total number of subjects.
#' @param p Numeric value. Total number of nodes (excluding moderator).
#' @param m If a value is provided, a moderator is generated and named \code{M}
#'   in the resultant data. If \code{TRUE}, then a normal distribution with a
#'   mean of 0 will be used to generate the initial value of \code{m}, which
#'   will serve as the population mean for \code{m} throughout the simulation.
#'   If a numeric value is provided, then this will serve as the population
#'   mean, and all subsequent draws will be taken from a normal distribution
#'   with that mean. If \code{m = "binary"}, then this will simply set the
#'   argument \code{mbinary = TRUE}. If \code{m = "ordinal"}, this will set
#'   \code{mord = TRUE}. To simulate \code{m} from a skewed distribution, there
#'   are two options: if \code{m = "skewed"}, then the \code{alpha} parameter of
#'   the \code{\link[sn:rmsn]{sn::rmsn}} will automatically be set to 3.
#'   Alternatively, a vector of length two can be supplied, containing the
#'   element \code{"skewed"} as well as the desired value of \code{alpha}.
#'   Lastly, a function can be provided for \code{m} if the user wishes to
#'   sample \code{m} from another distribution. The requirement is that the
#'   function have only one argument, and only returns a single numeric value.
#'   The input of the argument should be the location parameter of the desired
#'   sampling distribution.
#' @param m2 Numeric. If \code{m2 >= 1}, then this will determine the number of
#'   interaction effects between the moderator and some node in the network. If
#'   a value between 0 and 1 is provided, then this determines the probability
#'   of any given edge being moderated by the moderator.
#' @param b1 Can provide an adjacency matrix to use for generating data.
#' @param b2 Can provide an interaction matrix for generated moderated data.
#' @param sparsity Numeric value between 0 and 1. Determines the sparsity of
#'   sampled network matrices.
#' @param intercepts A vector of means for sampling node values.
#' @param nIter Number of iterations for generating each instance of a datapoint
#'   with the Gibbs sampler.
#' @param msym If \code{TRUE} then will force the interaction matrix to be
#'   symmetric.
#' @param onlyDat If \code{TRUE} then the function only returns the simulated
#'   data.
#' @param pbar If \code{TRUE} then a progress bar will be shown as samples are
#'   generated.
#' @param div A value to use as a sign that the sampler diverged. Can be
#'   increased based on expected range of values. If a datapoint is larger than
#'   \code{div}, then the sampler will stop.
#' @param gibbs If \code{TRUE}, then Gibbs sampling will be used. Otherwise,
#'   data are generated from the \code{\link[mvtnorm:rmvnorm]{mvtnorm::rmvnorm}}
#'   function based on the partial correlation matrix that is created.
#' @param ordinal Logical. Determines whether to generate ordinal values or not.
#' @param nLevels Number of levels for the ordinal variables. Only relevant if
#'   \code{ordinal} is not \code{FALSE}.
#' @param mord Logical. Determines whether the moderator variable should be
#'   simulated as ordinal.
#' @param time If \code{TRUE} then the time it takes to simulate the data is
#'   printed to screen at the end of the sampling.
#' @param mbinary Logical. Determines whether the moderator should be a binary
#'   variable.
#' @param minOrd The minimum number of unique values allowed for each variable.
#' @param m1 Functions similarly to \code{m2}, except that this argument refers
#'   to the number/probability of main effects of the moderator on any given
#'   node.
#' @param m1_range Numeric vector of length 2. The range of values for moderator
#'   main effect coefficients.
#' @param m2_range Numeric vector of length 2. The range of values for moderator
#'   interaction effect coefficients.
#' @param modType Determines the type of moderation to employ, such as
#'   \code{"none", "full", "partial"}. If \code{modType = "full"}, then for any
#'   interaction terms there will be full moderation, such that all pairwise
#'   relationships for moderated paths will be set to zero. If \code{modType =
#'   "partial"}, then pairwise edges for moderated paths will always be nonzero.
#'   If \code{modType = "none"}, no constraints will be applied (e.g., could
#'   produce a mix between full and partial moderation).
#' @param lags If \code{TRUE} or 1, then arguments are rerouted to the
#'   \code{\link{mlGVARsim}} function to simulate temporal data for a single
#'   individual.
#' @param V Numeric, either 1 or 2. Determines whether to randomize the order of
#'   simulating node values at each iteration of the Gibbs sampler. If \code{V =
#'   2}, then the order is randomized at each iteration. If \code{V = 1}, then
#'   the sampler moves through the nodes from the first to the last in order at
#'   each iteration.
#' @param skewErr The skewness parameter for the \code{alpha} argument in the
#'   \code{\link[sn:rmsn]{sn::rmsn}} function. Only relevant when \code{gibbs =
#'   FALSE} and no moderator is specified.
#' @param onlyNets If \code{TRUE} then only the network models are returned,
#'   without the data. Could be used to create random models and then simulate
#'   data by another method.
#' @param netArgs Only for use by the internal function
#'   \code{modnets:::simNet2}, which serves as a wrapper for the current
#'   function to prevent it from failing.
#' @param nCores Numeric value indicating the number of CPU cores to use for the
#'   resampling. If \code{TRUE}, then the
#'   \code{\link[parallel:detectCores]{parallel::detectCores}} function will be
#'   used to maximize the number of cores available.
#' @param cluster Character vector indicating which type of parallelization to
#'   use, if \code{nCores > 1}. Options include \code{"mclapply"} and
#'   \code{"SOCK"}.
#' @param getChains Logical. Determines whether to return the data-generating
#'   chains from the Gibbs sampler.
#' @param const Numeric. The constant to be used by the internal
#'   \code{modnets:::simPcor} function.
#' @param fixedPar Numeric. If provided, then this will be set as the
#'   coefficient value for all edges in the network. Provides a way to
#'   standardize the parameter values while varying the sparsity of the network.
#'   If \code{length(fixedPar) == 1}, then the same value will be used for all
#'   parameters. If \code{length(fixedPar) == 2}, then the first value will be
#'   for pairwise relationships, and the second value will be for interaction
#'   terms.
#' @param V2 If \code{V2 = 1} and \code{m2} is between 0 and 1, the number of
#'   interaction terms in the model will be determined by multiplying \code{m2}
#'   with the number of elements in the interaction matrix and taking the
#'   \code{ceiling}.
#' @param ... Additional arguments.
#'
#' @return Simulated network models as well as data generated from those models.
#'   For GGMs, model matrices are always symmetric. For temporal networks (when
#'   \code{lags = 1}), columns predict rows.
#' @export
#'
#' @seealso \code{\link{mlGVARsim}, \link{mnetPowerSim}, \link{plotNet},
#'   \link{net}, \link{netInts}, \link{plotBoot}, \link{plotCoefs}}
#'
#' @examples
#'
#' # Generate a moderated GGM along with data
#' set.seed(1)
#' x <- simNet(N = 100, p = 3, m = TRUE)
#'
#' net(x) # Get data-generating adjacency matrix
#' netInts(x) # Get data-generating interaction matrix
#'
#' plot(x) # Plot the moderated network that generated the data
#'
#' # Generate a single-subject GVAR model with data
#' set.seed(1)
#' x <- simNet(N = 500, p = 3, m = TRUE, lags = 1)
#'
#' net(x, n = 'temporal') # Get the data-generating time-lagged adjacency matrix
#' net(x, n = 'contemporaneous') # Get the data-generating standardized residual covariance matrix
#'
#' plot(x, which.net = 'beta') # 'beta' is another way of referring to the temporal network
#' plot(x, which.net = 'pcc') # 'pcc' is another way of referring to the contemporaneous network
simNet <- function(N = 100, p = 5, m = FALSE, m2 = .1, b1 = NULL, b2 = NULL,
                   sparsity = .5, intercepts = NULL, nIter = 250, msym = FALSE,
                   onlyDat = FALSE, pbar = TRUE, div = 10, gibbs = TRUE,
                   ordinal = FALSE, nLevels = 5, mord = FALSE, time = TRUE,
                   mbinary = FALSE, minOrd = 3, m1 = NULL, m1_range = NULL,
                   m2_range = c(.1, .3), modType = 'none', lags = NULL, V = 2,
                   skewErr = FALSE, onlyNets = FALSE, netArgs = NULL,
                   nCores = 1, cluster = 'SOCK', getChains = FALSE,
                   const = 1.5, fixedPar = NULL, V2 = 1, ...){
  t1 <- Sys.time()
  args <- tryCatch({list(...)}, error = function(e){list()})
  if(!is.null(netArgs)){list2env(netArgs[setdiff(names(netArgs), 'data')], envir = environment())}
  if(!is.null(lags) & !identical(as.numeric(lags), 0)){
    args1 <- append(list(nTime = N, nPerson = 1, nNode = p, lag = 1, GGMsparsity = sparsity), as.list(match.call())[-1])
    if(!identical(m, FALSE)){args1 <- replace(args1, 'm', ifelse(mbinary, 'binary', ifelse(mord, 'ordinal', m)))}
    if('nPerson' %in% names(args)){args1 <- replace(args1, 'nPerson', args$nPerson)}
    args1 <- append(args1, args[setdiff(names(args), names(args1))])
    out <- do.call(mlGVARsim, args1[intersect(names(args1), formalArgs('mlGVARsim'))])
    return(out)
  }
  skewErr <- ifelse(is.numeric(skewErr), skewErr, ifelse(!identical(skewErr, FALSE), 3, FALSE))
  gibbsFun <- switch(2 - identical(skewErr, FALSE), function(x){rnorm(1, x, 1)},
                     function(x){sn::rsn(1, x, alpha = skewErr)[1]})
  if(!is.null(m1_range)){if(!is.numeric(m1_range) | length(m1_range) != 2){m1_range <- NULL}}
  if(!is.numeric(m2_range) | length(m2_range) != 2){m2_range <- c(.1, .3)}
  modType <- match.arg(tolower(modType), c('none', 'full', 'partial', 'full2', 'partial2', 'zero'))
  if(grepl('2', modType) & is.null(m1)){m1 <- .5}
  if(is.null(intercepts)){
    intercepts <- rep(0, p)
  } else if('skewed' %in% intercepts){
    if(length(intercepts) == 1 | length(intercepts) > 2){intercepts <- '3'}
    intercepts <- replicate(p, sn::rsn(1, alpha = as.numeric(setdiff(intercepts, 'skewed'))))
  } else if(length(intercepts) != p){
    intercepts <- rnorm(p)
  }
  if(is.function(m)){
    mfun <- m
    m <- mfun()
  } else {
    if(is.character(m)){
      if('binary' %in% m){mbinary <- TRUE}
      if('ordinal' %in% m){mord <- TRUE}
      if(!'skewed' %in% m){m <- TRUE}
    }
    mfun <- switch(
      2 - mbinary,
      function(size = 1){sample(x = 0:1, size = size)},
      function(mean = 0){rnorm(n = 1, mean = mean)}
    )
  }
  if(isTRUE(m)){m <- ifelse(mbinary, 1, mfun())}
  if('skewed' %in% m){
    if(length(m) == 1 | length(m) > 2){m <- '3'}
    mskew <- as.numeric(setdiff(m, 'skewed'))
    mfun <- function(mean = 0){sn::rsn(n = 1, xi = mean, alpha = mskew)[1]}
    m <- mfun()
  }
  m <- ifelse(identical(m, FALSE), 0, m)
  if(is.null(b1)){
    b1 <- simPcor(p, sparsity, constant = const, finalRange = m1_range)
    if(!is.null(fixedPar)){b1[b1 != 0] <- fixedPar[1] * sign(b1[b1 != 0])}
  }
  if(is.null(b2)){
    mat <- b2 <- diag(0, p)
    pn <- (p * (p - 1))/2
    if(V2 == 1 & m2 < 1){m2 <- ceiling(m2 * pn)}
    while(all(b2 == 0) & modType != 'zero'){
      if(m2 >= 0 & m2 < 1){
        b2 <- sample(0:1, pn, TRUE, prob = c(1 - m2, m2))
      } else if(m2 >= 1){
        b2 <- numeric(pn)
        b2[sample(1:pn, ifelse(m2 > pn, pn, round(m2)))] <- 1
      }
      b2 <- b2 * sample(c(-1, 1), pn, TRUE, prob = c(.5, .5))
      mat[upper.tri(mat)] <- b2 * runif(pn, min(m2_range), max(m2_range))
      b2 <- as.matrix(Matrix::forceSymmetric(mat))
      if(msym){b2 <- simPcor(p, sparsity = 1 - m2, constant = const, finalRange = m2_range)}
      if(!is.null(fixedPar)){
        b2[b2 != 0] <- switch(length(fixedPar), fixedPar, fixedPar[2]) * sign(b2[b2 != 0])
      }
    }
    if(all(m == 0) | modType == 'zero'){b2 <- diag(0, p)}
  }
  if(!modType %in% c('none', 'zero') & !all(b2 == 0)){
    while(ifelse(grepl('full', modType), !all(b1[b2 != 0] == 0), any(b1[b2 != 0] == 0))){
      b1 <- simPcor(p, sparsity, constant = const, finalRange = m1_range)
      if(!is.null(fixedPar)){b1[b1 != 0] <- fixedPar[1] * sign(b1[b1 != 0])}
    }
  }
  diag(b1) <- diag(b2) <- 0
  if(is.null(m1) | identical(m1, 0) | all(m == 0)){
    m1 <- rep(0, p)
  } else {
    if(isTRUE(m1)){m1 <- .5}
    if(length(m1) == 1){
      if(is.null(m1_range)){m1_range <- c(.1, .4)}
      m10 <- runif(p, min(m1_range), max(m1_range)) * sample(c(-1, 1), p, TRUE, prob = c(.5, .5))
      if(m1 >= 0 & m1 < 1){
        m10 <- m10 * sample(0:1, p, TRUE, prob = c(1 - m1, m1))
      } else if(m1 >= 1){
        m10[sample(1:p, ifelse(m1 > p, 0, round(p - m1)))] <- 0
      }
      if(grepl('2', modType) & !all(b2 == 0)){
        mm01 <- apply(b2, 1, function(z) any(z != 0))
        while(ifelse(grepl('full', modType), !all(m10[mm01] == 0), any(m10[mm01] == 0))){
          m10 <- runif(p, min(m1_range), max(m1_range))
          if(m1 >= 0 & m1 < 1){
            m10 <- m10 * sample(0:1, p, TRUE, prob = c(1 - m1, m1))
          } else if(m1 >= 1){
            m10[sample(1:p, ifelse(m1 > p, 0, round(p - m1)))] <- 0
          }
        }
      }
      m1 <- m10
      if(!is.null(fixedPar) & !all(m1 == 0)){
        m1[m1 != 0] <- fixedPar[1] * sign(m1[m1 != 0])
      }
    }
  }
  if(onlyNets & !onlyDat){
    out <- structure(list(
      b1 = b1, b2 = b2, intercepts = intercepts, m = m, m1 = m1),
      m2 = m2, modType = modType, m1_range = m1_range, m2_range = m2_range)
    return(out)
  }
  if(nCores > 1 | isTRUE(nCores)){
    if(isTRUE(nCores)){nCores <- parallel::detectCores()}
    if(grepl('Windows', sessionInfo()$running)){cluster <- 'SOCK'}
    if(tolower(cluster) != 'mclapply'){
      cluster <- match.arg(toupper(cluster), c('SOCK', 'FORK'))
      cl <- parallel::makeCluster(nCores, type = cluster)
    } else {
      cl <- nCores
    }
    sampFun <- function(case, nIter, p, m, b1, b2, intercepts, div,
                        V, m1, skewErr, getChains){
      sampling <- matrix(NA, nrow = nIter, ncol = p)
      sampling[1, ] <- rnorm(p)
      m0 <- c(ifelse(m == 0, 0, mfun(m)), numeric(nIter - 1))
      for(iter in 2:nIter){
        m0[iter] <- ifelse(m == 0, 0, mfun(m))
        vv <- switch(V, 1:p, sample(1:p, p, replace = FALSE))
        for(v in vv){
          v_mu <- 0
          v_ps <- which(b1[v, ] != 0)
          if(length(v_ps) > 0){
            for(vp in v_ps){
              v_it <- iter - as.numeric(!vp %in% vv[1:which(vv == v)])
              v_mu <- c(v_mu, sampling[v_it, vp] * b1[v, vp])
            }
          }
          if(m1[v] != 0){v_mu <- c(v_mu, m0[iter] * m1[v])}
          v_ps2 <- which(b2[v, ] != 0)
          if(length(v_ps2) > 0){
            for(vp in v_ps2){
              v_it <- iter - as.numeric(!vp %in% vv[1:which(vv == v)])
              v_mu <- c(v_mu, ((sampling[v_it, vp] * m0[iter]) * b2[v, vp]))
            }
          }
          v_mu <- intercepts[v] + sum(v_mu)
          sampling[iter, v] <- gibbsFun(v_mu)
          if(any(abs(sampling[!is.na(sampling)]) > div)){stop('Sampler diverged')}
        }
      }
      out <- list(data = sampling[nIter, ], M = m0[nIter])
      if(getChains){
        out$chains <- structure(cbind(sampling, m0), dimnames = NULL)
        if(all(m == 0)){out$chains <- out$chains[, -ncol(out$chains)]}
      }
      return(out)
    }
    if(cluster == 'SOCK'){
      parallel::clusterExport(cl, c('mfun', 'sampFun', 'gibbsFun'), envir = environment())
    }
    out <- tryCatch({
      if(pbar){
        pbapply::pboptions(type = 'timer', char = '-')
        pbapply::pblapply(1:N, sampFun, nIter = nIter, p = p, m = m, b1 = b1,
                          b2 = b2, intercept = intercepts, div = div, V = V,
                          m1 = m1, skewErr = skewErr, getChains = getChains, cl = cl)
      } else if(tolower(cluster) != 'mclapply'){
        parallel::parLapply(
          cl, 1:N, sampFun, nIter = nIter, p = p, m = m, b1 = b1, b2 = b2,
          intercepts = intercepts, div = div, V = V, m1 = m1,
          skewErr = skewErr, getChains = getChains)
      } else {
        parallel::mclapply(
          1:N, sampFun, nIter = nIter, p = p, m = m, b1 = b1, b2 = b2,
          intercepts = intercepts, div = div, V = V, m1 = m1,
          skewErr = skewErr, getChains = getChains, mc.cores = nCores)
      }},
      error = function(e){TRUE}
    )
    if(tolower(cluster) != 'mclapply'){parallel::stopCluster(cl)}
    if(isTRUE(out)){stop('Sampler diverged')}
    data <- do.call(rbind, lapply(out, '[[', 'data'))
    M <- unlist(lapply(out, '[[', 'M'))
    if(getChains){chains <- do.call(abind::abind, c(lapply(out, '[[', 'chains'), along = 0))}
    rm(cl)
  } else if(isTRUE(gibbs) | !all(m == 0)){
    data <- matrix(NA, nrow = N, ncol = p); M <- c()
    chains <- array(NA, dim = c(N, nIter, p + as.numeric(all(m != 0))))
    if(pbar){pb <- txtProgressBar(max = N, style = 3)}
    for(case in 1:N){
      sampling <- matrix(NA, nrow = nIter, ncol = p)
      sampling[1, ] <- rnorm(p)
      m0 <- c(ifelse(m == 0, 0, mfun(m)), numeric(nIter - 1))
      for(iter in 2:nIter){
        m0[iter] <- ifelse(m == 0, 0, mfun(m))
        vv <- switch(V, 1:p, sample(1:p, p, replace = FALSE))
        for(v in vv){
          v_mu <- 0
          v_ps <- which(b1[v, ] != 0)
          if(length(v_ps) > 0){
            for(vp in v_ps){
              v_it <- iter - as.numeric(!vp %in% vv[1:which(vv == v)])
              v_mu <- c(v_mu, sampling[v_it, vp] * b1[v, vp])
            }
          }
          if(m1[v] != 0){v_mu <- c(v_mu, m0[iter] * m1[v])}
          v_ps2 <- which(b2[v, ] != 0)
          if(length(v_ps2) > 0){
            for(vp in v_ps2){
              v_it <- iter - as.numeric(!vp %in% vv[1:which(vv == v)])
              v_mu <- c(v_mu, ((sampling[v_it, vp] * m0[iter]) * b2[v, vp]))
            }
          }
          v_mu <- intercepts[v] + sum(v_mu)
          sampling[iter, v] <- gibbsFun(v_mu)
          if(any(abs(sampling[!is.na(sampling)]) > div)){stop('Sampler diverged')}
        }
      }
      if(getChains){
        chains[case, , ] <- structure(switch(
          2 - all(m == 0), sampling, cbind(sampling, m0)),
          dimnames = NULL)
      }
      data[case, ] <- sampling[nIter, ]
      M[case] <- m0[nIter]
      if(pbar){
        setTxtProgressBar(pb, case)
        if(case == N){close(pb)}
      }
    }
  } else {
    Sigma <- cov2cor(solve(diag(p) - b1))
    if(identical(skewErr, FALSE)){
      data <- mvtnorm::rmvnorm(n = N, mean = intercepts, sigma = Sigma)
    } else {
      data <- sn::rmsn(n = N, xi = intercepts, Omega = Sigma, alpha = rep(skewErr, p))
    }
  }
  out <- list(b1 = b1, b2 = b2, intercepts = intercepts)
  if(ordinal){
    for(i in 1:ncol(data)){
      ord <- c()
      while(length(unique(ord)) < minOrd){
        ord <- as.numeric(cut(data[, i], sort(c(-Inf, rnorm(nLevels - 1), Inf))))
      }
      data[, i] <- ord
    }
  }
  if(all(m == 0)){
    out$b2 <- NULL
  } else {
    if(mord & !mbinary){
      if(minOrd < 2){minOrd <- 2}
      while(length(unique(mord)) < minOrd){
        mord <- as.numeric(cut(M, sort(c(-Inf, rnorm(nLevels - 1), Inf))))
      }
      M <- mord
    }
    data <- cbind(data, M)
  }
  out <- append(list(data = data.frame(data)), out)
  if(onlyDat){out <- data.frame(data)} else if(getChains){out$chains <- chains}
  if(!all(m == 0)){
    if(onlyDat){
      attributes(out)[c('m', 'm1')] <- list(m = m, m1 = m1)
    } else {
      out <- append(out, list(m = m, m1 = m1))
    }
    attributes(out)[c('m2', 'modType')] <- list(m2, modType)
  }
  class(out) <- c(ifelse(onlyDat, 'data.frame', 'list'), 'ggmSim')
  attr(out, 'time') <- t2 <- Sys.time() - t1
  if(time){print(Sys.time() - t1)}
  return(out)
}

##### simNet2: wrapper for simNet to prevent failures
simNet2 <- function(..., nets = NULL, maxiter = 10, pbar = TRUE, div = 10,
                    divup = TRUE, maxdiv = 10^3, errors = FALSE, nostop = FALSE){
  args <- tryCatch({list(...)}, error = function(e){list()})
  defaults <- list(pbar = pbar, div = div, time = FALSE, onlyDat = FALSE)
  args2 <- args <- replace(args, names(defaults), defaults)
  # This was the adjustment made...
  if(!'lags' %in% names(args)){
    args$onlyNets <- TRUE
    args2$netArgs <- switch(2 - is.null(nets), do.call(simNet, args), nets)
  }
  tt <- t <- out <- 0
  div0 <- div
  while(identical(out, 0)){
    out <- tryCatch({do.call(simNet, args2)}, error = function(e){0})
    t <- t + 1
    if(t == 2){
      if(!is.null(nets)){
        if(divup & div < maxdiv){
          args$div <- args2$div <- div <- div * div0
          t <- 0
        } else {
          stop('Need to generate new network matrices')
        }
      } else {
        tt <- tt + 1; t <- 0
        if(!identical(errors, FALSE)){
          if(tt == 1 & div == div0){errors <- list()}
          errors <- append(errors, setNames(
            list(args2$netArgs), paste0('tt', tt, '_div', div)))
        }
        if(tt == maxiter){
          if(divup & div < maxdiv){
            args$div <- args2$div <- div <- div * div0
            tt <- 0
          } else {
            if(nostop){
              return(list())
            } else {
              stop('Parameters may be intractable')
            }
          }
        }
        args2$netArgs <- do.call(simNet, args)
      }
    }
  }
  attributes(out)[c('div', 't', 'tt')] <- list(div, t, tt)
  if(!identical(errors, FALSE)){attributes(out)$errors <- errors}
  return(out)
}
