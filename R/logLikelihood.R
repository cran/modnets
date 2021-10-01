#' Log-likelihood functions and Likelihood Ratio Tests for moderated networks
#'
#' Computes log-likelihood, AIC, and BIC for a whole network, or for each node
#' in the network. Also compares two or more networks using a likelihood ratio
#' test (LRT).
#'
#' Fits LRT to a list of network models to compare them all against each other.
#' Obtain all possible LRTs comparing a list of SUR models. Can include tests
#' comparing RMSEA values. The \code{nodes} argument determines whether to
#' perform these computations in an omnibus or nodewise fashion.
#'
#' One key thing to note is that when using \code{\link{modTable}} or
#' \code{\link{SURtable}}, the LRT column indicates the number of times that
#' each network was selected over others with respect to the pairwise LRTs.
#'
#' @param net0 Output from one of the main \code{modnets} functions.
#'   Alternatively, a list of network models can be provided. This list should
#'   be named for the easiest interpretation of results.
#' @param net1 For \code{\link{modLL}} and \code{\link{SURll}}, can be used to
#'   supply a second network to compare to \code{net0} via an LRT. Or if
#'   \code{lrt = FALSE}, then relevant statistics will be returned for both
#'   \code{net0} and \code{net1}. Importantly, if one network is provided for
#'   \code{net0}, and another is provided for \code{net1}, then the names in the
#'   output will reflect these argument names. This can be somewhat confusing at
#'   times, so ultimately it is not recommended to use this argument. Instead,
#'   try supplying both networks (or more) as a named list to the \code{net0}
#'   argument for the most customization.
#' @param nodes Logical. Determines whether to compute omnibus or nodewise
#'   statistics and tests. If \code{TRUE}, then LL values for nodewise models
#'   will be returned, and any LRTs requested will reflect nodewise tests.
#' @param lrt Logical. Determines whether to conduct an LRT or not. If
#'   \code{FALSE}, then only LL-related statistics will be returned for all
#'   models supplied. Only relevant for \code{\link{modLL}} and
#'   \code{\link{SURll}}
#' @param all Logical. If \code{TRUE}, then omnibus LL statistics as well as
#'   nodewise statistics are returned for either LL function.
#' @param d Number of decimal places to round outputted statistics to.
#' @param alpha Alpha level for LRTs. Defaults to .05.
#' @param orderBy Can be one of \code{"LL", "df", "AIC", "BIC"} to indicate
#'   which statistic to order the table by. If using \code{\link{modTable}} or
#'   \code{\link{SURtable}}, then a value of \code{TRUE} will organize the
#'   output by the LRT column, which indicates the number of times that a
#'   particular model performed better than another based on the pairwise LRTs.
#'   Higher values indicate that the model was selected more often in comparison
#'   with other models in the list.
#' @param decreasing Logical. Determines whether to organize output from highest
#'   to lowest, or vice versa, in accordance with the value of \code{orderBy}.
#' @param s Character string indicating which type of residual covariance matrix
#'   to compute for SUR models. Options include \code{"res", "dfres", "sigma"}.
#'   \code{"sigma"} uses the residual covariance matrix as computed by the
#'   \code{\link[systemfit:systemfit]{systemfit::systemfit}} function.
#'   \code{"res"} and \code{"dfres"} compute the matrix based directly on the
#'   residual values. \code{"dfres"} is the sample estimator that uses \code{N -
#'   1} in the denominator, while \code{"res"} just uses \code{N}.
#' @param sysfits Logical, only relevant to \code{\link{SURll}} when multiple
#'   networks are included in a list. Does not currently work when there are two
#'   networks in the list, but does work with 3 or more. Returns the omnibus
#'   model statistics based on functions available to output from the
#'   \code{\link[systemfit:systemfit]{systemfit::systemfit}} function. This
#'   allows for some additional statistics such as \code{SSR, detSigma, OLS.R2,
#'   McElroy.R2}.
#' @param names Character vector containing the names of the models being
#'   compared. Only relevant to the \code{\link{modTable}} and
#'   \code{\link{SURtable}}. Alternatively, models can be named by supplying a
#'   named list to the \code{net0} argument.
#' @param rmsea Logical. Relevant to \code{\link{modTable}} and
#'   \code{\link{SURtable}}. Determines whether to return RMSEA values, as well
#'   as tests comparing RMSEA across each pair of models.
#'
#' @return A table or list of results depending on which function is used.
#' @export
#' @name LogLikelihood
#'
#' @seealso \code{\link{fitNetwork}, \link{mlGVAR}}
#'
#' @examples
#' data <- na.omit(psychTools::msq[, c('hostile', 'lonely', 'nervous', 'sleepy', 'depressed')])
#'
#' ##### Use modLL() for GGMs
#' ggm1 <- fitNetwork(data[, -5])
#' ggm2 <- fitNetwork(data, covariates = 5)
#' ggm3 <- fitNetwork(data, moderators = 5)
#'
#' modLL(ggm1)
#' modLL(ggm2)
#'
#' modLL(ggm1, ggm2)
#' modLL(ggm1, ggm2, nodes = TRUE)
#'
#' modLL(list(ggm1 = ggm1, ggm2 = ggm2))
#' modLL(list(GGM1 = ggm1, GGM2 = ggm2), nodes = TRUE)
#'
#' ggms <- list(ggm1, ggm2, ggm3)
#'
#' modLL(ggms)
#' modTable(ggms)
#' modTable(ggms, names = c("GGM1", "GGM2", "GGM3"))
#'
#' names(ggms) <- c("GGM1", "GGM2", "GGM3")
#' modTable(ggms)
#' modLL(ggms)
#'
#' ##### Use SURll() for SUR networks
#' sur1 <- fitNetwork(data[, -5], lags = TRUE)
#' sur2 <- fitNetwork(data, covariates = 5, lags = TRUE)
#' sur3 <- fitNetwork(data, moderators = 5, lags = TRUE)
#'
#' SURll(sur1)
#' SURll(sur2)
#'
#' SURll(sur1, sur2)
#' SURll(sur1, sur2, nodes = TRUE)
#' SURll(list(SUR1 = sur1, SUR2 = sur2), nodes = TRUE)
#'
#' surs <- list(sur1, sur2, sur3)
#'
#' SURll(surs)
#' SURtable(surs, names = c('SUR1', "SUR2", "SUR3"))
#'
#' names(surs) <- c("SUR1", "SUR2", "SUR3")
#' SURll(surs)
#' SURtable(surs)
modLL <- function(net0, net1 = NULL, nodes = FALSE, lrt = NULL, all = FALSE,
                  d = 4, alpha = .05, orderBy = NULL, decreasing = TRUE){
  # Log-likelihood for whole network
  mll <- function(object){
    k <- length(object$fitobj)
    res <- as.matrix(do.call(cbind, lapply(object$fitobj, resid)))
    n <- nrow(res)
    sigma <- (t(res) %*% res)/n
    inv <- solve(sigma)
    ll <- sum(((-k/2) * log(2 * pi)) - (.5 * log(det(sigma))) - (.5 * diag(res %*% inv %*% t(res))))
    df <- sum(sapply(lapply(object$mods, '[[', "model"), nrow)) + 1
    return(c(LL = ll, df = df))
  }

  # Log-likelihood and information criteria for single regression
  uni_mll <- function(object){
    getInd <- function(ind, x){
      ind <- c("deviance", "LL_model", "df.residual", "AIC", "BIC")[ind]
      X <- ifelse(ind == "df.residual", list(x$fitobj), list(x$mods))[[1]]
      return(unname(sapply(X, '[[', ind)))
    }
    out <- do.call(cbind.data.frame, lapply(1:5, getInd, x = object))
    colnames(out) <- c("RSS", "LL", "df", "AIC", "BIC")
    rownames(out) <- names(object$mods)
    out
  }

  # Combines mll with information criteria for full network
  omni_mll <- function(object){
    k <- length(object$fitobj)
    n <- nrow(object$data) * k
    ll <- mll(object = object)
    aic <- (2 * ll[2]) - (2 * ll[1])
    bic <- (ll[2] * log(n)) - (2 * ll[1])
    out <- c(ll, AIC = unname(aic), BIC = unname(bic))
    out
  }

  # Likelihood ratio test for comparing two networks
  mod_lrt <- function(object, d = 4, alpha = .05, N = NULL){
    if(is.list(object)){
      if(length(object) > 2){object <- object[1:2]}
      nn <- names(object)
      ll0 <- object[[1]]$LL; df0 <- object[[1]]$df
      ll1 <- object[[2]]$LL; df1 <- object[[2]]$df
      omnibus <- FALSE
    } else {
      if(nrow(object) > 2){object <- object[1:2, ]}
      nn <- rownames(object)
      ll0 <- object[1, 1]; df0 <- object[1, 2]
      ll1 <- object[2, 1]; df1 <- object[2, 2]
      omnibus <- TRUE
    }
    lldiff <- abs(ll0 - ll1) * 2
    dfdiff <- abs(df0 - df1)
    ps <- pchisq(q = lldiff, df = dfdiff, lower.tail = FALSE)
    decision <- c()
    for(i in seq_along(ps)){
      if(ps[i] <= alpha){
        decision[i] <- ifelse(ll0[i] > ll1[i], nn[1], nn[2])
      } else if(ps[i] == 1){
        decision[i] <- "- "
      } else if(ps[i] > alpha){
        if(omnibus){
          decision[i] <- ifelse(df0 < df1, nn[1], nn[2])
        } else {
          decision[i] <- ifelse(df0[i] > df1[i], nn[1], nn[2])
        }
      }
    }
    if(!omnibus){
      if(!is.null(d)){ps <- round(ps, d)}
      out <- data.frame(LL_diff2 = lldiff, Df_diff = dfdiff, pval = ps, decision = decision)
      rownames(out) <- rownames(object[[1]])
    } else {
      RMSEA <- function(X2, df, N){
        rmsea <- sqrt(max(c(((X2/N)/df) - (1/N), 0)))
        lower.l <- function(lambda){(pchisq(q = X2, df = df, ncp = lambda) - .95)}
        lambda.l <- tryCatch({uniroot(f = lower.l, lower = 0, upper = X2)$root}, error = function(e){0})
        rmsea.lower <- sqrt(lambda.l/(N * df))
        upper.l <- function(lambda){(pchisq(q = X2, df = df, ncp = lambda) - .05)}
        lambda.u <- tryCatch({uniroot(f = upper.l, lower = 0, upper = max(N, X2 * 4))$root}, error = function(e){1})
        rmsea.upper <- sqrt(lambda.u/(N * df))
        rmsea.pvalue <- 1 - pchisq(q = X2, df = df, ncp = (N * df * (.05^2)))
        return(c(lower = rmsea.lower, RMSEA = rmsea, upper = rmsea.upper, p.value = rmsea.pvalue))
      }
      rmsea <- RMSEA(lldiff, dfdiff, N)
      if(!is.null(d)){
        ps <- round(ps, d)
        lldiff <- round(lldiff, d)
        rmsea <- round(rmsea, d)
      }
      if(object[2, 2] < object[1, 2]){object <- object[order(object[, 2]), ]}
      out0 <- data.frame(LL_diff2 = c("", lldiff), Df_diff = c("", dfdiff),
                         pval = c("", ps), decision = c("", decision))
      out <- data.frame(object[, 1:2], out0, object[, 3:4])
      attr(out, "RMSEA") <- rmsea
    }
    return(out)
  }

  # Helps simplify interface; you can put two networks in a list for net0
  nn <- paste0("net", 0:1)
  if(length(net0) == 2 & ifelse(
    is.null(net1), TRUE, ifelse(is.logical(net1), net1, FALSE))){
    if(is.logical(net1)){nodes <- net1}
    nn <- ifelse(!is.null(names(net0)), list(names(net0)), list(nn))[[1]]
    net1 <- net0[[2]]; net0 <- net0[[1]]
  }

  # Take the between-subjects network if mlGVAR
  if(isTRUE(attr(net0, "mlGVAR"))){net0 <- net0$betweenNet}

  # If we have a ggm...
  if("ggm" %in% names(attributes(net0))){
    omni0 <- switch(2 - isTRUE('modLL' %in% names(net0)), net0$modLL$omnibus, omni_mll(net0))
    uni0 <- switch(2 - isTRUE('modLL' %in% names(net0)), net0$modLL$nodes, uni_mll(net0))
    if(is.logical(net1)){nodes <- net1; net1 <- NULL}
    if(!is.null(net1)){
      if(isTRUE(attr(net1, "mlGVAR"))){net1 <- net1$betweenNet}
      if(is.null(lrt)){lrt <- TRUE}
      omni1 <- switch(2 - isTRUE('modLL' %in% names(net1)), net1$modLL$omnibus, omni_mll(net1))
      uni1 <- switch(2 - isTRUE('modLL' %in% names(net1)), net1$modLL$nodes, uni_mll(net1))
      omni0 <- rbind(omni0, omni1)
      uni0 <- list(uni0, uni1)
      rownames(omni0) <- names(uni0) <- nn
    }
  } else {
    if(all(sapply(net0, function(z) isTRUE(attr(z, "mlGVAR"))))){
      net0 <- lapply(net0, '[[', "betweenNet")
    }
    stopifnot(all(sapply(net0, function(z) isTRUE(attr(z, "ggm")))))
    yLL <- all(sapply(net0, function(z) 'modLL' %in% names(z)))
    if(is.null(lrt)){lrt <- FALSE}
    if(is.logical(net1)){nodes <- net1}
    omni0 <- t(switch(2 - yLL, sapply(lapply(net0, '[[', 'modLL'), '[[', 'omnibus'), sapply(net0, omni_mll)))
    if(!is.null(orderBy)){
      orderBy <- switch(match.arg(tolower(as.character(
        orderBy)), c("true", "ll", "loglik", "aic", "bic", "df", "rank")),
        ll =, loglik = "LL", aic = "AIC", bic = "BIC", "df")
      net0 <- net0[order(omni0[, orderBy], decreasing = decreasing)]
      #net0 <- net0[order(sapply(net0, function(z){
      #  omni_mll(z)[orderBy]}), decreasing = decreasing)]
    }
    nn <- ifelse(!is.null(names(net0)), list(names(net0)),
                 list(paste0("net", 0:(length(net0) - 1))))[[1]]
    uni0 <- switch(2 - yLL, lapply(lapply(net0, '[[', 'modLL'), '[[', 'nodes'), lapply(net0, uni_mll))
    #omni0 <- do.call(rbind, lapply(net0, omni_mll))
    #uni0 <- lapply(net0, uni_mll)
    rownames(omni0) <- names(uni0) <- nn
  }

  # Collect output
  out <- ifelse(all, list(list(nodes = uni0, omnibus = omni0)),
                ifelse(nodes, list(uni0), list(omni0)))[[1]]

  # Do you want to conduct an LRT?
  if(ifelse(!is.null(net1), lrt, FALSE)){
    lrt0 <- mod_lrt(object = omni0, d = d, alpha = alpha, N = prod(dim(net1$dat$Y)))
    lrt1 <- mod_lrt(object = uni0, d = d, alpha = alpha)
    out <- list(nodes = uni0, LRT = lrt1, omnibus = lrt0)
    if(!all){out <- ifelse(nodes, list(lrt1), list(lrt0))[[1]]}
  }

  # RETURN
  return(out)
}

#' @rdname LogLikelihood
#' @export
SURll <- function(net0, net1 = NULL, nodes = FALSE, lrt = NULL, all = FALSE,
                  d = 4, alpha = .05, s = "res", orderBy = NULL,
                  decreasing = TRUE, sysfits = FALSE){
  sll <- function(fit, s = "res"){
    if("SURfit" %in% names(fit)){fit <- fit$SURfit}
    s <- match.arg(tolower(s), choices = c("res", "dfres", "sigma"))
    resid <- residuals(fit)
    residCov <- getCoefs(fit = fit, mat = s)
    residCovInv <- solve(residCov)
    resid <- as.matrix(resid)
    nEq <- ncol(resid)
    ll <- 0
    for(i in 1:nrow(resid)){
      ll <- ll - (nEq/2) * log(2 * pi) - .5 * log(det(residCov)) - .5 * resid[i, , drop = FALSE] %*% residCovInv %*% t(resid[i, , drop = FALSE])
    }
    df <- fit$rank + (nEq * (nEq + 1))/2
    out <- c(LL = as.numeric(ll), df = df)
    out
  }
  uni_sll <- function(fit){
    if("SURnet" %in% names(fit)){fit <- append(fit$SURnet, list(fitobj = fit$SURfit$eq))}
    getInd <- function(ind, x){
      ind <- c("deviance", "LL_model", "df.residual", "AIC", "BIC")[ind]
      X <- ifelse(ind == "df.residual", list(x$fitobj), list(x$mods))[[1]]
      return(unname(sapply(X, '[[', ind)))
    }
    out <- do.call(cbind.data.frame, lapply(1:5, getInd, x = fit))
    colnames(out) <- c("RSS", "LL", "df", "AIC", "BIC")
    rownames(out) <- names(fit$mods)
    out
  }
  omni_sll <- function(fit, s = "res"){
    k <- length(fit$SURnet$mods)
    n <- nrow(fit$SURnet$data$X) * k
    ll <- sll(fit = fit, s = s)
    aic <- (2 * ll[2]) - (2 * ll[1])
    bic <- (ll[2] * log(n)) - (2 * ll[1])
    out <- c(ll, AIC = unname(aic), BIC = unname(bic))
    out
  }
  sur_lrt <- function(object, d = 4, alpha = .05, N = NULL){
    if(is.list(object)){
      if(length(object) > 2){object <- object[1:2]}
      nn <- names(object)
      ll0 <- object[[1]]$LL; df0 <- object[[1]]$df
      ll1 <- object[[2]]$LL; df1 <- object[[2]]$df
      omnibus <- FALSE
    } else {
      if(nrow(object) > 2){object <- object[1:2, ]}
      nn <- rownames(object)
      ll0 <- object[1, 1]; df0 <- object[1, 2]
      ll1 <- object[2, 1]; df1 <- object[2, 2]
      omnibus <- TRUE
    }
    lldiff <- abs(ll0 - ll1) * 2
    dfdiff <- abs(df0 - df1)
    ps <- pchisq(q = lldiff, df = dfdiff, lower.tail = FALSE)
    decision <- c()
    for(i in seq_along(ps)){
      if(ps[i] <= alpha){
        decision[i] <- ifelse(ll0[i] > ll1[i], nn[1], nn[2])
      } else if(ps[i] == 1){
        decision[i] <- "- "
      } else if(ps[i] > alpha){
        if(omnibus){
          decision[i] <- ifelse(df0 < df1, nn[1], nn[2])
        } else {
          decision[i] <- ifelse(df0[i] > df1[i], nn[1], nn[2])
        }
      }
    }
    if(!omnibus){
      if(!is.null(d)){ps <- round(ps, d)}
      out <- data.frame(LL_diff2 = lldiff, Df_diff = dfdiff, pval = ps, decision = decision)
      rownames(out) <- rownames(object[[1]])
    } else {
      RMSEA <- function(X2, df, N){
        rmsea <- sqrt(max(c(((X2/N)/df) - (1/N), 0)))
        lower.l <- function(lambda){(pchisq(q = X2, df = df, ncp = lambda) - .95)}
        lambda.l <- tryCatch({uniroot(f = lower.l, lower = 0, upper = X2)$root}, error = function(e){0})
        rmsea.lower <- sqrt(lambda.l/(N * df))
        upper.l <- function(lambda){(pchisq(q = X2, df = df, ncp = lambda) - .05)}
        lambda.u <- tryCatch({uniroot(f = upper.l, lower = 0, upper = max(N, X2 * 4))$root}, error = function(e){1})
        rmsea.upper <- sqrt(lambda.u/(N * df))
        rmsea.pvalue <- 1 - pchisq(q = X2, df = df, ncp = (N * df * (.05^2)))
        return(c(lower = rmsea.lower, RMSEA = rmsea, upper = rmsea.upper, p.value = rmsea.pvalue))
      }
      rmsea <- RMSEA(lldiff, dfdiff, N)
      if(!is.null(d)){
        ps <- round(ps, d)
        lldiff <- round(lldiff, d)
        rmsea <- round(rmsea, d)
      }
      if(object[2, 2] < object[1, 2]){object <- object[order(object[, 2]), ]}
      out0 <- data.frame(LL_diff2 = c("", lldiff), Df_diff = c("", dfdiff),
                         pval = c("", ps), decision = c("", decision))
      out <- data.frame(object[, 1:2], out0, object[, 3:4])
      attr(out, "RMSEA") <- rmsea
    }
    return(out)
  }
  nn <- paste0("net", 0:1)
  if(length(net0) == 2 & ifelse(
    is.null(net1), TRUE, ifelse(is.logical(net1), net1, FALSE))){
    if(isTRUE(net1)){nodes <- net1}
    nn <- ifelse(!is.null(names(net0)), list(names(net0)), list(nn))[[1]]
    net1 <- net0[[2]]; net0 <- net0[[1]]
  }
  if(isTRUE(attr(net0, "mlGVAR"))){net0 <- net0$fixedNets}
  if("SURnet" %in% names(net0)){
    yLL <- isTRUE('SURll' %in% names(net0))
    omni0 <- switch(2 - yLL, net0$SURll$omnibus, omni_sll(fit = net0, s = s))
    uni0 <- switch(2 - yLL, net0$SURll$nodes, uni_sll(fit = net0))
    if(is.logical(net1)){nodes <- net1; net1 <- NULL}
    if(!is.null(net1)){
      if(isTRUE(attr(net1, "mlGVAR"))){net1 <- net1$fixedNets}
      if(is.null(lrt)){lrt <- TRUE}
      yLL <- isTRUE('SURll' %in% names(net1))
      omni1 <- switch(2 - yLL, net1$SURll$omnibus, omni_sll(fit = net1, s = s))
      uni1 <- switch(2 - yLL, net1$SURll$nodes, uni_sll(fit = net1))
      omni0 <- rbind(omni0, omni1)
      uni0 <- list(uni0, uni1)
      rownames(omni0) <- names(uni0) <- nn
    }
  } else {
    if(all(sapply(net0, function(z) isTRUE(attr(z, "mlGVAR"))))){
      net0 <- lapply(net0, '[[', "fixedNets")
    }
    stopifnot(all(sapply(net0, function(z) "SURnet" %in% names(z))))
    yLL <- all(sapply(net0, function(z) 'SURll' %in% names(z)))
    if(is.null(lrt)){lrt <- FALSE}
    if(is.logical(net1) & !sysfits){nodes <- net1}
    omni0 <- t(switch(2 - yLL, sapply(lapply(net0, '[[', 'SURll'), '[[', 'omnibus'), sapply(net0, omni_sll)))
    if(!is.null(orderBy)){
      orderBy <- switch(match.arg(tolower(as.character(
        orderBy)), c("true", "ll", "loglik", "aic", "bic", "df", "rank")),
        ll =, loglik = logLik, aic = AIC, bic = BIC, TRUE)
      if(is.function(orderBy)){
        orderBy <- c('LL', 'AIC', 'BIC')[which(sapply(
          list(logLik, AIC, BIC), function(z) identical(z, orderBy)))]
      }
      net0 <- net0[order(sapply(net0, function(z){
        if(isTRUE(orderBy)){
          return(sum(sapply(lapply(z$SURnet$mods, '[[', 'model'), nrow)))
        } else {
          return(omni0[, orderBy])
        }
      }), decreasing = decreasing)]
    }
    nn <- ifelse(!is.null(names(net0)), list(names(net0)),
                 list(paste0("net", 0:(length(net0) - 1))))[[1]]
    uni0 <- switch(2 - yLL, lapply(lapply(net0, '[[', 'SURll'), '[[', 'nodes'),
                   lapply(net0, uni_sll))
    rownames(omni0) <- names(uni0) <- nn
    sysgo <- all(sapply(net0, function(z) 'SURfit' %in% names(z)))
    if(sysfits & sysgo){
      net00 <- lapply(lapply(net0, '[[', "SURfit"), summary)
      ssr <- colSums(sapply(uni0, '[[', "RSS"))
      detrc <- sapply(net00, '[[', "detResidCov")
      olsr2 <- sapply(net00, '[[', "ols.r.squared")
      mcelroy <- sapply(net00, '[[', "mcelroy.r.squared")
      omni0 <- cbind(omni0, SSR = ssr, detSigma = detrc,
                     OLS.R2 = olsr2, McElroy.R2 = mcelroy)
      if(!is.null(d)){omni0 <- round(omni0, d)}
      return(omni0)
    }
  }
  out <- ifelse(all, list(list(nodes = uni0, omnibus = omni0)),
                ifelse(nodes, list(uni0), list(omni0)))[[1]]
  if(ifelse(!is.null(net1), lrt, FALSE)){
    lrt0 <- sur_lrt(object = omni0, d = d, alpha = alpha,
                    N = prod(dim(net1$SURnet$data$Y)))
    lrt1 <- sur_lrt(object = uni0, d = d, alpha = alpha)
    out <- list(nodes = uni0, LRT = lrt1, omnibus = lrt0)
    if(!all){out <- ifelse(nodes, list(lrt1), list(lrt0))[[1]]}
  }
  return(out)
}

#' @rdname LogLikelihood
#' @export
modTable <- function(net0, nodes = FALSE, orderBy = TRUE, d = 4, alpha = .05,
                     decreasing = TRUE, names = NULL, rmsea = FALSE){
  n <- length(net0)
  if(is.list(net0) & n <= 2){
    return(modLL(net0, nodes = nodes, orderBy = orderBy, d = d,
                 alpha = alpha, decreasing = decreasing))
  }
  if(!is.null(names)){names(net0) <- names}
  if(is.null(names(net0))){names(net0) <- paste0("fit", 1:n)}
  if(all(sapply(net0, function(z) isTRUE(attr(z, "mlGVAR"))))){
    net0 <- lapply(net0, "[[", "betweenNet")
  }
  stopifnot(all(sapply(net0, function(z) isTRUE(attr(z, "ggm")))))
  orderBy <- ifelse(is.null(orderBy), FALSE, ifelse(identical(tolower(orderBy), "lrt"), TRUE, orderBy))
  if(!is.logical(orderBy)){
    orderBy <- switch(match.arg(tolower(as.character(
      orderBy)), c("ll", "loglik", "aic", "bic", "df", "rank")),
      ll =, loglik = "LL", aic = "AIC", bic = "BIC", "df")
    net0 <- net0[order(sapply(net0, function(z){
      modLL(z)[orderBy]}), decreasing = decreasing)]
  }
  tt <- combn(n, 2)
  lls0 <- modLL(net0)
  lrts0 <- lapply(seq_len(ncol(tt)), function(i){
    modLL(net0[tt[, i]], d = d, alpha = alpha)})
  out1 <- t(sapply(lrts0, rownames))
  out2 <- t(sapply(lrts0, function(z) as.numeric(z[2, 3:5])))
  out3 <- sapply(lrts0, function(z) z[2, 6])
  out4 <- cbind.data.frame(out1, out2, out3)
  colnames(out4) <- c("net0", "net1", "Chisq", "Df", "pval", "decision")
  rmsea0 <- data.frame(out4[, 1:2], do.call(rbind, lapply(lrts0, attr, "RMSEA")))
  select <- table(out4$decision)
  if(any(names(select) == "- ")){select <- select[names(select) != "- "]}
  lls0 <- cbind(lls0, LRT = numeric(nrow(lls0)))
  lls0[match(names(select), rownames(lls0)), "LRT"] <- unname(select)
  if(isTRUE(orderBy)){
    lls0 <- lls0[order(lls0[, "LRT"], decreasing = TRUE), ]
  }
  out <- list(LRT = out4, omnibus = lls0, RMSEA = rmsea0)
  if(!rmsea){out$RMSEA <- NULL}
  if(nodes != FALSE){
    lls1 <- modLL(net0, nodes = TRUE)
    stopifnot(length(unique(lapply(lls1, rownames))) == 1)
    nodenames <- unique(lapply(lls1, rownames))[[1]]
    lls2 <- lapply(nodenames, function(fit){
      nodemod <- matrix(unlist(sapply(lls1, function(node){
        node[fit, -1]})), nrow = length(lls1), ncol = 4, byrow = TRUE)
      dimnames(nodemod) <- list(names(net0), c("LL", "df", "AIC", "BIC"))
      return(nodemod)
    })
    names(lls2) <- nodenames
    lrts1 <- lapply(seq_len(ncol(tt)), function(i){
      modLL(net0[tt[, i]], nodes = TRUE, d = d, alpha = alpha)})
    netLRTs <- lapply(colnames(lrts1[[1]]), function(nn){
      n1 <- data.frame(out4[, 1:2], "|", do.call(rbind, lapply(lrts1, "[[", nn)))
      colnames(n1)[-c(1:2)] <- c("|", nodenames)
      return(n1)
    })
    names(netLRTs) <- colnames(lrts1[[1]])
    deci <- netLRTs$decision[, -c(1:3)]
    deci2 <- do.call(cbind.data.frame, lapply(1:ncol(deci), function(z){
      z <- table(deci[, z])
      if(any(names(z) == "- ")){z <- z[names(z) != "- "]}
      zz <- setNames(numeric(length(net0)), names(net0))
      zz[match(names(z), names(zz))] <- unname(z)
      return(zz)
    }))
    colnames(deci2) <- colnames(deci)
    out <- list(nodes = lls2, LRT = netLRTs, counts = deci2)
    if(is.character(nodes)){
      out$decision <- netLRTs$decision
      out$LRT <- NULL
    }
  }
  attr(out, "alpha") <- alpha
  return(out)
}

#' @rdname LogLikelihood
#' @export
SURtable <- function(net0, nodes = FALSE, orderBy = TRUE, d = 4, alpha = .05,
                     decreasing = TRUE, names = NULL, rmsea = FALSE, s = "res"){
  n <- length(net0)
  if(is.list(net0) & n <= 2){
    return(SURll(net0, nodes = nodes, orderBy = orderBy, d = d,
                 alpha = alpha, decreasing = decreasing, s = s))
  }
  if(!is.null(names)){names(net0) <- names}
  if(is.null(names(net0))){names(net0) <- paste0("fit", 1:n)}
  if(all(sapply(net0, function(z) isTRUE(attr(z, "mlGVAR"))))){
    net0 <- lapply(net0, '[[', "fixedNets")
  }
  stopifnot(all(sapply(net0, function(z) "SURnet" %in% names(z))))
  if(!is.null(orderBy)){
    oms <- SURll(net0)
    orderBy <- switch(match.arg(tolower(as.character(
      orderBy)), c("true", "ll", "loglik", "aic", "bic", "df", "rank")),
      ll =, loglik = logLik, aic = AIC, bic = BIC, TRUE)
    if(is.function(orderBy)){
      orderBy <- c('LL', 'AIC', 'BIC')[which(sapply(
        list(logLik, AIC, BIC), function(z) identical(z, orderBy)))]
    }
    net0 <- net0[order(sapply(net0, function(z){
      if(isTRUE(orderBy)){
        return(sum(sapply(lapply(z$SURnet$mods, '[[', 'model'), nrow)))
      } else {
        return(oms[, orderBy])
      }
    }), decreasing = decreasing)]
  }
  tt <- combn(n, 2)
  lls0 <- SURll(net0, s = s)
  lrts0 <- lapply(seq_len(ncol(tt)), function(i){
    SURll(net0[tt[, i]], d = d, alpha = alpha, s = s)})
  out1 <- t(sapply(lrts0, rownames))
  out2 <- t(sapply(lrts0, function(z) as.numeric(z[2, 3:5])))
  out3 <- sapply(lrts0, function(z) z[2, 6])
  out4 <- cbind.data.frame(out1, out2, out3)
  colnames(out4) <- c("net0", "net1", "Chisq", "Df", "pval", "decision")
  rmsea0 <- data.frame(out4[, 1:2], do.call(rbind, lapply(lrts0, attr, "RMSEA")))
  select <- table(out4$decision)
  if(any(names(select) == "- ")){select <- select[names(select) != "- "]}
  lls0 <- cbind(lls0, LRT = numeric(nrow(lls0)))
  lls0[match(names(select), rownames(lls0)), "LRT"] <- unname(select)
  out <- list(LRT = out4, omnibus = lls0, RMSEA = rmsea0)
  if(!rmsea){out$RMSEA <- NULL}
  if(nodes != FALSE){
    lls1 <- SURll(net0, nodes = TRUE, s = s)
    stopifnot(length(unique(lapply(lls1, rownames))) == 1)
    nodenames <- unique(lapply(lls1, rownames))[[1]]
    lls2 <- lapply(nodenames, function(fit){
      nodemod <- matrix(unlist(sapply(lls1, function(node){
        node[fit, -1]})), nrow = length(lls1), ncol = 4, byrow = TRUE)
      dimnames(nodemod) <- list(names(net0), c("LL", "df", "AIC", "BIC"))
      return(nodemod)
    })
    names(lls2) <- nodenames
    lrts1 <- lapply(seq_len(ncol(tt)), function(i){
      SURll(net0[tt[, i]], nodes = TRUE, d = d, alpha = alpha, s = s)})
    netLRTs <- lapply(colnames(lrts1[[1]]), function(nn){
      n1 <- data.frame(out4[, 1:2], "|", do.call(rbind, lapply(lrts1, '[[', nn)))
      colnames(n1)[-c(1:2)] <- c("|", nodenames)
      return(n1)
    })
    names(netLRTs) <- colnames(lrts1[[1]])
    deci <- netLRTs$decision[, -c(1:3)]
    deci2 <- do.call(cbind.data.frame, lapply(1:ncol(deci), function(z){
      z <- table(deci[, z])
      if(any(names(z) == "- ")){z <- z[names(z) != "- "]}
      zz <- setNames(numeric(length(net0)), names(net0))
      zz[match(names(z), names(zz))] <- unname(z)
      return(zz)
    }))
    colnames(deci2) <- colnames(deci)
    out <- list(nodes = lls2, LRT = netLRTs, counts = deci2)
    if(is.character(nodes)){
      out$decision <- netLRTs$decision
      out$LRT <- NULL
    }
  }
  attr(out, "alpha") <- alpha
  return(out)
}
