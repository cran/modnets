#' Get adjacency matrices from fit objects
#'
#' \code{\link{net}} returns the adjacency matrix for any network model fit
#' using functions from the \code{modnets} package. \code{\link{netInts}}
#' returns a matrix of interaction terms associated with a moderated network.
#'
#' For GGMs when a non-symmetric matrix is requested, columns will represent
#' outcomes and rows will represent predictors. For temporal networks, columns
#' represent predictors and rows represent outcomes.
#'
#' Can also be used with output from the \code{\link{resample}} and
#' \code{\link{bootNet}} functions.
#'
#' @param fit A fitted network model. Can be the output from
#'   \code{\link{fitNetwork}}, \code{\link{mlGVAR}}, \code{\link{lmerVAR}},
#'   \code{\link{bootNet}}, \code{\link{resample}}, \code{\link{simNet}}, or
#'   \code{\link{mlGVARsim}}.
#' @param n When multiple networks exist for a single object, this allows the
#'   user to indicate which adjacency matrix to return. For a GGM, all values of
#'   this argument return the same adjacency matrix. For a SUR network,
#'   \code{"beta"} and \code{"temporal"} return the coefficients associated with
#'   the temporal network, while \code{"pdc"} returns the Partial Directed
#'   Correlations, or the standardized temporal network.
#'   \code{"contemporaneous"} and \code{"pcc"} return the standardized
#'   contemporaneous network (Partial Contemporaneous Correlations).
#'   \code{"kappa"} returns the unstandardized residual covariance matrix. All
#'   of these terms apply for multilevel networks, but \code{"between"} can also
#'   return the between-subjects network. If a numeric or logical value is
#'   supplied, however, this argument will function as the \code{threshold}
#'   argument. A numeric value will set a threshold at the supplied value, while
#'   \code{TRUE} will set a threshold of .05.
#' @param threshold A numeric or logical value to set a p-value threshold.
#'   \code{TRUE} will automatically set the threshold at .05.
#' @param rule Only applies to GGMs (including between-subjects networks) when a
#'   threshold is supplied. The \code{"AND"} rule will only preserve edges when
#'   both corresponding coefficients have p-values below the threshold, while
#'   the \code{"OR"} rule will preserve an edge so long as one of the two
#'   coefficients have a p-value below the supplied threshold.
#' @param binary Logical. If \code{TRUE} then the weighted adjacency matrix will
#'   be converted into an unweighted adjacency matrix.
#' @param nodewise Logical, only applies to GGMs (including between-subjects
#'   networks). If \code{TRUE} then the adjacency matrix will retain all
#'   coefficients in their original form. In this case, values in rows represent
#'   the coefficients predicting the columns.
#' @param d Numeric. Only used for output of \code{\link{mlGVARsim}}, or
#'   \code{\link{simNet}} when \code{lags = 1}. Sets the number of decimal
#'   places to round the output to.
#' @param r Numeric. Chooses which rows/columns to remove from the output, if
#'   desired.
#' @param avg Logical. For \code{\link{netInts}}, determines whether to take the
#'   average two corresponding interaction terms.
#' @param empty Logical. Determines the output of \code{\link{netInts}} when
#'   \code{fit} is not a moderated network. If \code{TRUE} then an empty list
#'   will be returned. If \code{FALSE} then a matrix of zeros will be returned.
#' @param mselect Only used for \code{\link{netInts}} when there is more than
#'   one exogenous moderator. Allows the user to indicate which moderator should
#'   be used to construct the interaction matrix.
#'
#' @return An adjacency matrix representing a network or a matrix of interaction
#'   terms.
#' @export
#'
#' @seealso \code{\link{fitNetwork}, \link{mlGVAR}, \link{lmerVAR},
#'   \link{bootNet}, \link{resample}, \link{simNet}, \link{mlGVARsim}}
#'
#' @examples
#' x <- fitNetwork(ggmDat, 'M')
#'
#' net(x, threshold = .05)
#' netInts(x, threshold = TRUE)
#'
#' \donttest{
#' y <- mlGVAR(mlgvarDat, 'M')
#'
#' net(y, n = 'beta')
#' net(y, n = 'pcc')
#' net(y, n = 'between')
#'
#' netInts(y)
#' }
net <- function(fit, n = "beta", threshold = FALSE, rule = "OR",
                binary = FALSE, nodewise = FALSE, d = 14, r = NULL){
  if(inherits(fit, c('splitNets', 'try-error'))){return(NULL)}
  if(is(fit, 'bootNet') | (isTRUE(attr(fit, 'resample')) & 'fit0' %in% names(fit))){fit <- fit$fit0}
  if(is(fit, 'ggmSim')){return(fit[[grep('^trueNet$|^b1$|^adjMat$', names(fit))]])}
  if(is(fit, "matrix")){return(fit)}
  if(is(fit, "mgm")){
    out <- fit$pairwise$wadj * replace(fit$pairwise$signs, is.na(fit$pairwise$signs), 0)
    if(!is.null(r)){out <- out[-r, -r]}
    return(out)
  }
  if(isTRUE(n) | is.numeric(n)){
    if(tolower(threshold) %in% c('and', 'or')){rule <- threshold}
    if(is.numeric(n)){
      threshold <- n
      n <- "beta"
    } else {
      n <- threshold <- "beta"
    }
  }
  n <- match.arg(tolower(n), c(
    "beta", "contemporaneous", "between", "pdc", "pcc", "kappa", "temporal"))
  n1 <- switch(n, between = "", beta = , pdc = , temporal = "$temporal", "$contemporaneous")
  n2 <- switch(n, kappa = "kappa", pdc = "PDC$adjMat", "adjMat")
  n3 <- paste0("fit", n1, "$", n2)
  if(isTRUE(attr(fit, "mlGVAR"))){fit <- if(n == "between"){fit$betweenNet} else {fit$fixedNets}}
  if("SURnet" %in% names(fit)){fit <- fit$SURnet}
  if(isTRUE(attr(fit, "ggm"))){n3 <- ifelse(nodewise, "fit$nodewise$adjNW", "fit$adjMat")}
  if("fixedPCC" %in% names(fit)){
    threshold <- FALSE
    n <- switch(n, contemporaneous = "pcc", temporal = "beta", n)
    if(n != "between" & "fixedResults" %in% names(fit)){fit <- fit$fixedResults}
    out <- fit[[grep(n, tolower(names(fit)))[1]]]
    if(n == "beta" & ncol(out) != nrow(out)){out <- out[, -1]}
    if(n == "between"){diag(out) <- 0}
    if(n == "pdc"){out <- t(out)}
    if(!is.null(d)){out <- round(out, d)}
  } else {
    if(isTRUE(attr(fit, "lmerVAR")) & n == "between"){n3 <- paste0("fit$between$adjMat")}
    out <- eval(parse(text = n3))
  }
  if(threshold != FALSE){
    if(!is.numeric(threshold)){threshold <- .05}
    rule <- match.arg(tolower(rule), c("or", "and"))
    atts <- names(attributes(fit))
    if("SURnet" %in% atts){
      stopifnot(!n == "between")
      n4 <- ifelse(n %in% c("beta", "temporal", "pdc"), "temporal", "contemporaneous")
      n5 <- ifelse(n4 == "temporal", "fit$temporal$coefs", "fit$contemporaneous")
      pvals <- eval(parse(text = n5))$pvals
      if(ncol(pvals) != nrow(pvals)){
        if(any(grepl("lag", colnames(pvals)))){
          colnames(pvals) <- gsub("[.]lag1[.]$", "", colnames(pvals))
        }
        pvals <- pvals[, intersect(colnames(pvals), gsub("[.]y$", "", rownames(pvals)))]
      }
    } else {
      n4 <- ifelse(n %in% c("beta", "temporal", "pdc"), "Beta",
                   ifelse(n == "between", "gammaOmega", "gammaTheta"))
      pvals <- if("ggm" %in% atts){fit$nodewise$pvalsNW} else {fit$coefs[[n4]]$Pvals}
      if("ggm" %in% atts){n4 <- ifelse(nodewise & rule == "or", "Beta", "")}
    }
    if(any(is.na(pvals))){pvals[is.na(pvals)] <- 1}
    if(n4 %in% c("temporal", "Beta")){
      out <- out * ifelse(pvals <= threshold, 1, 0)
    } else if(rule == "or"){
      out <- out * ifelse(pvals <= threshold | t(pvals) <= threshold, 1, 0)
    } else if(rule == "and"){
      out <- out * ifelse(pvals <= threshold & t(pvals) <= threshold, 1, 0)
    }
  }
  if(binary){out <- abs(sign(out))}
  if(!is.null(r)){out <- out[-r, -r]}
  return(out)
}

#' @rdname net
#' @export
netInts <- function(fit, n = 'temporal', threshold = FALSE, avg = FALSE,
                    rule = 'none', r = NULL, empty = TRUE, mselect = NULL){
  rules <- c('none', 'or', 'and')
  eout <- function(fit, empty = TRUE){
    n <- tryCatch({ncol(net(fit))}, error = function(e){TRUE})
    cond <- isTRUE(empty | isTRUE(n) | is.null(n) | is.na(n))
    return(switch(2 - cond, list(), diag(0, n)))
  }
  if(is(fit, 'splitNets')){return(fit$ints)}
  if(is(fit, 'mgm')){
    m <- ifelse(!is.null(r), r, fit$call$moderators[length(fit$call$moderators)])
    mgmInt <- function(x, m){
      ni <- out <- x$interactions$indicator
      if(length(ni) == 2){
        ni <- t(apply(ni[[2]], 1, function(z) setdiff(z, m)))
        vals <- unlist(x$interactions$weightsAgg[[2]]) * x$interactions$signs[[2]]
        out <- matrix(0, ncol(net(x)) - 1, ncol(net(x)) - 1)
        for(i in seq_along(vals)){out[ni[i, 1], ni[i, 2]] <- vals[i]}
        out <- out + t(out)
      }
      return(out)
    }
    out <- mgmInt(fit, m)
    if(!is.null(r) & !is(out, 'list')){out <- out[-r, -r]}
    if(is(out, 'list')){return(eout(fit, empty))} else {return(out)}
  }
  if(is(fit, 'ggmSim')){if('b2' %in% names(fit)){return(fit$b2)} else {return(eout(fit, empty))}}
  if(is(fit, 'bootNet') | (isTRUE(attr(fit, 'resample')) & 'fit0' %in% names(fit))){fit <- fit$fit0}
  if(isTRUE(n) & isTRUE(threshold)){avg <- TRUE}
  if(isTRUE(n) | is.numeric(n)){
    if(tolower(as.character(threshold)) %in% rules){rule <- threshold}
    if(is.numeric(n)){
      threshold <- n
      n <- "temporal"
    } else {
      n <- threshold <- "temporal"
    }
  }
  n <- match.arg(tolower(n), c("temporal", "between"))
  atts <- names(attributes(fit))
  if("mlGVAR" %in% atts){fit <- switch(2 - isTRUE(n == "temporal"), fit$fixedNets, fit$betweenNet)}
  if("SURnet" %in% names(fit)){fit <- fit$SURnet}
  if("lmerVAR" %in% atts){
    stopifnot("ints" %in% names(fit$coefs))
    out <- fit$coefs$ints$coefs
    pvals <- fit$coefs$ints$Pvals
  } else if(any(c('mlGVARsim', 'GVARsim') %in% atts)){ # NEW
  #} else if(any(c('mlVARsim', 'simMLgvar') %in% atts)){
    #stopifnot("mm" %in% names(fit))
    #out <- fit$mm$mb2
    stopifnot('interactions' %in% names(fit)) # NEW
    out <- fit$interactions$mb2 # NEW
  } else if(!'interactions' %in% names(fit)){
    out <- tryCatch({mmat(fit = fit, m = mselect)}, error = function(e){TRUE})
    if(isTRUE(out)){return(eout(fit, empty))}
    pvals <- out$pvals
    out <- out$betas
  } else {
    out <- fit$interactions[[1]]
    if(isTRUE(attr(fit, "ggm"))){
      pvals <- out[[2]][[2]]
      out <- out[[2]][[1]]
    } else {
      rule <- 'none'
      pvals <- fit$interactions$pvals
      if(ncol(out) != nrow(out) & !is.null(mselect)){
        if(isTRUE(mselect)){mselect <- fit$call$moderators[1]}
        m1 <- paste0(':', mselect); m2 <- paste0(mselect, ':')
        mm <- paste0(c(m1, m2), collapse = '|')
        my <- paste0(gsub('[.]y$', '', rownames(out)), collapse = '|')
        mb1 <- out[, grep(mm, colnames(out))]
        mb1 <- mb1[, grep(my, colnames(mb1))]
        mp1 <- pvals[, grep(mm, colnames(pvals))]
        mp1 <- mp1[, grep(my, colnames(mp1))]
        out <- mb1
        pvals <- mp1
      }
    }
  }
  if(threshold != FALSE & !"mlGVARsim" %in% atts){ # NEW
  #if(threshold != FALSE & !"mlVARsim" %in% atts){
    rule <- match.arg(tolower(rule), rules)
    if(isTRUE(attr(fit, 'ggm')) & rule == 'none'){rule <- 'or'}
    if(!is.numeric(threshold)){threshold <- .05}
    if(rule == 'none' | isTRUE(attr(fit, 'SURnet'))){
      out <- out * ifelse(pvals <= threshold, 1, 0)
    } else if(rule == 'or'){
      out <- out * ifelse(pvals <= threshold | t(pvals) <= threshold, 1, 0)
    } else if(rule == 'and'){
      out <- out * ifelse(pvals <= threshold & t(pvals) <= threshold, 1, 0)
    }
  }
  if(avg){out <- (t(out) + out)/2}
  if(!is.null(r)){out <- out[-r, -r]}
  if(is(out, 'list')){return(eout(fit, empty))}
  return(out)
}
