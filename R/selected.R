#' Shows which variables were selected for each node of a network
#'
#' Provides a quick representation showing which variables were selected as
#' predictors of each node in a network, both for unmoderated and moderated
#' networks. Especially useful as a way to see which variables were selected in
#' a variable selection procedure, such as through the \code{\link{varSelect}}
#' and \code{\link{resample}} functions.
#'
#' The \code{threshold} argument allows the user to set a threshold for
#' p-values, such that the output only reflects the predictors that are
#' significant at that threshold. This argument can be utilized whether or not a
#' variable selection procedure has been employed.
#'
#' @param object Output from either \code{\link{fitNetwork}} or
#'   \code{\link{mlGVAR}}
#' @param threshold Can be a numeric value between \code{0} and \code{1}, or
#'   defaults to \code{.05} when set to \code{TRUE}
#' @param mod Only relevant to models fit with \code{\link{mlGVAR}}
#'
#' @return A table where the columns represent nodes, and the rows show which
#'   variables were selected in predicting each node. For moderated networks,
#'   the output is a list that separates main effects (\code{mods}) from
#'   interaction effects (\code{ints}).
#' @export
#'
#' @seealso \code{\link{fitNetwork}, \link{mlGVAR}}
#'
#' @examples
#' \donttest{
#' fit1 <- fitNetwork(ggmDat)
#' selected(fit1)
#'
#' fit2 <- mlGVAR(mlgvarDat, m = 'M', verbose = FALSE)
#' selected(fit2, threshold = TRUE, mod = 'temporal') # Can also set to 'between'
#'
#' fit3 <- fitNetwork(gvarDat, moderators = 'M', type = 'varSelect', lags = 1)
#' selected(fit3)
#' }
selected <- function(object, threshold = FALSE, mod = c('temporal', 'between')){
  ints <- NULL
  if(threshold != FALSE & !is.numeric(threshold)){threshold <- .05}
  if(isTRUE(attr(object, "mlGVAR"))){
    mod <- match.arg(tolower(mod), c('temporal', 'between'))
    object <- object[[ifelse(mod == 'temporal', 2, 3)]]
  }
  if(isTRUE(attr(object, "SURnet"))){
    if('SURnet' %in% names(object)){object <- object$SURnet}
    mods0 <- object$temporal$coefs[[ifelse(threshold == FALSE, 1, 2)]]
    mods <- lapply(seq_len(nrow(mods0)), function(z){
      if(threshold == FALSE){
        z <- colnames(mods0)[mods0[z, ] != 0]
      } else {
        z <- colnames(mods0)[mods0[z, ] <= threshold]
      }
      z <- replace(z, length(z) == 0, "")
      z <- replace(z, is.null(z), "")
      return(z[z != "(Intercept)"])
    })
    names(mods) <- rownames(mods0)
  } else {
    if(threshold != FALSE & "fitobj" %in% names(object)){
      mods <- lapply(object$fitobj, function(z){
        z <- summary(z)$coefficients[, 4]
        z <- ifelse(z <= threshold, 1, 0)
        z <- names(z)[which(z == 1)]
        z <- replace(z, length(z) == 0, "")
        z <- replace(z, is.null(z), "")
        return(z[z != "(Intercept)"])
      })
    } else {
      mods <- lapply(object$mods, function(z){
        z <- rownames(z$model)[-1]
        return(replace(z, length(z) == 0, ""))
      })
    }
  }
  if(any(grepl(":", mods))){
    ints <- lapply(mods, function(z){
      z <- z[grepl(":", z)]
      return(replace(z, length(z) == 0, ""))
    })
    ix <- max(sapply(ints, length))
    ints <- do.call(cbind.data.frame, lapply(ints, function(z){
      if(length(z) < ix){z <- c(z, rep("", ix - length(z)))}
      return(z)
    }))
    mods <- lapply(mods, function(z){
      z <- z[!grepl(":", z)]
      return(replace(z, length(z) == 0, ""))
    })
  }
  mx <- max(sapply(mods, length))
  mods <- do.call(cbind.data.frame, lapply(mods, function(z){
    if(length(z) < mx){z <- c(z, rep("", mx - length(z)))}
    return(z)
  }))
  out <- mods
  if(!is.null(ints)){out <- list(mods = mods, ints = ints)}
  return(out)
}
