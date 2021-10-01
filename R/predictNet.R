#' Calculate prediction error from network models
#'
#' See the prediction error based on different statistics for either GGMs or
#' SURs. Also can compare and find the change values (such as R-squared change)
#' between two networks of the same size (i.e., with the same nodes).
#'
#' @param object Output from \code{\link{fitNetwork}} or \code{\link{mlGVAR}}.
#'   If using output from \code{\link{mlGVAR}}, then one of the two networks
#'   must be provided (i.e., either \code{fixedNets} or \code{betweenNet}).
#' @param data The dataset used to fit the network model, or another network of
#'   the same type and size to be compared with the network specified in the
#'   first argument. If the prediction error for only one network is desired,
#'   and the dataset is included as an element of the relevant object, then this
#'   can be left as \code{NULL}.
#' @param all if \code{TRUE} then returns a list containing the observed
#'   outcomes used to fit the models, their predicted values, and the prediction
#'   error for each outcome.
#' @param scale Logical; determines whether or not to standardize the data
#'   before computing prediction error. This argument will be removed.
#'
#' @return A table showing different measures of prediction error associated
#'   with each node of the network. Or, if two networks are provided, a table
#'   that shows the difference in prediction error for each node across the two
#'   networks. Specifically, this is computed by taking the statistics for
#'   \code{data} and subtracting them from those for \code{object}.
#'
#'   If \code{all = TRUE}, then the following output is returned: \describe{
#'   \item{Y}{The observed values of the outcome variables based on the data
#'   provided.} \item{preds}{The predicted values of the outcomes based on the
#'   models provided.} \item{errors}{Table containing prediction error
#'   statistics for each node.} }
#'
#' @export
#'
#' @seealso \code{\link{fitNetwork}}
#'
#' @examples
#' fit1 <- fitNetwork(ggmDat, covariates = 'M')
#' fit2 <- fitNetwork(ggmDat, moderators = 'M')
#'
#' predictNet(fit1)
#' predictNet(fit1, all = TRUE)
#'
#' predictNet(fit2, fit1) # Find the differences in prediction error across the two models
predictNet <- function(object, data = NULL, all = FALSE, scale = FALSE){
  # Use predictDelta if there are two networks
  if(is(data, 'ggm') | is(data, 'SURnet')){
    return(do.call(predictDelta, list(mod1 = object, mod0 = data, scale = scale)))
  }

  # Function to create design matrix for all models
  interactionMatrix <- function(data, y, type, moderators = NULL,
                                lags = NULL, threeWay = TRUE){
    p <- ncol(data)
    mainMod <- paste(colnames(data)[-y], collapse = " + ")
    if(!is.null(lags)){
      if(!is.null(moderators)){
        p <- ((p - 1)/length(lags)) + 1
        moderators <- moderators + 1
        ints <- t(combn((1:p)[-y], 2))
        ints_m <- c()
        for(i in 1:nrow(ints)){ints_m[i] <- any(ints[i,] %in% moderators)}
        ints <- ints[ints_m,]
        if(length(lags) == 1){
          intMods <- paste0(colnames(data)[ints[,1]], "*", colnames(data)[ints[,2]], collapse = "+")
        } else {
          lagInts <- list(); intMods <- list()
          for(j in 1:length(lags)){
            lagInts[[j]] <- data.frame(y = data[,"y"], data[,grep(paste0("lag", lags[j]), colnames(data))])
            intMods[[j]] <- paste0(colnames(lagInts[[j]])[ints[,1]], "*", colnames(lagInts[[j]])[ints[,2]], collapse = "+")
          }
          for(i in 1:(length(lags) - 1)){intMods[[i]] <- paste0(intMods[[i]], " + ")}
          intMods <- do.call(paste0, intMods)
        }
        form <- as.formula(paste0(colnames(data)[y], " ~ ", mainMod, " + ", intMods))
      } else {
        form <- as.formula(paste0(colnames(data)[y], " ~ ", mainMod))
      }
    } else {
      if(!is.null(moderators) & (if((y %in% moderators) & threeWay == FALSE){length(moderators) > 1} else {TRUE})){
        if((y %in% moderators) & threeWay == TRUE){
          ints <- t(combn((1:p)[-y], 2))
          intMods <- paste0("V", ints[,1], ".*V", ints[,2], ".", collapse = " + ")
        } else {
          if(y %in% moderators){moderators <- moderators[-which(moderators == y)]}
          nm <- length(moderators)
          intMods <- list()
          for(i in 1:nm){
            intMods[[i]] <- paste0(colnames(data)[-c(y, moderators[i])], "*V", moderators[i], ".", collapse = " + ")
          }
          if(nm > 1){for(i in 1:(nm - 1)){intMods[[i]] <- paste0(intMods[[i]], " + ")}}
          intMods <- do.call(paste0, intMods)
        }
        form <- as.formula(paste0(colnames(data)[y], " ~ ", mainMod, " + ", intMods))
      } else {
        form <- as.formula(paste0(colnames(data)[y], " ~ ", mainMod))
      }
    }
    X <- model.matrix(form, data)[,-1]
    X
  }

  # SURnet vs. ggm
  if("SURnet" %in% c(names(object), names(attributes(object)))){
    if("SURnet" %in% names(object)){object <- object$SURnet}
    data <- object$data
    mods <- object$mods
    type <- object$call$type
    moderators <- object$call$moderators
    if(is.character(moderators)){moderators <- which(colnames(data$X) %in% moderators)}
    lags <- object$call$lags
    predObj <- Y <- list()
    for(i in 1:length(mods)){
      predObj[[i]] <- cbind(1, data$X) %*% mods[[i]]$model
      Y[[i]] <- as.numeric(data$Y[, i])
    }
    Y <- data.frame(do.call(cbind, Y))
    yhat <- data.frame(as.matrix(do.call(cbind, predObj)))
    colnames(Y) <- colnames(data$Y)
    colnames(yhat) <- paste0(colnames(data$Y), "hat")
    data <- data$Y
    probObj <- vector("list", length = ncol(data))
  } else {
    if(is.null(data)){
      x <- try(data <- object$data)
      if(class(x) == "try-error"){stop("Must supply dataset")}
    }
    mods <- object$mods
    type <- object$call$type
    moderators <- object$call$moderators
    if(is.character(moderators)){
      if("mods0" %in% names(object)){
        moderators <- which(colnames(object$mods0$dat) %in% moderators)
      } else {
        moderators <- which(colnames(data) %in% moderators)
      }
    }
    if("lags" %in% names(object$call)){lags <- object$call$lags} else {lags <- NULL}
    if(!is.null(moderators)){if(length(moderators) == 1){if(moderators == 0){moderators <- NULL}}}
    if(is.null(colnames(data))){colnames(data) <- paste0("V", 1:ncol(data))}
    predObj <- list(); Y <- list(); probObj <- vector("list", length = ncol(data))
    if("ggm" %in% names(attributes(object))){
      mods <- lapply(mods, function(z) list(model = z$model))
      type <- unname(sapply(object$fitobj, attr, "family"))
      if("gaussian" %in% type){type[type == "gaussian"] <- "g"}
      if("binomial" %in% type){type[type == "binomial"] <- "c"}
      for(i in 1:length(mods)){
        if(type[i] == "c"){
          dat <- setup(data = object$mods0$dat, type = type, y = i, lags = lags, scale = scale)
          if(!is.null(lags)){y <- 1} else {y <- i}
          X <- interactionMatrix(data = dat, y = y, type = type, moderators = moderators, lags = lags)
          n <- nrow(X); p <- ncol(X)
          coefs <- list(mods[[i]]$model * -1, mods[[i]]$model)
          n_cats <- length(coefs)
          potentials <- matrix(NA, nrow(data), n_cats)
          for(cat in 1:n_cats){potentials[, cat] <- exp(cbind(1, X) %*% coefs[[cat]])[, 1]}
          probs <- potentials[, 1:n_cats]/rowSums(potentials[, 1:n_cats])
          probObj[[i]] <- probs
          predObj[[i]] <- sort(unique(dat[, y]))[apply(probs, 1, which.max)]
          Y[[i]] <- as.numeric(dat[, y])
        } else {
          predObj[[i]] <- predict(object = object$fitobj[[i]])
          Y[[i]] <- as.numeric(data[, i])
        }
      }
    } else {
      for(i in 1:length(mods)){
        dat <- setup(data = data, type = type, y = i, lags = lags, scale = scale)
        if(!is.null(lags)){y <- 1} else {y <- i}
        X <- interactionMatrix(data = dat, y = y, type = type, moderators = moderators, lags = lags)
        n <- nrow(X); p <- ncol(X)
        if(type[i] == "c"){
          coefs <- mods[[i]]$model
          n_cats <- length(coefs)
          potentials <- matrix(NA, n, n_cats)
          for(cat in 1:n_cats){potentials[,cat] <- exp(cbind(1, X) %*% coefs[[cat]])[,1]}
          probs <- potentials[,1:n_cats]/rowSums(potentials[,1:n_cats])
          probObj[[i]] <- probs
          predObj[[i]] <- sort(unique(dat[,y]))[apply(probs, 1, which.max)]
        } else {
          predObj[[i]] <- cbind(1, X) %*% mods[[i]]$model
        }
        Y[[i]] <- as.numeric(dat[,y])
      }
    }
    Y <- data.frame(do.call(cbind, Y))
    yhat <- data.frame(as.matrix(do.call(cbind, predObj)))
    colnames(Y) <- paste0(colnames(data), ".y")
    colnames(yhat) <- paste0(colnames(data), ".yhat")
    names(probObj) <- paste0(colnames(data), ".probs")
  }

  # Calculate the type of error/predictive statistic
  calcError <- function(y, yhat, type, mod){
    if(type == "g"){
      R2 <- 1 - sum((y - yhat)^2)/sum((y - mean(y))^2)
      adjR2 <- 1 - (1 - R2) * ((length(y) - 1)/(length(y) - sum(mod != 0)))
      MSE <- sum((y - yhat)^2)/(length(y) - sum(mod != 0))
      RMSE <- sqrt(MSE)
      predError <- list(R2 = R2, adjR2 = adjR2, MSE = MSE, RMSE = RMSE)
    }
    if(type == "c"){
      CC <- sum((y == yhat)/length(y))
      nCC <- (CC - max(table(y)/sum(table(y))))/(1 - max(table(y)/sum(table(y))))
      CCmarg <- (nCC - CC)/(nCC - 1)
      predError <- list(CC = CC, nCC = nCC, CCmarg = CCmarg)
    }
    lapply(predError, round, 3)
  }

  # Compute errors
  errors <- lapply(1:ncol(data), function(z) calcError(Y[,z], yhat[,z], type[z], mods[[z]]$model))

  # Make everything into a pretty table
  if(all(c("g", "c") %in% type)){
    errors <- lapply(errors, function(z) do.call(cbind, z))
    for(i in 1:length(errors)){
      if("R2" %in% colnames(errors[[i]])){
        errors[[i]] <- c(errors[[i]], rep(NA, 3))
      } else {
        errors[[i]] <- c(rep(NA, 4), errors[[i]])
      }
    }
    errors <- data.frame(Variable = colnames(data), do.call(rbind, errors))
    colnames(errors)[2:8] <- c("R2", "adjR2", "MSE", "RMSE", "CC", "nCC", "CCmarg")
  } else {
    errors <- data.frame(Variable = colnames(data), do.call(rbind, errors))
  }

  # Decide what to return
  if(all == TRUE){
    output <- list(Y = Y, preds = yhat, probs = probObj, errors = errors)
    if(all(type == "g")){output$probs <- NULL}
  } else {
    output <- data.frame(Variable = errors[, 1], apply(errors[, -1], 2, unlist))
  }

  # RETURN
  return(output)
}

##### predictDelta: Calculate change in prediction for across nested models
predictDelta <- function(mod1, mod0, scale = FALSE){
  # SURnet vs. ggm
  if("SURnet" %in% names(mod0)){
    stopifnot(ncol(mod0$SURnet$data$Y) == ncol(mod1$SURnet$data$Y))
    type <- rep("g", ncol(mod0$SURnet$data$Y))
    r0 <- attr(mod0, 'rank')
    r1 <- attr(mod1, 'rank')
    if(r1 < r0){mod00 <- mod1; mod1 <- mod0; mod0 <- mod00}
  } else {
    type <- mod1$call$type
    if("ggm" %in% names(attributes(mod1))){
      if(is.list(type)){
        type <- unname(sapply(mod1$fitobj, attr, "family"))
        if("gaussian" %in% type){type[type == "gaussian"] <- "g"}
        if("binomial" %in% type){type[type == "binomial"] <- "c"}
      } else if(length(type) == 1){
        type <- rep(substr(type, 1, 1), ncol(data))
      }
    } else {
      stopifnot(mod1$call$type == mod0$call$type)
    }
  }

  # Compute errors
  errors <- list(predictNet(mod0, all = FALSE, scale = scale), predictNet(mod1, all = FALSE, scale = scale))

  # See how many nodes are in each network, and whether they can be compared
  nodes <- unique(sapply(errors, nrow))
  if(length(nodes) > 1){stop('Cannot compare networks of different sizes')}

  # Collect variable types
  p <- ifelse(all(type == "c"), 4, ifelse(all(type == "g"), 5, 8))

  # Make matrix of comparisons
  delta <- data.frame(matrix(NA, nrow = length(nodes), ncol = p))
  for(i in seq_len(nodes)){
    if(type[i] == "g"){
      delta[i,2] <- unlist(errors[[2]][i, "R2"]) - unlist(errors[[1]][i, "R2"])
      delta[i,3] <- unlist(errors[[2]][i, "adjR2"]) - unlist(errors[[1]][i, "adjR2"])
      delta[i,4] <- unlist(errors[[2]][i, "MSE"]) - unlist(errors[[1]][i, "MSE"])
      delta[i,5] <- unlist(errors[[2]][i, "RMSE"]) - unlist(errors[[1]][i, "RMSE"])
    }
    if(type[i] == "c"){
      delta[i,(p-2)] <- unlist(errors[[2]][i, "CC"]) - unlist(errors[[1]][i, "CC"])
      delta[i,(p-1)] <- unlist(errors[[2]][i, "nCC"]) - unlist(errors[[1]][i, "nCC"])
      delta[i,p] <- unlist(errors[[2]][i, "CCmarg"]) - unlist(errors[[1]][i, "CCmarg"])
    }
  }
  delta[,1] <- errors[[1]][,1]
  colnames(delta) <- colnames(errors[[1]])

  # RETURN
  return(delta)
}

##### setup: Prepare data for model fitting
setup <- function(data, type, y = NULL, lags = NULL, scale = TRUE, allNum = FALSE){
  if(length(type) == 1 & ncol(data) > 1){type <- rep(type, ncol(data))}
  data <- data.frame(data)
  for(i in which(type == "c")){
    data[,i] <- as.factor(data[,i])
    levels(data[,i]) <- paste(1:length(levels(data[,i])))
  }
  if(scale == TRUE){for(i in which(type == "g")){data[,i] <- scale(data[,i])}}
  data <- data.frame(data)
  names(data) <- paste0("V", 1:ncol(data), ".")
  # PUTTING THIS IN JUST FOR BELOW
  lagIt <- function(data, y, lags = 1){
    if(max(lags) > 1){
      lagged <- list()
      for(i in 1:length(lags)){
        lagged[[i]] <- data[-c((nrow(data) - lags[i] + 1):nrow(data)),]
        names(lagged[[i]]) <- paste0(names(lagged[[i]]), "lag", lags[i], ".")
      }
      check <- sapply(lagged, nrow)
      if(any(check != check[length(check)])){
        check2 <- which(check != check[length(check)])
        for(j in 1:length(check2)){
          lagged[[j]] <- lagged[[j]][-c(1:(nrow(lagged[[j]]) - nrow(lagged[[length(lagged)]]))),]
          rownames(lagged[[j]]) <- c(1:nrow(lagged[[j]]))
        }
      }
      lagged <- do.call(cbind, lagged)
    } else {
      lagged <- data[-nrow(data),]
      names(lagged) <- paste0(names(lagged), "lag1.")
    }
    response <- data[-c(1:max(lags)),]
    new <- data.frame(y = response[,y], lagged)
    new
  }
  if(!is.null(lags)){data <- lagIt(data = data, y = y, lags = lags)}
  if(allNum == TRUE & scale == TRUE){
    for(j in 1:ncol(data)){data[,j] <- as.numeric(data[,j])}
  }
  data
}
