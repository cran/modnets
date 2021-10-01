#' Variable selection for moderated networks
#'
#' Perform variable selection via the LASSO, best subsets selection, forward
#' selection, backward selection, or sequential replacement on unmoderated
#' networks. Or, perform variable selection via the hierarchical LASSO for
#' moderated networks. Can be used for both GGMs and SUR networks.
#'
#' The primary value of the output is to be used as input when fitting the
#' selected model with the \code{\link{fitNetwork}} function. Specifically, the
#' output of \code{\link{varSelect}} can be assigned to the \code{type} argument
#' of \code{\link{fitNetwork}} in order to fit the constrained models that were
#' selected across nodes.
#'
#' For moderated networks, the only variable selection approach available is
#' through the \code{glinternet} package, which implements the hierarchical
#' LASSO. The criterion for model selection dictates which function from the
#' package is used, where information criteria use the
#' \code{\link[glinternet:glinternet]{glinternet::glinternet}} function to
#' compute models, and cross-validation calls the
#' \code{\link[glinternet:glinternet.cv]{glinternet::glinternet.cv}} function.
#'
#' @param data \code{n x k} dataframe or matrix.
#' @param m Character vector or numeric vector indicating the moderator(s), if
#'   any. Can also specify \code{"all"} to make every variable serve as a
#'   moderator, or \code{0} to indicate that there are no moderators. If the
#'   length of \code{m} is \code{k - 1} or longer, then it will not be possible
#'   to have the moderators as exogenous variables. Thus, \code{exogenous} will
#'   automatically become \code{FALSE}.
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
#' @param lags Numeric or logical. Can only be 0, 1 or \code{TRUE} or
#'   \code{FALSE}. \code{NULL} is interpreted as \code{FALSE}. Indicates whether
#'   to fit a time-lagged network or a GGM.
#' @param exogenous Logical. Indicates whether moderator variables should be
#'   treated as exogenous or not. If they are exogenous, they will not be
#'   modeled as outcomes/nodes in the network. If the number of moderators
#'   reaches \code{k - 1} or \code{k}, then \code{exogenous} will automatically
#'   be \code{FALSE}.
#' @param type Determines whether to use gaussian models \code{"g"} or binomial
#'   models \code{"c"}. Can also just use \code{"gaussian"} or
#'   \code{"binomial"}. Moreover, a vector of length \code{k} can be provided
#'   such that a value is given to every variable. Ultimately this is not
#'   necessary, though, as such values are automatically detected.
#' @param center Logical. Determines whether to mean-center the variables.
#' @param scale Logical. Determines whether to standardize the variables.
#' @param gamma Numeric value of the hyperparameter for the \code{"EBIC"}
#'   criterion. Only relevant if \code{criterion = "EBIC"}. Recommended to use a
#'   value between 0 and .5, where larger values impose a larger penalty on the
#'   criterion.
#' @param nfolds Only relevant if \code{criterion = "CV"}. Determines the number
#'   of folds to use in cross-validation.
#' @param varSeed Numeric value providing a seed to be set at the beginning of
#'   the selection procedure. Recommended for reproducible results.
#' @param useSE Logical. Only relevant if \code{method = "glinternet"} and
#'   \code{criterion = "CV"}. Indicates whether to use the standard error of the
#'   estimates across folds, if \code{TRUE}, or to use the standard deviation,
#'   if \code{FALSE}.
#' @param nlam if \code{method = "glinternet"}, determines the number of lambda
#'   values to evaluate in the selection path.
#' @param covs Numeric or character string indicating a variable to be used as a
#'   covariate. Currently not working properly.
#' @param verbose Logical. Determines whether to provide output to the console
#'   about the status of the procedure.
#' @param beepno Character string or numeric value to indicate which variable
#'   (if any) encodes the survey number within a single day. Must be used in
#'   conjunction with \code{dayno} argument.
#' @param dayno Character string or numeric value to indicate which variable (if
#'   any) encodes the survey number within a single day. Must be used in
#'   conjunction with \code{beepno} argument.
#'
#' @return List of all models, with the selected variables for each along with
#'   model coefficients and the variable selection models themselves. Primarily
#'   for use as input to the \code{type} argument of the
#'   \code{\link{fitNetwork}} function.
#' @export
#'
#' @seealso \code{\link{resample}, \link{fitNetwork}, \link{bootNet},
#'   \link{mlGVAR}, \link[glinternet:glinternet]{glinternet::glinternet},
#'   \link[glinternet:glinternet.cv]{glinternet::glinternet.cv},
#'   \link[glmnet:glmnet]{glmnet::glmnet},
#'   \link[glmnet:cv.glmnet]{glmnet::cv.glmnet},
#'   \link[leaps:regsubsets]{leaps::regsubsets}}
#'
#' @examples
#' \donttest{
#' vars1 <- varSelect(ggmDat, criterion = 'BIC', method = 'subset')
#' fit1 <- fitNetwork(ggmDat, type = vars1)
#'
#' vars2 <- varSelect(ggmDat, criterion = 'CV', method = 'glmnet')
#' fit2 <- fitNetwork(ggmDat, type = vars2, which.lam = 'min')
#'
#' # Add a moderator
#' vars3 <- varSelect(ggmDat, m = 'M', criterion = 'EBIC', gamma = .5)
#' fit3 <- fitNetwork(ggmDat, moderators = 'M', type = vars3)
#' }
varSelect <- function(data, m = NULL, criterion = "AIC", method = "glmnet",
                      lags = NULL, exogenous = TRUE, type = "g", center = TRUE,
                      scale = FALSE, gamma = .5, nfolds = 10, varSeed = NULL,
                      useSE = TRUE, nlam = NULL, covs = NULL, verbose = TRUE,
                      beepno = NULL, dayno = NULL){
  dat <- data
  ALL <- FALSE
  dmnames <- NULL # VARSELECT START
  mall <- which(sapply(c(0, 'all'), identical, as.character(m)))
  if(any(mall)){
    m <- switch(mall, NULL, 1:ncol(dat))
    if(mall == 2){ALL <- TRUE}
  } else if(isTRUE(is.character(m))){
    m <- which(colnames(data) %in% m)
  }
  if(is.null(lags)){
    if(length(m) >= ncol(dat) - 1){exogenous <- FALSE}
    if(!exogenous){ALL <- TRUE}
  } else if(any(!sapply(c(beepno, dayno), is.null))){
    stopifnot(!is.null(beepno) & !is.null(dayno))
    dat <- getConsec(data = dat, beepno = beepno, dayno = dayno)
  }
  if(!is(dat, 'list')){
    dat <- structure(data.frame(dat), samp_ind = attr(data, "samp_ind"))
    if(is.null(m) | ALL | !is.null(lags)){
      if(!is.null(lags)){
        vs <- colnames(dat)
        if(!is.null(m) & !ALL){mname <- vs[m]}
        if(!is.null(covs)){
          if(is.character(covs)){covs <- which(colnames(dat) %in% covs)}
          dat <- dat[, -covs]
          if(!is.null(m) & !ALL){m <- which(colnames(dat) %in% mname)}
        }
        dmnames <- colnames(lagMat(data = dat, m = m)$X)
        dat <- lagMat(data = dat, type = type, center = center, scale = scale)
        dat$full <- cbind.data.frame(dat$X, dat$Y) # FULLFIX
        if(ALL | (!is.null(m) & all(m == 0))){exogenous <- FALSE; m <- 1:ncol(dat$Y)}
        #if(ALL | is.null(m)){exogenous <- FALSE; m <- 1:ncol(dat$Y)}
        if(!is.null(m)){
          if(exogenous & length(m) < ncol(dat$Y)){
            dat$Y <- dat$Y[, -m]
            dat$full <- dat$full[, -m]
            if(!is(dat$Y, 'matrix')){ # NEW
            #if(class(dat$Y) != "matrix"){
              dat$Y <- as.matrix(dat$Y, ncol = 1)
              colnames(dat$Y) <- colnames(dat$full)[1]
            }
          } else if(length(m) == ncol(dat$Y)){
            mname <- gsub("[.]y$", "", colnames(dat$Y))
            exogenous <- FALSE
            ALL <- TRUE
          }
        }
      } else {
        mname <- colnames(dat)[m]
        dat <- list(Y = dat, X = dat)
      }
    } else if(length(m) >= ncol(dat) - 1){
      mname <- colnames(dat)[m]
      exogenous <- FALSE; ALL <- TRUE
    } else if(class(m) == "list"){
      dat <- list(Y = dat, X = data.frame(dat, m))
      m <- ncol(dat$Y)
    } else if(class(m) %in% c("numeric", "integer")){
      dn <- colnames(dat)
      dat <- list(Y = dat[, -m], X = data.frame(dat[, -m], dat[, m]))
      colnames(dat$X) <- c(dn[-m], dn[m])
      m <- ncol(dat$Y)
    }
  }
  if(!is.null(m)){
    if(class(m) == "list"){
      stopifnot(!ALL)
      mname <- names(m)
      m <- ncol(dat$Y)
      if(m < (ncol(dat$X) - 1)){
        covariates <- TRUE
        dn <- colnames(dat$X)
        m <- which(dn == mname)
        dat$X <- data.frame(dat$X[, -m], dat$X[, m])
        colnames(dat$X) <- c(dn[-m], dn[m])
        m <- ncol(dat$X) - 1
      }
    }
    if(all(m != 0)){m0 <- m}
    if(!method %in% c('hiernet', 'glinternet')){method <- 'glinternet'}
  }
  if(!is.null(varSeed)){set.seed(varSeed)}
  p <- ncol(dat$Y); data <- dat$X; Y <- dat$Y
  method <- match.arg(tolower(method), c(
    "hiernet", "glinternet", "subset", "backward",
    "forward", "seqrep", "glmnet", "lasso"))
  criterion <- toupper(match.arg(tolower(criterion), c(
    "cv", "aic", "bic", "ebic", "cp", "rss", "adjr2", "rsq", "r2")))
  ### VARSELECT
  if(method %in% c("hiernet", "glinternet")){
    if(is.null(nlam)){nlam <- ifelse(method == "hiernet", 20, 50)}
    hierMods <- list(); t <- c()
    for(i in 1:p){
      if(all(colnames(dat$Y) %in% colnames(dat$X))){
        data <- cbind(y = Y[, i], dat$X[, -i])
      } else {
        data <- cbind(y = Y[, i], dat$X)
      }
      if(verbose == TRUE){cat("Fitting model ", i, "/", p, "...", sep = "")}
      if(verbose == "pbar"){if(i == 1){pb <- txtProgressBar(min = 0, max = p, style = 2, width = 43)}}
      if(ALL & !is.null(m)){
        if(all(m != 0)){
          m <- switch(2 - (i %in% m), NULL, which(colnames(data)[-1] %in% mname))
        }
      }
      tx <- Sys.time()
      hierMods[[i]] <- fitHierLASSO(data = data, yvar = i, type = type, m = m,
                                    nlam = nlam, nfolds = nfolds, gamma = gamma,
                                    method = method, useSE = useSE, lags = lags,
                                    criterion = criterion, dmnames = dmnames,
                                    verbose = ifelse(is.logical(verbose), verbose, FALSE))
      t[i] <- tx <- Sys.time() - tx
      names(t)[i] <- attr(tx, "units")
      if(ALL & exists("m0", inherits = FALSE)){m <- m0}
      if(verbose == TRUE){
        cat("  Complete! ", "(", round(t[i], 2), " ", attr(tx, "units"), ")\n", sep = "")}
      if(verbose == "pbar"){setTxtProgressBar(pb, i); if(i == p){close(pb)}}
    }
    if(verbose == TRUE){
      if(length(unique(names(t))) == 1){
        cat("####### Total time:", round(sum(t), 2), names(t)[1], "\n\n")
      } else {
        tt <- t
        if(length(unique(names(tt))) == 2){
          if(all(sort(unique(names(tt))) == c("mins", "secs"))){
            tt[names(tt) == "mins"] <- 60 * tt[names(tt) == "mins"]
            cat("####### Total time:", round(sum(tt)/60, 2), "mins\n\n")
          } else if(all(sort(unique(names(tt))) == c("hours", "mins"))){
            tt[names(tt) == "hours"] <- 60 * tt[names(tt) == "hours"]
            cat("####### Total time:", round(sum(tt)/60, 2), "hours\n\n")
          }
        } else if(all(sort(unique(names(tt))) == c("hours", "mins", "secs"))){
          tt[names(tt) == "hours"] <- 360 * tt[names(tt) == "hours"]
          tt[names(tt) == "mins"] <- 60 * tt[names(tt) == "mins"]
          cat("####### Total time:", round(sum(tt)/360, 2), "hours\n\n")
        }
      }
    }
    names(hierMods) <- colnames(Y)
    attributes(hierMods)$method <- method
    attributes(hierMods)$criterion <- criterion
    if(criterion == "EBIC"){attributes(hierMods)$gamma <- gamma}
    attributes(hierMods)$time <- t
    if(exists("covariates", inherits = FALSE)){attributes(hierMods)$covariates <- TRUE}
    if(!is.null(covs)){attributes(hierMods)$covs <- covs}
    if(!is.null(lags)){
      attributes(hierMods)$moderators <- mname
      attributes(hierMods)$exogenous <- exogenous
    }
    return(hierMods)
  }
  ### LASSO SELECTION
  if(method %in% c("lasso", "glmnet")){
    lassoMods <- list(); t <- c()
    if(is.null(nlam)){nlam <- 100}
    for(i in 1:p){
      if(all(colnames(dat$Y) %in% colnames(dat$X))){
        data <- cbind(y = Y[, i], dat$X[, -i])
      } else {
        data <- cbind(y = Y[, i], dat$X)
      }
      if(verbose != FALSE){if(i == 1){pb <- txtProgressBar(min = 0, max = p, style = 2, width = 43)}}
      tx <- Sys.time()
      lassoMods[[i]] <- lassoSelect(data = data, yvar = i, criterion = criterion,
                                    type = type, gamma = gamma, nfolds = nfolds,
                                    nlam = nlam)
      t[i] <- tx <- Sys.time() - tx
      names(t)[i] <- attr(tx, "units")
      if(verbose != FALSE){setTxtProgressBar(pb, i); if(i == p){close(pb)}}
    }
    names(lassoMods) <- colnames(Y)
    attributes(lassoMods)[c("method", "criterion", "time")] <- list("glmnet", criterion, t)
    if(criterion %in% c("EBIC", "CV")){attr(lassoMods, "gamma") <- gamma}
    if(exists("covariates", inherits = FALSE)){attributes(lassoMods)$covariates <- TRUE}
    if(!is.null(covs)){attributes(lassoMods)$covs <- covs}
    return(lassoMods)
  }
  if(method %in% c("subset", "backward", "forward", "seqrep")){
    if(criterion == "CV"){criterion <- "Cp"}
    if(tolower(criterion) %in% c("aic", "ebic")){criterion <- "bic"}
    ind <- match.arg(tolower(criterion), c("cp", "bic", "adjr2", "rsq", "r2", "rss"))
    if(method == "subset"){method <- "exhaustive"}
    if(ind == "r2"){ind <- "rsq"}
    best <- ifelse(ind %in% c("cp", "bic", "rss"), which.min, which.max)
    regMods <- bestMods <- list()
    for(i in 1:p){
      if(all(colnames(dat$Y) %in% colnames(dat$X))){data <- dat$X[,-i]}
      regMods[[i]] <- summary(leaps::regsubsets(data, Y[,i], nvmax = ncol(data), method = method))
      bestMods[[i]] <- ifelse(regMods[[i]]$which, 1, 0)[best(regMods[[i]][[ind]]), -1]
    }
    bestMods <- lapply(bestMods, function(z) names(z)[z == 1])
    bestMods <- lapply(1:p, function(z){
      bestOut <- list(mod0 = bestMods[[z]], fitobj = regMods[[z]])
      attr(bestOut, "family") <- "g"
      return(bestOut)
    })
    names(bestMods) <- colnames(Y)
    attributes(bestMods)$method <- "regsubsets"
    attributes(bestMods)$criterion <- ind
    return(bestMods)
  }
}

##### lassoSelect: performs variable selection using the LASSO (glmnet)
lassoSelect <- function(data, yvar, type = "g", criterion = "EBIC",
                        gamma = .5, nfolds = 10, nlam = 100, alpha = 1){
  if(any(grepl(":", colnames(data)))){data <- data[,!grepl(":", colnames(data))]}
  criterion <- match.arg(criterion, c("CV", "EBIC", "BIC", "AIC"))
  y <- as.numeric(data[, 1])
  x <- data <- as.matrix(data[, -1])
  if(length(type) > 1){type <- type[yvar]}
  fam <- ifelse(type %in% c("g", "gaussian"), "gaussian", "binomial")
  if(criterion == "CV"){
    fit <- glmnet::cv.glmnet(x = x, y = y, family = fam, type.measure = "deviance",
                             nfolds = nfolds, nlambda = nlam, alpha = alpha)
  } else {
    fit <- glmnet::glmnet(x = x, y = y, family = fam, alpha = alpha, nlambda = nlam)
  }
  getIndices <- function(fit, y, x, fam = "gaussian", criterion = "EBIC",
                         gamma = .5, lam = "null", keepFit = FALSE){
    n <- length(y); p <- ncol(x)
    fam <- match.arg(fam, c("gaussian", "binomial", "multinomial"))
    lam <- match.arg(lam, c("null", "lambda.min", "lambda.1se"))
    if(criterion != "CV"){n_lambdas <- length(fit$lambda)}
    if(fam == "gaussian"){
      if(criterion != "CV"){
        beta0 <- matrix(coef(fit, s = 1)[1], ncol = 1)
        yhat <- rep(1, n) * as.vector(beta0)
        n_neighbors <- sapply(1:n_lambdas, function(z){
          colSums(as.matrix(coef(fit, s = fit$lambda[z])[-1,]) != 0)
        })
        LL_model <- sum(dnorm(y, mean = yhat, sd = sqrt(sum((y - yhat)^2)/n), log = TRUE))
      } else {
        mods <- lapply(c("lambda.min", "lambda.1se"), function(z){
          betas <- matrix(coef(fit, s = z), ncol = 1)
          yhat <- cbind(1, x) %*% as.vector(betas)
          n_neighbors <- colSums(matrix(coef(fit, s = z)[-1, ], ncol = 1) != 0)
          LL_model <- sum(dnorm(y, mean = yhat, sd = sqrt(sum((y - yhat)^2)/n), log = TRUE))
          ic_lambda <- -2 * LL_model + n_neighbors * log(n) + 2 * gamma * n_neighbors * log(p)
          return(list(betas = matrix(betas[-1, ], ncol = 1), EBIC = ic_lambda))
        })
      }
    } else if(fam %in% c("multinomial", "binomial")){
      lam <- ifelse(criterion == "CV", list(c("lambda.min", "lambda.1se")), list(1))[[1]]
      mods <- lapply(lam, function(lam0){
        cats <- unique(y)
        n_cats <- length(cats)
        m_respdum <- matrix(NA, n, n_cats)
        m_coefs <- matrix(NA, n, n_cats)
        m_LL_parts <- matrix(NA, nrow = n, ncol = n_cats + 1)
        X <- cbind(rep(1, n), x)
        for(catIter in 1:n_cats){
          m_respdum[, catIter] <- (y == cats[catIter]) * 1
          if(fam == "multinomial"){
            m_coefs[, catIter] <- X %*% matrix(coef(fit, s = lam0)[[catIter]], ncol = 1)
          } else {
            m_coefs[, catIter] <- X %*% ((matrix(coef(fit, s = lam0), ncol = 1) * ifelse(catIter == 1, -1, 1))/2)
          }
          m_LL_parts[, catIter] <- m_respdum[, catIter] * m_coefs[, catIter]
        }
        m_LL_parts[, (n_cats + 1)] <- -log(rowSums(exp(m_coefs)))
        LL_model <- sum(rowSums(m_LL_parts))
        if(lam0 == 1){
          n_lambdas <- length(fit$lambda)
          n_neighbors <- c()
          for(NN in 1:n_lambdas){
            coefs_bin <- vector("list", length = n_cats)
            for(ca in 1:n_cats){
              if(fam == "multinomial"){
                coefs_bin[[ca]] <- as.matrix(coef(fit, s = fit$lambda[NN])[[ca]][-1, ]) != 0
              } else {
                coefs_bin[[ca]] <- as.matrix(coef(fit, s = fit$lambda[NN])) != 0
              }
            }
            n_neighbors[NN] <- colSums(Reduce("+", coefs_bin) != 0)
            if(fam == "binomial"){n_neighbors[NN] <- n_neighbors[NN] - 1}
          }
          return(list(LL_model = LL_model, n_neighbors = n_neighbors))
        } else {
          coefs_bin <- vector("list", length = n_cats)
          if(fam == "multinomial"){
            for(ca in 1:n_cats){coefs_bin[[ca]] <- as.matrix(coef(fit, s = lam0)[[ca]][-1, ]) != 0}
          } else {
            for(ca in 1:n_cats){coefs_bin[[ca]] <- as.matrix(coef(fit, s = lam0)[-1, ]) != 0}
          }
          n_neighbors <- colSums(Reduce("+", coefs_bin) != 0)
          ic_lambda <- -2 * LL_model + n_neighbors * log(n) + 2 * gamma * n_neighbors * log(p)
          betas <- matrix(coef(fit, s = lam0), ncol = 1)
          return(list(betas = matrix(betas[-1, ], ncol = 1), EBIC = ic_lambda))
        }
      })
      if(criterion != "CV"){
        LL_model <- mods[[1]]$LL_model
        n_neighbors <- mods[[1]]$n_neighbors
      }
    }
    if(criterion != "CV"){
      LL_sat <- 1/2 * fit$nulldev + LL_model
      deviance <- (1 - fit$dev.ratio) * fit$nulldev
      LL_lambda_models <- -1/2 * deviance + LL_sat
      ic_lambda <- -2 * LL_lambda_models + n_neighbors * ifelse(
        criterion == "AIC", 2, log(n)) + ifelse(
          criterion == "EBIC", list(2 * gamma * n_neighbors * log(p)), list(0))[[1]]
      allCoefs <- lapply(seq_len(n_lambdas), function(z) coef(fit)[, z])
      betas <- allCoefs[[which.min(ic_lambda)]][-1]
      coefs <- Matrix::Matrix(betas, sparse = TRUE)
      rownames(coefs) <- names(betas)
      fitobj <- list(fit = fit, fit0 = NA, crit = ic_lambda)
      if(keepFit){fitobj$fit0 <- glmnet::glmnet(x, y, fam, lambda = fit$lambda[which.min(ic_lambda)])}
      names(fitobj)[3] <- criterion
      output <- list(mod0 = names(betas)[betas != 0], coefs = coefs,
                     fitobj = fitobj, allCoefs = allCoefs)
      if(length(output$mod0) == 0){output$mod0 <- 1}
    } else {
      coefs <- Matrix::Matrix(do.call(cbind, lapply(mods, '[[', "betas")), sparse = TRUE)
      rownames(coefs) <- colnames(x)
      colnames(coefs) <- paste0("mod", c("0", "1se"))
      fitobj <- list(fitCV = fit, fit0 = NA, fit1se = NA)
      if(keepFit){
        fitobj$fit0 <- glmnet::glmnet(x, y, fam, lambda = fit$lambda.min)
        fitobj$fit1se <- glmnet::glmnet(x, y, fam, lambda = fit$lambda.1se)
      }
      attr(fitobj$fit0, "EBIC") <- mods[[1]]$EBIC
      attr(fitobj$fit1se, "EBIC") <- mods[[2]]$EBIC
      output <- list(mod0 = colnames(x)[coefs[, 1] != 0],
                     mod1se = colnames(x)[coefs[, 2] != 0],
                     coefs = coefs, fitobj = fitobj)
      if(length(output$mod0) == 0){output$mod0 <- 1}
      if(length(output$mod1se) == 0){output$mod1se <- 1}
    }
    attr(output, "family") <- fam
    return(output)
  }
  out <- getIndices(fit, y, x, fam, criterion, gamma)
  return(out)
}

##### fitHierLASSO: performs variable selection using the hierarchical LASSO
fitHierLASSO <- function(data, yvar, type = "g", m = NULL, criterion = "CV",
                         method = "glinternet", gamma = .5, nfolds = 10,
                         nlam = 50, lags = NULL, useSE = TRUE, diag = FALSE,
                         outMsgs = FALSE, dmnames = NULL, verbose = TRUE){
  if(any(grepl(":", colnames(data)))){data <- data[,!grepl(":", colnames(data))]}
  method <- match.arg(tolower(method), c("hiernet", "glinternet"))
  criterion <- match.arg(criterion, c("CV", "EBIC", "BIC", "AIC"))
  y <- as.numeric(data[, 1])
  x <- data <- as.matrix(data[, -1])
  if(method == "hiernet"){
    out1 <- capture.output({fitPath <- hierNet::hierNet.path(x, y, nlam = nlam, strong = TRUE, diagonal = diag)})
    out2 <- capture.output({fitCV <- hierNet::hierNet.cv(fitPath, x, y, nfolds = nfolds)})
    out3 <- capture.output({fit0 <- hierNet::hierNet(x, y, lam = fitCV$lamhat, strong = TRUE, diagonal = diag)})
    out4 <- capture.output({fit1se <- hierNet::hierNet(x, y, lam = fitCV$lamhat.1se, strong = TRUE, diagonal = diag)})
    mod0 <- c(fit0$bp - fit0$bn, fit0$th[lower.tri(fit0$th)])
    mod1se <- c(fit1se$bp - fit1se$bn, fit1se$th[lower.tri(fit1se$th)])
    coefs <- Matrix::Matrix(cbind(mod0, mod1se), sparse = TRUE)
  } else if(method == "glinternet"){
    if(length(type) > 1){type <- type[yvar]}
    fam <- ifelse(type %in% c("g", "gaussian"), "gaussian", "binomial")
    type <- rep(1, ncol(x))
    if(criterion == "CV"){
      fitCV <- tryCatch({glinternet::glinternet.cv(x, y, type, nFolds = nfolds, nLambda = nlam,
                                                   interactionCandidates = m, family = fam)},
                        error = function(e){
                          failed <- TRUE
                          take <- 1
                          if(verbose){cat("\n")}
                          while(failed == TRUE){
                            if(take <= 5){
                              if(verbose){cat("  Failed.. trying again, take =", take, "\n")}
                              fitCV <- try(glinternet::glinternet.cv(x, y, type, nFolds = nfolds, nLambda = nlam,
                                                                     interactionCandidates = m, family = fam), silent = TRUE)
                              if(class(fitCV) == "try-error"){
                                failed <- TRUE
                                take <- take + 1
                              } else {
                                failed <- FALSE
                              }
                            } else if(take <= 10){
                              if(verbose){cat("  Failed.. trying nlam = 20, take =", take, "\n")}
                              fitCV <- try(glinternet::glinternet.cv(x, y, type, nFolds = nfolds, nLambda = 20,
                                                                     interactionCandidates = m, family = fam), silent = TRUE)
                              if(class(fitCV) == "try-error"){
                                failed <- TRUE
                                take <- take + 1
                              } else {
                                failed <- FALSE
                              }
                            } else {
                              if(verbose){cat("  Failed.. trying nlam = 20 & nFolds = 3, take =", take, "\n")}
                              fitCV <- try(glinternet::glinternet.cv(x, y, type, nFolds = 3, nLambda = 20,
                                                                     interactionCandidates = m, family = fam), silent = TRUE)
                              if(class(fitCV) == "try-error"){
                                failed <- TRUE
                                take <- take + 1
                                if(take == 20){break}
                              } else {
                                failed <- FALSE
                              }
                            }
                          }
                          fitCV
                        })
      if(useSE == TRUE){
        lamlist <- fitCV$lambda
        errm <- fitCV$cvErr
        errse <- fitCV$cvErrStd <- fitCV$cvErrStd/sqrt(nfolds)
        o <- which.min(errm)
        lamhat <- lamlist[o]
        oo <- errm <= errm[o] + errse[o]
        fitCV$lambdaHat1Std <- lamlist[oo & lamlist >= lamhat][1]
      }
      which.lam0 <- which(fitCV$lambda == fitCV$lambdaHat)
      while(is.null(fitCV$glinternetFit$activeSet[[which.lam0]])){
        if(verbose){cat("\n  Mod0 empty.. choosing new lambda\n")}
        which.lam0 <- which.lam0 + 1
      }
      fitCV$lambdaHat <- fitCV$lambda[which.lam0]
      which.lam1se <- which(fitCV$lambda == fitCV$lambdaHat1Std)
      while(is.null(fitCV$glinternetFit$activeSet[[which.lam1se]])){
        if(verbose){cat("\n  Mod1SE empty.. choosing new lambda\n")}
        which.lam1se <- which.lam1se + 1
      }
      fitCV$lambdaHat1Std <- fitCV$lambda[which.lam1se]
      fit0 <- glinternet::glinternet(x, y, type, lambda = fitCV$lambdaHat,
                                     interactionCandidates = m, family = fam)
      fit1se <- glinternet::glinternet(x, y, type, lambda = fitCV$lambdaHat1Std,
                                       interactionCandidates = m, family = fam)
      attributes(fit1se)$useSE <- attributes(fitCV)$useSE <- useSE
      mod0 <- coef(fit0)[[2]]
      mod1se <- coef(fit1se)[[2]]
    } else {
      fit <- glinternet::glinternet(x, y, type, interactionCandidates = m,
                                    family = fam, nLambda = nlam)
      coefs <- coef(fit)[-1]
      mains <- 1:ncol(x)
      ints <- t(combn(mains, 2))
      ints2 <- as.numeric(apply(ints, 1, paste, collapse = ""))
      if(is.null(lags) & !is.null(m)){
        vs <- colnames(x)
        vs1 <- vs[m]
        if(length(vs1) > 1){vs1 <- paste0("(", paste(vs1, collapse = " + "), ")")}
        vs2 <- as.formula(paste0("~ . * ", vs1))
        dmnames <- colnames(model.matrix(vs2, data.frame(x)))[-1]
      }
      allCoefs <- lapply(coefs, function(z){
        zmain1 <- z$mainEffects$cont
        zmain2 <- z$mainEffectsCoef$cont
        if(length(zmain1) != 0){
          if(any(!mains %in% zmain1)){
            zmiss1 <- mains[!mains %in% zmain1]
            zcoefs1 <- c(zmain2, rep(0, length(zmiss1)))[order(c(zmain1, zmiss1))]
          } else {
            zcoefs1 <- zmain2[order(zmain1)]
          }
        } else {
          zcoefs1 <- rep(0, length(mains))
        }
        zint1 <- z$interactions$contcont
        zint2 <- z$interactionsCoef$contcont
        if(length(zint1) != 0){
          zints1 <- as.numeric(apply(zint1, 1, paste, collapse = ""))
          if(nrow(ints) != nrow(zint1)){
            zcoefs2 <- rep(0, nrow(ints))
            zcoefs2[which(ints2 %in% zints1)] <- zint2
          } else {
            zcoefs2 <- zint2[match(zints1, ints2)]
          }
        } else {
          zcoefs2 <- rep(0, nrow(ints))
        }
        betas <- unlist(c(zcoefs1, zcoefs2))
        names(betas) <- c(colnames(x), apply(combn(colnames(x), 2), 2, paste, collapse = ":"))
        if(!is.null(m)){
          betas <- betas[which(names(betas) %in% dmnames)]
          #if(is.null(lags)){
          #  x2 <- c(colnames(x), paste0(colnames(x)[-m], ":", colnames(x)[m]))
          #  betas <- betas[which(names(betas) %in% x2)]
          #} else {
          #  betas <- betas[which(names(betas) %in% dmnames)]
          #}
        }
        return(betas)
      })
      n_neighbors <- sapply(allCoefs, function(z) sum(z != 0))
      LL_models <- sapply(1:length(allCoefs), function(z){
        s2 <- sum((y - fit$fitted[, z + 1])^2)/length(y)
        sum(dnorm(y, mean = fit$fitted[, z + 1], sd = sqrt(s2), log = TRUE))
      })
      p <- length(mains) + nrow(ints)
      if(!is.null(m)){
        if(all(m == 0)){
          p <- length(mains)
        } else {
          p <- length(c(colnames(x), paste0(colnames(x)[-m], ":", colnames(x)[m])))
        }
      }
      ic_lambda <- -2 * LL_models + n_neighbors * ifelse(
        criterion == "AIC", 2, log(nrow(x))) + ifelse(
          criterion == "EBIC", list(2 * gamma * n_neighbors * log(p)), list(0))[[1]]
      betas <- allCoefs[[which.min(ic_lambda)]]
      lambda_min <- fit$lambda[which.min(ic_lambda) + 1]
      coefs <- Matrix::Matrix(betas, sparse = TRUE)
      rownames(coefs) <- names(betas)
      fitobj <- list(fit = fit, fit0 = NA, crit = ic_lambda)
      if(ifelse(!is.null(m), ifelse(all(m == 0), FALSE, TRUE), TRUE) & method == "hiernet"){
        fitobj$fit0  <- glinternet::glinternet(x, y, type, interactionCandidates = m,
                                               lambda = lambda_min, family = fam)
      }
      names(fitobj)[3] <- criterion
      output <- list(mod0 = names(betas)[betas != 0], coefs = coefs,
                     fitobj = fitobj, allCoefs = allCoefs)
      attr(output, "family") <- ifelse(fam == "gaussian", "g", "c")
      return(output)
    }
    mains <- 1:ncol(x)
    ints <- t(combn(mains, 2))
    ints2 <- as.numeric(apply(ints, 1, paste, collapse = ""))
    if(length(mod0$mainEffects$cont) != 0){
      if(any(!mains %in% mod0$mainEffects$cont)){
        mod0miss1 <- mains[!mains %in% mod0$mainEffects$cont]
        mod0coefs1 <- c(mod0$mainEffectsCoef$cont,
                        rep(0, length(mod0miss1)))[order(c(mod0$mainEffects$cont, mod0miss1))]
      } else {
        mod0coefs1 <- mod0$mainEffectsCoef$cont[order(mod0$mainEffects$cont)]
      }
    } else {
      mod0coefs1 <- rep(0, length(mains))
    }
    if(length(mod1se$mainEffects$cont) != 0){
      if(any(!mains %in% mod1se$mainEffects$cont)){
        mod1semiss1 <- mains[!mains %in% mod1se$mainEffects$cont]
        mod1secoefs1 <- c(mod1se$mainEffectsCoef$cont,
                          rep(0, length(mod1semiss1)))[order(c(mod1se$mainEffects$cont, mod1semiss1))]
      } else {
        mod1secoefs1 <- mod1se$mainEffectsCoef$cont[order(mod1se$mainEffects$cont)]
      }
    } else {
      mod1secoefs1 <- rep(0, length(mains))
    }
    if(length(mod0$interactions$contcont) != 0){
      mod0ints1 <- as.numeric(apply(mod0$interactions$contcont, 1, paste, collapse = ""))
      if(nrow(ints) != nrow(mod0$interactions$contcont)){
        mod0coefs2 <- rep(0, nrow(ints))
        mod0coefs2[which(ints2 %in% mod0ints1)] <- mod0$interactionsCoef$contcont
      } else {
        mod0coefs2 <- mod0$interactionsCoef$contcont[match(mod0ints1, ints2)]
      }
    } else {
      mod0coefs2 <- rep(0, nrow(ints))
    }
    if(length(mod1se$interactions$contcont) != 0){
      mod1seints1 <- as.numeric(apply(mod1se$interactions$contcont, 1, paste, collapse = ""))
      if(nrow(ints) != nrow(mod1se$interactions$contcont)){
        mod1secoefs2 <- rep(0, nrow(ints))
        mod1secoefs2[which(ints2 %in% mod1seints1)] <- mod1se$interactionsCoef$contcont
      } else {
        mod1secoefs2 <- mod1se$interactionsCoef$contcont[match(mod1seints1, ints2)]
      }
    } else {
      mod1secoefs2 <- rep(0, nrow(ints))
    }
    mod0 <- unlist(c(mod0coefs1, mod0coefs2))
    mod1se <- unlist(c(mod1secoefs1, mod1secoefs2))
    coefs <- Matrix::Matrix(cbind(mod0, mod1se), sparse = TRUE)
  }
  allNames <- c(colnames(data), apply(combn(colnames(data), 2), 2, paste, collapse = ":"))
  rownames(coefs) <- allNames
  if(outMsgs == TRUE & method == "hiernet"){
    output <- list(mod0 = allNames[mod0 != 0], mod1se = allNames[mod1se != 0], coefs = coefs,
                   fitobj = list(fitCV = fitCV, fit0 = fit0, fit1se = fit1se),
                   outMsgs = list(outPath = out1, outCV = out2, out0 = out3, out1se = out4))
  } else {
    output <- list(mod0 = allNames[mod0 != 0], mod1se = allNames[mod1se != 0], coefs = coefs,
                   fitobj = list(fitCV = fitCV, fit0 = fit0, fit1se = fit1se))
  }
  attr(output, "family") <- ifelse(method == "glinternet", ifelse(fam == "gaussian", "g", "c"), "g")
  return(output)
}
