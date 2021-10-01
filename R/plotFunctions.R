#' Plot moderated and unmoderated network models
#'
#' Core function for plotting various types of network models. Accessible
#' through the \code{plot()} S3 generic function.
#'
#' @param x Output from any of the \code{modnets} model fitting or simulation
#'   functions.
#' @param which.net When multiple networks exist for a single object, this
#'   allows the user to indicate which network to plot. For a GGM, all values of
#'   this argument return the same adjacency matrix. For a SUR network,
#'   \code{"beta"} and \code{"temporal"} plot the temporal network, while
#'   \code{"pdc"} plots the Partial Directed Correlations, or the standardized
#'   temporal network. \code{"contemporaneous"} and \code{"pcc"} plot the
#'   standardized contemporaneous network (Partial Contemporaneous
#'   Correlations). All of these terms apply for multilevel networks, but
#'   \code{"between"} can also plot the between-subjects network. Additionally,
#'   the value \code{"coef"} will plot the model coefficients and confidence
#'   intervals, defaulting to the \code{\link{plotCoefs}} function. Moreover,
#'   with GGMs or outputs from \code{\link{mlGVAR}} with a moderated
#'   between-subjects network, the value \code{"ints"} will call the
#'   \code{\link{intsPlot}} function. If a numeric or logical value is supplied,
#'   however, this argument will function as the \code{threshold} argument. A
#'   numeric value will set a threshold at the supplied value, while \code{TRUE}
#'   will set a threshold of .05.
#' @param threshold A numeric or logical value to set a p-value threshold.
#'   \code{TRUE} will automatically set the threshold at .05.
#' @param layout Character. Corresponds to the \code{layout} argument in the
#'   \code{\link[qgraph:qgraph]{qgraph::qgraph}} function.
#' @param predict If \code{TRUE}, then prediction error associated with each
#'   node will be plotted as a pie graph around the nodes. For continuous
#'   variables, the type of prediction error is determined by the \code{con}
#'   argument. For categorical variables, the type of error is determined by the
#'   \code{cat} argument. The desired value of \code{con} or \code{can} can be
#'   supplied directly into the present argument as well. Alternatively, another
#'   network model constituted by the same nodes can be supplied in order to
#'   plot the difference in prediction error, such as R-squared change.
#' @param mnet Logical. If \code{TRUE}, the moderator will be plotted as a
#'   square "node" in the network, along with main effects represented as
#'   directed edges.
#' @param names If \code{TRUE}, then the variable names associated with the
#'   model will be plotted as labels on the nodes. If \code{FALSE}, then nodes
#'   will be labeled with numbers rather than names. Alternatively, a character
#'   vector can be provided to serve as custom labels for the nodes.
#' @param nodewise Only applies to GGMs. If \code{TRUE}, then nodewise edges
#'   will be plotted rather than the undirected averages of corresponding edges.
#' @param scale Logical. Only applies when \code{predict} does not equal
#'   \code{FALSE}. The value of this argument is sent to the
#'   \code{\link{predictNet}} function. This argument will be removed.
#' @param lag This argument will be removed. The function will automatically
#'   detect whether the network is based on time-lagged data.
#' @param con Character string indicating which type of prediction error to plot
#'   for continuous variables, if \code{predict} does not equal \code{FALSE}.
#'   Options are: \code{"R2", "adjR2", "MSE", "RMSE"}
#' @param cat Character string indicating which type of prediction error to plot
#'   for categorical variables, if \code{predict} does not equal \code{FALSE}.
#'   Options are: \code{"nCC", "CC", "CCmarg"}
#' @param covNet Logical. Only applies when a covariate is modeled. Allows the
#'   covariate to be plotted as a separate square "node".
#' @param plot Logical. If \code{FALSE}, then a \code{qgraph} object will be
#'   returned rather than plotted.
#' @param elabs Logical. If \code{TRUE}, the values of the edges will be plotted
#'   as labels on the edges.
#' @param elsize numeric
#' @param rule Only applies to GGMs (including between-subjects networks) when a
#'   threshold is supplied. The \code{"AND"} rule will only preserve edges when
#'   both corresponding coefficients have p-values below the threshold, while
#'   the \code{"OR"} rule will preserve an edge so long as one of the two
#'   coefficients have a p-value below the supplied threshold.
#' @param binarize Logical. If \code{TRUE}, the network will be plotted as an
#'   unweighted network. Only applies to GGMs.
#' @param mlty Logical. If \code{FALSE}, then moderated edges are displayed as
#'   solid lines. If \code{TRUE}, then moderated edges are shown as dashed
#'   lines.
#' @param mselect If the model contains more than one moderator, input the
#'   character string naming which moderator you would like the plot to reflect.
#'   Only affects which lines are dashed or solid. Not compatible with the
#'   \code{mnet} argument.
#' @param ... Additional arguments.
#'
#' @return Displays a network plot, or returns a \code{qgraph} object if
#'   \code{plot = FALSE}.
#' @export
#'
#' @seealso \code{\link{fitNetwork}, \link{predictNet}, \link{mlGVAR},
#'   \link{lmerVAR}, \link{simNet}, \link{mlGVARsim}, \link{plotCoefs},
#'   \link{intsPlot}, \link{resample}}
#'
#' @examples
#' fit1 <- fitNetwork(ggmDat)
#'
#' plot(fit1)
#' plotNet(fit1) # This and the command above produce the same result
#'
#' fit2 <- fitNetwork(gvarDat, moderators = 'M', lags = 1)
#'
#' plot(fit2, 'pdc') # Partial Directed Correlations
#' plot(fit2, 'pcc') # Partial Contemporaneous Correlations
plotNet <- function(x, which.net = 'temporal', threshold = FALSE, layout = 'spring',
                    predict = FALSE, mnet = FALSE, names = TRUE, nodewise = FALSE,
                    scale = FALSE, lag = NULL, con = 'R2', cat = 'nCC', covNet = FALSE,
                    plot = TRUE, elabs = FALSE, elsize = 1, rule = 'OR',
                    binarize = FALSE, mlty = TRUE, mselect = NULL, ...){
  getEdgeColors <- function(adjMat){
    obj <- sign(as.vector(adjMat))
    colMat <- rep(NA, length(obj))
    if(any(obj == 1)){colMat[obj == 1] <- "darkgreen"}
    if(any(obj == 0)){colMat[obj == 0] <- "darkgrey"}
    if(any(obj == -1)){colMat[obj == -1] <- "red"}
    colMat <- matrix(colMat, ncol = ncol(adjMat), nrow = nrow(adjMat))
    if(all(adjMat %in% 0:1)){colMat[colMat != 'darkgrey'] <- 'darkgrey'}
    dimnames(colMat) <- dimnames(adjMat)
    colMat
  }
  if(is.logical(which.net)){threshold <- which.net; which.net <- 'temporal'}
  if(is.numeric(which.net)){if(which.net < 1){threshold <- which.net; which.net <- 'temporal'}}
  if(isTRUE(covNet)){
    check <- isTRUE('covariates' %in% names(x$mods0))
    if(!check){stop('No covariates detected in fitted object')}
    object <- covNet(object = x, mnet = mnet, threshold = threshold)
    call <- replace(as.list(match.call())[-1], c('x', 'covNet', 'directed', 'mnet'),
                    list(x = object, covNet = FALSE, directed = object$d, mnet = FALSE))
    p <- attr(object, "cov")
    px <- ncol(object$data) - p
    call$shape <- c(rep("circle", px), rep("square", p))
    if(any(startsWith(names(call), 'border'))){
      border <- which(startsWith(names(call), 'border'))
      call <- append(call[-border], list(border.width = c(rep(1, px), rep(call[[border]], p))))
    }
    if('color' %in% names(call)){call$color <- c(rep('white', px), rep(call$color, p))}
    return(do.call(plotNet, call))
  }
  if(is.character(which.net) & length(which.net) == 1){
    if(startsWith(tolower(which.net), 'coef')){
      return(plotCoefs(x, plot = plot, ...))
    } else if(startsWith(tolower(which.net), 'int')){
      return(intsPlot(x, ...))
    }
  }
  if(any(class(x) %in% c('qgraph', 'bootnetResult', 'bootnet'))){
    if(is(x, 'bootnet')){x <- x$sample}
    if(is(x, 'bootnetResult')){
      if(x$default == 'graphicalVAR'){
        xx <- unname(which(sapply(c('t', 'c', 'b', 'p'), function(z) startsWith(which.net, z))))
        if(length(xx) == 0){stop('Must specify which.net')}
        which.net <- switch(xx, 'beta', 'kappa', 'beta', match.arg(toupper(which.net), c('PCC', 'PDC')))
        x <- x$results[[match(which.net, names(x$results))]]
        if(which.net == 'beta'){x <- x[, -1]}
      }
    }
    x <- getWmat(x)
  }
  if(is(x, 'bootNet') | (isTRUE(attr(x, 'resample')) & 'fit0' %in% names(x))){
    x <- x$fit0
  }
  if(is(x, 'ggmSim')){
    names(x)[grep('^trueNet$|^b1$', names(x))] <- 'adjMat'
    dimnames(x$adjMat) <- rep(list(1:ncol(x$adjMat)), 2)
    x$edgeColors <- getEdgeColors(x$adjMat)
    if('b2' %in% names(x)){x$modEdges <- 1 * (x$b2 != 0) + 1}
  }
  if(is(x, 'mgm')){x <- net(x, which.net, threshold)}
  atts <- names(attributes(x))
  if(!is(x, 'list') | 'GVARsim' %in% atts){ # NEW
  #if(!is(x, 'list') | 'simMLgvar' %in% atts){
    predict <- NULL; nodewise <- FALSE
    if(is(x, "mlGVARsim") | 'GVARsim' %in% atts){ # NEW
    #if(is(x, "mlVARsim") | 'simMLgvar' %in% atts){
      if('interactions' %in% names(x) & which.net %in% c('temporal', 'beta')){ # NEW
      #if('mm' %in% names(x) & which.net %in% c('temporal', 'beta')){
        nodewise <- TRUE
        #modEdges <- t(1 * (x$mm$mb2 != 0) + 1)
        modEdges <- t(1 * (x$interactions$mb2 != 0) + 1) # NEW
      }
      x <- t(net(x, which.net))
    }
    stopifnot(ncol(x) == nrow(x))
    if(isTRUE(names)){names <- colnames(x)}
    x <- list(adjMat = x, edgeColors = getEdgeColors(x))
    if(nodewise){
      nodewise <- FALSE
      x$modEdges <- modEdges
    }
  }
  if(any(c("SURnet", "mlGVAR", "lmerVAR") %in% atts)){
    if(is.numeric(which.net)){which.net <- c("t", "c", "b", "pdc")[which.net]}
    which.net <- match.arg(tolower(which.net), c(
      "temporal", "contemporaneous", "between", "beta", "pdc", "pcc"))
    which.net <- switch(which.net, pcc = "contemporaneous", beta = "temporal", which.net)
    if(isTRUE(attr(x, "mlGVAR"))){
      x <- x[[switch(which.net, between = "betweenNet", "fixedNets")]]
    }
    if("SURnet" %in% names(x)){
      if(!is.null(mselect) & length(x$call$moderators) > 1 & isTRUE(x$call$exogenous)){
        if(isTRUE(mselect)){mselect <- x$call$moderators[1]}
        adjm <- netInts(fit = x, n = 'temporal', threshold = TRUE,
                        avg = FALSE, rule = 'none', empty = FALSE,
                        mselect = mselect)
        x$SURnet$temporal$modEdges <- abs(sign(adjm)) + 1
      }
      x <- x$SURnet
    }
    if(!"adjMat" %in% names(x)){
      x[c("adjMat", "edgeColors")] <- eval(parse(text = paste0("x$", switch(
        which.net, pdc = "temporal$PDC", which.net))))[c("adjMat", "edgeColors")]
      if(!startsWith(which.net, "c") & "modEdges" %in% names(x$temporal)){
        x$modEdges <- x$temporal$modEdges
      }
    }
  }
  obmods <- x$call$moderators
  exo <- ifelse('exogenous' %in% names(x$call), x$call$exogenous, TRUE)
  if(isTRUE(attr(x, 'ggm')) & ifelse(
    all(!sapply(list(obmods, mselect), is.null)),
    length(obmods) > 1 & exo, FALSE)){
    if(isTRUE(mselect)){mselect <- obmods[1]}
    adjm <- netInts(fit = x, n = 'between', threshold = TRUE,
                    avg = !nodewise, rule = rule, empty = FALSE,
                    mselect = mselect)
    if(nodewise){
      x$nodewise$modEdgesNW <- abs(sign(adjm)) + 1
      diag(x$nodewise$modEdgesNW) <- 0
    } else {
      x$modEdges <- abs(sign(adjm)) + 1
      diag(x$modEdges) <- 0
    }
  }
  if(threshold != FALSE){
    if(!is.numeric(threshold)){threshold <- .05}
    if(mnet & "mnet" %in% names(x)){
      mn <- x$call$moderators
      adj1 <- net(x, "beta", threshold, rule)
      if(isTRUE(attr(x, "SURnet"))){
        rn <- gsub("[.]y$", "", rownames(x$mnet$adjMat))
        ind <- which(rn == mn)
        x$mnet$adjMat[-ind, -ind] <- adj1
        x$mnet$adjMat[-ind, ind] <- x$mnet$adjMat[-ind, ind] * ifelse(
          x$temporal$coefs$pvals[, mn] <= threshold, 1, 0)
      } else {
        ind <- nrow(x$mnet$adjMat)
        x$mnet$adjMat[-ind, -ind] <- adj1
        mps <- x$mods0$Bm
        mps2 <- matrix(0, nrow = nrow(adj1), ncol = 4)
        mps2[match(rownames(mps), rownames(adj1)), ] <- mps
        x$mnet$adjMat[ind, -ind] <- x$mnet$adjMat[ind, -ind] * ifelse(mps2[, 4] <= threshold, 1, 0)
        #x$mnet$adjMat[ind, -ind] <- x$mnet$adjMat[ind, -ind] * ifelse(x$mods0$Bm[, 4] <= threshold, 1, 0)
      }
    } else {
      if("ggm" %in% atts & isTRUE(nodewise)){
        x$nodewise$adjNW <- net(x, "beta", threshold, rule, TRUE)
      } else {
        x$adjMat <- net(x, which.net, threshold, rule)
      }
    }
  }
  if(isTRUE(names)){
    names <- gsub("[.]y$", "", rownames(x$adjMat))
    if(length(names) == 0){names <- 1:nrow(x$adjMat)}
  } else if(is.null(names) | all(names == FALSE)){
    names <- 1:nrow(x$adjMat)
  }
  names <- names[1:nrow(x$adjMat)]
  if(!is.null(lag)){warning('Lagged networks are automatically detected.')}
  lag <- NULL ### FUNCTION ARGUMENT WILL BE REMOVED
  if(is.null(lag) & "adjMats" %in% names(x)){
    stop("More than one lag modeled; need to specify which to plot")
  } else if(!"adjMats" %in% names(x)){
    if(any(grepl("lag", colnames(x$adjMat)))){lag <- 1}
  }
  if(!is.null(predict) & !mnet & !identical(predict, FALSE)){
    concat <- c('R2', 'adjR2', 'MSE', 'RMSE', 'nCC', 'CC', 'CCmarg')
    con <- match.arg(con, choices = c("R2", "adjR2", "MSE", "RMSE"))
    cat <- match.arg(cat, choices = c("nCC", "CC", "CCmarg"))
    type <- x$call$type
    if("ggm" %in% names(attributes(x))){
      type <- unname(sapply(x$fitobj, attr, "family"))
      if("gaussian" %in% type){type[type == "gaussian"] <- "g"}
      if("binomial" %in% type){type[type == "binomial"] <- "c"}
    }
    tt <- length(type)
    if(is(predict, 'list')){
      if(isTRUE(attr(predict, "mlGVAR"))){
        predict <- predict[[switch(which.net, between = "betweenNet", "fixedNets")]]
      }
      predict <- list(predictNet(predict, all = FALSE, scale = scale),
                      predictNet(x, all = FALSE, scale = scale))
      stopifnot(length(predict) == 2)
      pie <- list(); pieColor <- list()
      for(i in 1:tt){
        if(type[i] == "g"){
          pie[[i]] <- c(predict[[1]][i,con][[1]],
                        unlist(predict[[2]][i,con]) - unlist(predict[[1]][i,con]))
          pieColor[[i]] <- c("lightblue", "lightblue4")
          if(pie[[i]][1] < 0 & pie[[i]][2] < 0){
            pie[[i]][1] <- abs(pie[[i]][1])
            pie[[i]][2] <- abs(pie[[i]][2])
            pieColor[[i]] <- c("tomato", "tomato4")
          }
          if(pie[[i]][2] < 0){
            pie[[i]][1] <- pie[[i]][1] + pie[[i]][2]
            pie[[i]][2] <- abs(pie[[i]][2])
            pieColor[[i]][2] <- "tomato"
          }
          if(pie[[i]][1] < 0){
            pie[[i]][1] <- abs(pie[[i]][1])
            pie[[i]][2] <- pie[[i]][2] - pie[[i]][1]
            pieColor[[i]][1] <- "tomato"
          }
        }
        if(type[i] == "c"){
          pie[[i]] <- c(predict[[1]][i,cat][[1]],
                        unlist(predict[[2]][i,cat]) - unlist(predict[[1]][i,cat]))
          pieColor[[i]] <- c("peachpuff", "peachpuff4")
          if(pie[[i]][2] < 0){
            pie[[i]][1] <- pie[[i]][1] + pie[[i]][2]
            pie[[i]][2] <- abs(pie[[i]][2])
            pieColor[[i]][2] <- "tomato"
          }
        }
      }
    } else {
      if(is.character(predict)){
        predict <- toupper(intersect(tolower(predict[1]), tolower(concat)))
        if(length(predict) > 0){
          if(startsWith(predict, 'ADJ')){predict <- 'adjR2'}
          if(startsWith(predict, 'N')){predict <- 'nCC'}
          if(endsWith(predict, 'G')){predict <- 'CCmarg'}
          if(grepl('CC', predict)){cat <- predict} else {con <- predict}
        }
      }
      predict <- predictNet(x, all = FALSE, scale = scale)
      pie <- c(); pieColor <- c()
      for(i in 1:tt){
        if(type[i] == "g"){
          pie[i] <- predict[i,con]
          pieColor[i] <- "lightblue"
        }
        if(type[i] == "c"){
          pie[i] <- predict[i,cat]
          pieColor[i] <- "peachpuff"
        }
      }
      if(any(pie < 0)){
        pieColor[pie < 0] <- "tomato"
        pie[pie < 0] <- abs(pie[pie < 0][[1]])
      }
    }
  } else {
    pie <- NULL; pieColor <- NULL
  }
  if(!is.null(lag)){
    if(lag != 1 | lag == 1 & "adjMats" %in% names(x)){
      if(lag > length(x$adjMats)){
        lag <- which(sapply(lapply(lapply(x$adjMats, colnames),
                                   function(z) grep(paste0("lag", lag), z)), length) != 0)
      }
      qgraph(input = t(x$adjMats[[lag]]), layout = layout, labels = names,
             edge.color = t(x$edgeColors[[lag]]), edge.labels = elabs,
             edge.label.cex = elsize, DoNotPlot = !plot, pie = pie,
             pieColor = pieColor, ...)
    } else if(!mnet){
      if("modEdges" %in% names(x) & mlty){lty <- t(x$modEdges)} else {lty <- 1}
      qgraph(input = t(x$adjMat), layout = layout, edge.color = t(x$edgeColors),
             labels = names, DoNotPlot = !plot, pie = pie, pieColor = pieColor,
             edge.labels = elabs, lty = lty, edge.label.cex = elsize, ...)
    } else {
      if(class(names) %in% c("numeric", "integer")){
        names <- 1:ncol(x$mnet$adjMat)
      } else {
        names <- gsub(".lag1.", "", colnames(x$mnet$adjMat))
      }
      lty <- t(x$mnet$modEdges)
      qgraph(input = t(x$mnet$adjMat), layout = layout, labels = names,
             edge.color = t(x$mnet$edgeColors), DoNotPlot = !plot,
             pie = pie, pieColor = pieColor, shape = x$mnet$shape,
             edge.labels = elabs, lty = lty, edge.label.cex = elsize, ...)
    }
  } else {
    if(!nodewise){
      if(!mnet){
        if(binarize){x$adjMat[x$adjMat != 0] <- 1; x$edgeColors <- NA}
        if("modEdges" %in% names(x) & mlty){lty <- x$modEdges} else {lty <- 1}
        qgraph(input = x$adjMat, layout = layout, edge.color = x$edgeColors,
               labels = names, DoNotPlot = !plot, pie = pie, pieColor = pieColor,
               edge.labels = elabs, lty = lty, edge.label.cex = elsize, ...)
      } else {
        pp <- ncol(x$mnet$adjMat) - 1
        if(class(names) %in% c("numeric", "integer")){
          names <- 1:(pp + 1)
        } else {
          names <- c(names, colnames(x$mnet$adjMat)[pp + 1])
        }
        lty <- x$mnet$modEdges
        qgraph(input = x$mnet$adjMat, layout = layout, labels = names,
               edge.color = x$mnet$edgeColors, DoNotPlot = !plot, pie = pie,
               pieColor = pieColor, edge.labels = elabs, edge.label.cex = elsize,
               lty = lty, directed = x$mnet$d, shape = c(rep("circle", pp), "square"), ...)
      }
    } else {
      if("modEdgesNW" %in% names(x$nodewise) & mlty){lty <- x$nodewise$modEdgesNW} else {lty <- 1}
      qgraph(input = x$nodewise$adjNW, layout = layout, labels = names,
             edge.color = x$nodewise$edgeColsNW, DoNotPlot = !plot,
             pie = pie, pieColor = pieColor, edge.labels = elabs,
             edge.label.cex = elsize, lty = lty, ...)
    }
  }
}

# @describeIn plotNet For ggms!
#' @rdname plotNet
#' @export
plot.ggm <- function(x, which.net = 'temporal', threshold = FALSE, layout = 'spring',
                     predict = FALSE, mnet = FALSE, names = TRUE, nodewise = FALSE,
                     scale = FALSE, lag = NULL, con = 'R2', cat = 'nCC', covNet = FALSE,
                     plot = TRUE, elabs = FALSE, elsize = 1, rule = 'OR',
                     binarize = FALSE, mlty = TRUE, mselect = NULL, ...){
  args <- as.list(match.call())[-1]
  do.call(plotNet, args)
}

#' @rdname plotNet
#' @export
plot.SURnet <- function(x, which.net = 'temporal', threshold = FALSE, layout = 'spring',
                        predict = FALSE, mnet = FALSE, names = TRUE, nodewise = FALSE,
                        scale = FALSE, lag = NULL, con = 'R2', cat = 'nCC', covNet = FALSE,
                        plot = TRUE, elabs = FALSE, elsize = 1, rule = 'OR',
                        binarize = FALSE, mlty = TRUE, mselect = NULL, ...){
  args <- as.list(match.call())[-1]
  do.call(plotNet, args)
}

#' @rdname plotNet
#' @export
plot.mlGVAR <- function(x, which.net = 'temporal', threshold = FALSE, layout = 'spring',
                        predict = FALSE, mnet = FALSE, names = TRUE, nodewise = FALSE,
                        scale = FALSE, lag = NULL, con = 'R2', cat = 'nCC', covNet = FALSE,
                        plot = TRUE, elabs = FALSE, elsize = 1, rule = 'OR',
                        binarize = FALSE, mlty = TRUE, mselect = NULL, ...){
  args <- as.list(match.call())[-1]
  do.call(plotNet, args)
}

#' @rdname plotNet
#' @export
plot.lmerVAR <- function(x, which.net = 'temporal', threshold = FALSE, layout = 'spring',
                         predict = FALSE, mnet = FALSE, names = TRUE, nodewise = FALSE,
                         scale = FALSE, lag = NULL, con = 'R2', cat = 'nCC', covNet = FALSE,
                         plot = TRUE, elabs = FALSE, elsize = 1, rule = 'OR',
                         binarize = FALSE, mlty = TRUE, mselect = NULL, ...){
  args <- as.list(match.call())[-1]
  do.call(plotNet, args)
}

#' @rdname plotNet
#' @export
plot.ggmSim <- function(x, which.net = 'temporal', threshold = FALSE, layout = 'spring',
                        predict = FALSE, mnet = FALSE, names = TRUE, nodewise = FALSE,
                        scale = FALSE, lag = NULL, con = 'R2', cat = 'nCC', covNet = FALSE,
                        plot = TRUE, elabs = FALSE, elsize = 1, rule = 'OR',
                        binarize = FALSE, mlty = TRUE, mselect = NULL, ...){
  args <- as.list(match.call())[-1]
  do.call(plotNet, args)
}

#' @rdname plotNet
#' @export
plot.mlGVARsim <- function(x, which.net = 'temporal', threshold = FALSE, layout = 'spring',
                           predict = FALSE, mnet = FALSE, names = TRUE, nodewise = FALSE,
                           scale = FALSE, lag = NULL, con = 'R2', cat = 'nCC', covNet = FALSE,
                           plot = TRUE, elabs = FALSE, elsize = 1, rule = 'OR',
                           binarize = FALSE, mlty = TRUE, mselect = NULL, ...){
  args <- as.list(match.call())[-1]
  do.call(plotNet, args)
}

#' @rdname plotNet
#' @export
plot.GVARsim <- function(x, which.net = 'temporal', threshold = FALSE, layout = 'spring',
                         predict = FALSE, mnet = FALSE, names = TRUE, nodewise = FALSE,
                         scale = FALSE, lag = NULL, con = 'R2', cat = 'nCC', covNet = FALSE,
                         plot = TRUE, elabs = FALSE, elsize = 1, rule = 'OR',
                         binarize = FALSE, mlty = TRUE, mselect = NULL, ...){
  args <- as.list(match.call())[-1]
  do.call(plotNet, args)
}

#' Plot conditional networks at different levels of the moderator
#'
#' An easy wrapper for plotting the same network at different levels of a
#' moderator. Using the \code{mval} argument of the \code{\link{fitNetwork}}
#' function, you can create multiple models---conditional networks---wherein the
#' same model is fit at different values of the moderator.
#'
#' Importantly, this function will fix a common layout across all conditional
#' networks so that the network can be easily compared (visually) at different
#' levels of the moderator.
#'
#' @param nets List of network models fit with \code{\link{fitNetwork}}, where
#'   \code{mval} has been specified.
#' @param nodewise See corresponding argument in \code{\link{plotNet}}.
#' @param elsize Numeric value to indicate the size of the edge labels.
#' @param vsize Numeric value to indicate the size of the nodes. If \code{NULL},
#'   then a default value will be determined based on the number of nodes in the
#'   network.
#' @param elabs If \code{TRUE}, then edges will be labeled with their numeric
#'   values.
#' @param predict See corresponding argument in \code{\link{plotNet}}.
#' @param layout Can be a character string, corresponding to the options in
#'   \code{\link[qgraph:qgraph]{qgraph::qgraph}}, or can be a matrix that
#'   defines the layout (e.g., based on the
#'   \code{\link[qgraph:averageLayout]{qgraph::averageLayout}} function).
#'   Recommended to leave as \code{NULL}, so that the layout will be based on
#'   the list of networks provided.
#' @param which.net See corresponding argument in \code{\link{plotNet}}.
#' @param ... Additional arguments.
#'
#' @return Returns a plot where multiple conditional networks are plotted side
#'   by side.
#' @export
#'
#' @seealso \code{\link{fitNetwork}}
#'
#' @examples
#' data <- na.omit(psychTools::msq[, c('hostile', 'lonely', 'nervous', 'sleepy', 'depressed')])
#'
#' fit0 <- fitNetwork(data, moderators = 'depressed', mval = 0)
#' fit1 <- fitNetwork(data, moderators = 'depressed', mval = 1)
#' fit2 <- fitNetwork(data, moderators = 'depressed', mval = 2)
#'
#' fits <- list(fit0, fit1, fit2)
#' plotMods(fits)
plotMods <- function(nets, nodewise = FALSE, elsize = 2, vsize = NULL,
                     elabs = TRUE, predict = NULL, layout = NULL,
                     which.net = "temporal", ...){
  if(any(c("SURnet", "temporal") %in% unlist(lapply(nets, names)))){
    if("SURnet" %in% names(nets[[1]])){nets <- lapply(nets, '[[', "SURnet")}
    which.net <- match.arg(tolower(
      which.net), c("temporal", "contemporaneous", "pdc"))
    if(is.null(vsize)){vsize <- 10}
    nodewise <- FALSE
    ggm <- FALSE
  } else {
    ggm <- TRUE
  }
  if(is.null(vsize)){
    vsize <- (exp(-ncol(nets[[1]]$data)/80) * 8) + 1
    vsize <- vsize + (exp(-length(nets)/80) * 8)/length(nets)
  }
  getLayout <- function(x, which.net = "temporal"){
    stopifnot(class(x) == "list")
    averageLayout(lapply(x, function(z){
      plotNet(z, plot = FALSE, which.net = which.net)}))
  }
  getMax <- function(x, n = FALSE, which.net = "temporal"){
    if(!"ggm" %in% names(attributes(x[[1]]))){
      which.net <- match.arg(
        tolower(which.net), c("temporal", "contemporaneous", "pdc"))
      x <- lapply(x, '[[', ifelse(which.net == "pdc", "temporal", which.net))
      if(which.net == "pdc"){x <- lapply(x, '[[', "PDC")}
    }
    if(n == FALSE){
      max(unlist(lapply(x, function(z) max(abs(z$adjMat)))))
    } else {
      max(unlist(lapply(x, function(z) max(abs(z$nodewise$adjNW)))))
    }
  }
  if(is.null(layout)){layout <- getLayout(nets, which.net)}
  mx <- getMax(nets, nodewise, which.net)
  if(ggm){
    moderator <- capitalize(attr(nets[[1]], "moderator"))
    vals <- unname(sapply(nets, attr, "mval"))
  } else {
    moderator <- capitalize(nets[[1]]$call$moderators)
    vals <- unname(sapply(lapply(nets, '[[', "call"), '[[', "mval"))
  }
  layout(t(1:length(vals)))
  for(i in 1:length(vals)){
    plotNet(nets[[i]], predict = predict, layout = layout, elabs = elabs,
            elsize = elsize, nodewise = nodewise, maximum = mx,
            title = paste0(moderator, " = ", vals[i]),
            vsize = vsize, which.net = which.net, ...)
  }
}

#' Plot \code{bootNet} outputs
#'
#' Creates various types of plot to visualize \code{bootNet} objects.
#'
#' @param x Output from \code{\link{bootNet}}. Also some compatiblity with
#'   \code{resample} objects (when sampMethod != 'stability').
#' @param type The outcome measure to plot. Options include: \code{"edges",
#'   "strength", "ei", "outstrength", "instrength", "outei", "inei"}. The "out-"
#'   and "in-" options are only available for temporal networks. Moreover, both
#'   related options can be used together in temporal networks, by setting
#'   either \code{type = c("outstrength", "instrength")} or \code{type =
#'   c("outei", "inei")}.
#' @param net Determines which network to plot coefficients for. Options
#'   include: \code{"ggm", "temporal", "contemporaneous", "between"}. Only
#'   relevant to SUR networks or \code{mlGVAR} objects.
#' @param plot Primary use is to set as \code{"none"} or \code{FALSE} in order
#'   to return a table containing the constituents of the plot rather than
#'   create the plot itself. The options \code{"all"} and \code{"both"} each
#'   essentially indicate that both pairwise and interaction terms are plotted.
#'   Can also specify \code{"pairwise"} to only plot the pairwise terms, or
#'   \code{"interactions"} to only plot the interaction terms.
#' @param cor Numeric value to indicate the correlation stability value to be
#'   plotted. Only applies to the case-drop bootstrapping output.
#' @param order Determines how to arrange the predictors displayed in the plot.
#'   If \code{TRUE}, then defaults to \code{"mean"}. If \code{FALSE} then
#'   defaults to \code{"id"}. The \code{"mean"} option will arrange the values
#'   by the bootstrapped sample means. The \code{"sample"} option will arrange
#'   the values according to the statistics from the model fitted to the full
#'   sample. The \code{"id"} option will keep the variables in the same order
#'   that they appear in the dataframe. Not relevant to the case-drop bootstrap.
#' @param ci Numeric value between 0 and 1 to specify the confidence level.
#' @param pairwise Logical. Whether to plot pairwise relationships. Defaults to
#'   \code{TRUE}. If \code{FALSE}, this will override the \code{"all"} option of
#'   the \code{plot} argument.
#' @param interactions Logical. Whether to plot interactions. Defaults to
#'   \code{TRUE}. If \code{FALSE}, this will override the \code{"all"} option of
#'   the \code{plot} argument. Only relevant to moderated networks.
#' @param labels Logical. Determines whether to plot names of the predictors.
#' @param title Character vector the title label.
#' @param cis Either \code{"quantile"} or \code{"se"}. If \code{"quantile"},
#'   then confidence bands will be computed based on quantiles (specified by the
#'   \code{ci} argument) of the bootstrapped resamples. If \code{"se"}, then the
#'   confidence bands will be computed based on the standard errors associated
#'   with the sample statistics. Thus, the \code{"se"} argument will always
#'   produce a symmetric confidence band, whereas for \code{"quantile"} argument
#'   this is not necessary. Not relevant to outputs for the case-drop bootstrap.
#' @param true Defaults to \code{NULL}, not relevant for the case-drop
#'   bootstrap. Can supply another output from \code{\link{fitNetwork}}, or an
#'   adjacency matrix, to serve as the true network in the plot. If there are
#'   interactions in the model, then a \code{\link{fitNetwork}} object is
#'   recommended. Alternatively, this argument can be extremely useful for
#'   simulated data -- especially anything created with \code{\link{simNet}}.
#'   For whatever outcome (e.g., \code{edges, strength, EI}) is plotted,
#'   supplying another object to \code{true} will plot the values related to the
#'   true network, i.e., the data-generating model.
#' @param errbars Logical. Not relevant to the case-drop bootstrap. If
#'   \code{TRUE}, then error bars are used rather than confidence bands. Can be
#'   useful to home in on specific variables and see their confidence interval.
#' @param vline Logical or numeric. Not relevant to the case-drop bootstrap. If
#'   \code{TRUE}, then a dashed vertical line will be plotted at 0. If numeric,
#'   then the line will be plotted at the supplied intercept on the x-axis.
#' @param threshold Numeric or logical. Not relevant to the case-drop bootstrap.
#'   Has a significant effect on the bootstrapped coefficient distributions. If
#'   \code{TRUE}, then the default p-value threshold is set to .05. A numeric
#'   value can specify a different threshold. Causes the \code{\link{bootNet}}
#'   function to run the object again, only to re-compute the bootstrapped
#'   distributions after applying a p-value threshold to the results of each
#'   model iteration. If \code{NULL}, all coefficient estimates are used in
#'   estimating the posterior distribution of each parameter.
#' @param difference Logical. Not relevant to the case-drop bootstrap. If
#'   \code{TRUE}, then a difference plot is provided rather than a coefficient
#'   plot. In the difference plot, the diagonal squares reflect the fitted
#'   network coefficients for the the original sample. Black boxes indicate that
#'   the difference between the two edges, coefficients, or centrality values
#'   being compared is significantly different from 0. The significance level
#'   will have already been determined by the parameters used to fit the
#'   \code{bootNet} object. Gray boxes indicate the difference is not
#'   significantly different from 0.
#' @param color Logical. Only applies when \code{difference = TRUE}. Determines
#'   whether to add colors that reflect the sign of the sample values. Reflected
#'   in the diagonal of the difference plot.
#' @param text Logical. For difference plots, if \code{TRUE} then the statistics
#'   based on the full sample will be labeled in the diagonal boxes. For
#'   coefficient plots, setting this to \code{TRUE} will plot a label for each
#'   variable to reflect the proportion of times that it was selected across all
#'   bootstrapped iterations. Only relevant if a threshold was set for the
#'   fitted bootstrap models, either specified in the current function or was
#'   specified in creating the \code{\link{bootNet}} object. If a numeric value
#'   is provided, this will determine the size of the text label. Defaults to
#'   1.2 when \code{text = TRUE}.
#' @param textPos Supports the \code{text} argument for coefficient plots.
#'   Indicates the x-axis position of where to plot the coefficient labels.
#'   Generally will be numeric, but defaults to \code{"value"}, which means that
#'   the text will be labeled on top each point on the plot.
#' @param multi Useful when there are interactions in a model. If \code{TRUE},
#'   the single plot with a facet for both pairwise and interaction terms is
#'   split into two separate plots. Allows for a more elegant side-by-side plot,
#'   and allows arguments that are restricted for plots of either pairwise or
#'   interactions (such as \code{text}) are plotted. This argument will
#'   eventually be expanded to allow one to plot combinations of edge and
#'   centrality plots.
#' @param directedDiag See corresponding argument in the \code{\link{bootNet}}.
#'   function.
#' @param ... Additional arguments.
#'
#' @return A coefficient plot, difference plot, or correlation-stability plot.
#'   When \code{plot %in% c('none', FALSE)}, the table used to construct the
#'   relevant plot will be returned as output instead.
#' @export
#'
#' @seealso \code{\link{bootNet}, \link{resample}}
#'
#' @examples
#' \donttest{
#' boot1 <- bootNet(ggmDat, caseDrop = TRUE)
#'
#' plot(boot1)
#' plotBoot(boot1) # This functions the same as the command above
#'
#' boot2 <- bootNet(ggmDat)
#'
#' plot(boot2)
#' plot(boot2, difference = TRUE)
#' }
plotBoot <- function(x, type = 'edges', net = 'temporal', plot = 'all', cor = .7,
                     order = 'mean', ci = .95, pairwise = TRUE, interactions = TRUE,
                     labels = NULL, title = NULL, cis = 'quantile', true = NULL,
                     errbars = FALSE, vline = FALSE, threshold = FALSE,
                     difference = FALSE, color = FALSE, text = FALSE,
                     textPos = 'value', multi = NULL, directedDiag = FALSE, ...){
  nonzero <- is.character(difference)
  difference <- ifelse(is.character(difference), TRUE, difference)
  if(!is.null(multi)){
    call <- replace(as.list(match.call())[-1], 'multi', NULL)
    if(isTRUE(multi)){multi <- list(plot = c('pairwise', 'interactions'))}
    n <- names(multi)
    p1 <- invisible(do.call(plotBoot, replace(call, n, setNames(list(multi[[1]][1]), n))))
    legend <- switch(2 - ('legend' %in% names(call)), eval(call$legend), g_legend(p1))
    p2 <- invisible(do.call(plotBoot, replace(call, n, setNames(list(multi[[1]][2]), n))))
    return(gridExtra::grid.arrange(legend, gridExtra::arrangeGrob(
      p1 + theme(legend.position = 'none'),
      p2 + theme(legend.position = 'none'), nrow = 1),
      nrow = 2, heights = c(1, 10))
    )
  }
  args <- tryCatch({list(...)}, error = function(e){list()})
  fit0 <- switch(2 - ('fit0' %in% names(x)), x$fit0, list())
  if(isTRUE(attr(x, 'resample'))){
    if(!'data' %in% names(x)){return(plotCoefs(fit = x, title = title, ...))}
    if(length(fit0) == 0){
      fit0 <- switch(2 - ('fit0' %in% names(args)), args$fit0, TRUE)
    }
    x <- do.call(bootNet, append(list(data = x, fit0 = fit0, directedDiag = directedDiag),
                                 replace(args, 'fit0', NULL)))
  }
  if(!identical(threshold, FALSE) & 'bootFits' %in% names(x)){
    dat <- paste0('x$fit0$', ifelse('temporal' %in% names(x), 'SURnet$data', 'data'))
    if(isTRUE(attr(x$bootFits, 'resample'))){attributes(x$bootFits)$resample <- NULL}
    x <- bootNet(data = eval(parse(text = dat)), fits = x$bootFits,
                 threshold = threshold, fit0 = x$fit0,
                 directedDiag = directedDiag)
  }
  runonce <- TRUE
  plots <- c('none', 'pairwise', 'interactions', 'all', 'both')
  nets <- c('ggm', 'temporal', 'contemporaneous', 'between')
  types <- c('edges', 'strength', 'outstrength', 'instrength', 'ei', 'outei', 'inei')
  if(length(type) > 1){
    type <- match.arg(tolower(type), types, several.ok = TRUE)[1:2]
    stopifnot(all(type %in% c('outstrength', 'instrength')) | all(type %in% c('outei', 'inei')))
    runonce <- type[2]
    type <- type[1]
  }
  nty <- list(type, net, plot)
  plot <- which(!sapply(nty, is.character))
  plot <- ifelse(length(plot) == 0, unlist(nty)[nty %in% plots], nty[[plot]])
  nty <- unlist(nty)[-match(plot, unlist(nty))]
  nty <- match.arg(tolower(nty), c(nets, types), several.ok = TRUE)
  type <- nty[nty %in% types][1]
  net <- nty[nty %in% nets][1]
  net <- switch(net, between = 'ggm', 'NA' = 'temporal', net)
  type <- switch(type, outstrength = 'outStrength', instrength = 'inStrength',
                 outei = 'outEI', inei = 'inEI', ei = 'EI', 'NA' = 'edges', type)
  cis <- match.arg(tolower(cis), c('quantile', 'se'))
  #invisible(suppressMessages(require(ggplot2)))
  if(isTRUE(attr(x, 'mlGVAR')) & 'varMods' %in% names(x)){
    mlnet <- switch(net, ggm = 'between|means', 'fixed')
    if(is.null(title) & startsWith(mlnet, 'b')){title <- 'Between-subjects network\t'}
    fits <- x$varMods[[grep(mlnet, names(x$varMods))]]
    fit0 <- switch(2 - ('fit0' %in% names(args)), args$fit0, x[[grep(mlnet, names(x))]])
    data <- x$netData[[grep(mlnet, names(x$netData))]]
    args <- append(list(data = data, fits = fits, fit0 = fit0), replace(args, 'fit0', NULL))
    x <- do.call(bootNet, args)
  }
  lags <- any(c('temporal', 'contemporaneous') %in% names(x))
  if(lags){
    net <- switch(net, ggm = 'temporal', net)
    x <- x[[net]]
    if(net == 'contemporaneous'){boots <- x$boots}
    title <- switch(2 - is.null(title), paste(capitalize(net), 'Network\t'), title)
    if(identical(title, FALSE)){title <- NULL}
    if(type == 'strength' & net == 'temporal'){
      type <- 'outStrength'
    } else if(!type %in% c('edges', 'strength', 'EI') & net == 'contemporaneous'){
      type <- 'edges'
    }
    if(type == 'EI' & net == 'temporal'){
      type <- 'outEI'
    }
  } else {
    if(!type %in% c('edges', 'strength', 'EI')){type <- 'edges'}
    net <- 'ggm'
    boots <- x$boots
  }
  lags2 <- lags & (net == 'temporal')
  if(is.logical(plot)){plot <- ifelse(plot, 'all', 'none')}
  if(is.numeric(plot)){plot <- c('none', 'pairwise', 'interactions', 'all')[plot + 1]}
  pp <- match.arg(tolower(plot), c('all', 'both', 'pairwise', 'interactions', 'none'))
  pp <- switch(pp, both = 'all', none = FALSE, pp)
  if(pp == 'all' & (!pairwise | !interactions)){pp <- ifelse(!pairwise, 'interactions', 'pairwise')}
  if(pp == 'all' & difference | lags & !lags2 & difference){pp <- 'pairwise'}
  if(pp == 'interactions'){interactions <- TRUE; pairwise <- FALSE}
  if(pp == 'pairwise'){interactions <- FALSE; pairwise <- TRUE}
  type0 <- ifelse(type %in% c('outStrength', 'inStrength'), paste0('strength$', type), type)
  type0 <- ifelse(type %in% c('outEI', 'inEI'), paste0('EI$', type), type0)
  if(!"interactions" %in% names(x)){interactions <- FALSE}
  obj1 <- list()
  if(interactions | !pairwise){
    obj1 <- x$interactions
    if(lags){boots <- obj1$boots}
    p1 <- dat <- eval(parse(text = paste0('x$interactions$', type0)))
    if(pairwise){
      text <- FALSE
      obj1 <- x$pairwise
      obj2 <- x$interactions
      if('diffs' %in% names(obj1) & 'diffs' %in% names(obj2)){
        obj1$diffs <- setNames(lapply(seq_along(obj1$diffs), function(z){
          rbind(obj1$diffs[[z]], obj2$diffs[[z]])
        }), names(obj1$diffs))
      }
      p2 <- eval(parse(text = paste0('x$pairwise$', type0)))
      dat <- rbind.data.frame(p2, p1)
      dat$type <- rep(c("Pairwise", "Interactions"), each = nrow(dat)/2)
    } else {
      dat <- data.frame(dat)
      dat$type <- rep("Interactions", nrow(dat))
    }
  } else {
    boots <- x$boots
    if("pairwise" %in% names(x)){
      obj1 <- x$pairwise
      if(lags){boots <- obj1$boots}
      type0 <- paste0('pairwise$', type0)
    } else {
      obj1 <- x # THIS IS A NEW ADDITION
    }
    dat <- eval(parse(text = paste0('x$', type0)))
    dat$type <- rep("Pairwise", nrow(dat))
  }
  if('diffs' %in% names(obj1)){
    obj1 <- obj1$diffs
  } else if(lags & difference){
    obj1 <- x$diffs
  }
  if(!identical(text, FALSE) & type == 'edges'){ # had !lags
    boots <- boots[[switch(2 - interactions, ifelse(
      !lags, 'boot_int_edges', 'boot_edges'), 'boot_edges')]]
  }
  if(grepl('EI', type)){type <- gsub('I', 'Influence', gsub('E', 'Expected', type))}
  if(isTRUE(attr(x, "caseDrop"))){
    dat <- dat2 <- dat[order(dat$subN), ]
    ci <- paste0("q", c((1 - ci)/2, 1 - (1 - ci)/2) * 100)
    if(length(unique(dat$type)) == 1 & FALSE){ # MAY DELETE THIS WHOLE PART
      css <- dat[dat[, ci[1]] > cor, 'drop']
      css <- attr(dat2, 'CS') <- ifelse(length(css) == 0, 0, max(css))
      if(pp != FALSE){cat(paste0('CS: ', css, ' (cor = ', cor, ')'))}
    } else if(FALSE){ # MAY DELETE
      cssPair <- dat[dat$type == 'Pairwise', 'drop'][dat[dat$type == 'Pairwise', ci[1]] > cor]
      cssPair <- attr(dat2, 'CS_Pair') <- ifelse(length(cssPair) == 0, 0, max(cssPair))
      cssInt <- dat[dat$type == 'Interactions', 'drop'][dat[dat$type == 'Interactions', ci[1]] > cor]
      cssInt <- attr(dat2, 'CS_Int') <- ifelse(length(cssInt) == 0, 0, max(cssInt))
      if(pp != FALSE){cat(paste0('CS_Pair: ', cssPair, ' (cor = ', cor, ')'), '\n')}
      if(pp != FALSE){cat(paste0('CS_Int: ', cssInt, ' (cor = ', cor, ')'))}
    }
    legLab <- capitalize(type)
    N <- dat$N[1]
    p <- ggplot(dat, aes(x = subN, y = mean, group = type, colour = type, fill = type))
    p <- p + geom_line(lwd = 1) + geom_point()
    p <- p + geom_ribbon(colour = NA, alpha = .1, aes_string(ymin = ci[1], ymax = ci[2]))
    p <- p + scale_x_reverse(breaks = seq(.9, .1, by = -.1) * N,
                             labels = paste0(seq(90, 10, by = -10), "%"))
    p <- p + ylim(-1, 1) + xlab("Sampled cases") + ylab("Average correlation with original sample")
    p <- p + labs(fill = legLab, colour = legLab, group = legLab) + theme_bw()
    p <- p + geom_hline(yintercept = 0) + geom_hline(yintercept = cor, linetype = 2, alpha = .3)
    if(!is.null(title)){p <- p + ggtitle(title)}
  } else if(difference){
    if(order %in% c(6, 'true')){order <- TRUE}
    if(is.logical(order)){order <- ifelse(order, 'mean', 'id')}
    if(is.numeric(order)){order <- switch(order, 'id', 'sample', 'mean', 'v1', 'v2', 'true')}
    order <- match.arg(order, c('id', 'sample', 'mean', 'v1', 'v2', 'true'))
    if(!order %in% c('id', 'sample', 'mean')){order <- 'id'}
    order <- switch(order, mean = 'boot_mean', id = ifelse(type == 'edges', 'edge', 'node'), order)
    dat0 <- dat
    if(all(startsWith(names(obj1), 'int_'))){names(obj1) <- gsub('^int_', '', names(obj1))}
    if(startsWith(type, 'out') | startsWith(type, 'in')){
      tt <- ifelse(gsub('out|in', '', type) == 'Strength', 's', 'E')
      dat <- obj1[[which(sapply(names(obj1), startsWith, tt))]]
      dat <- dat[[which(sapply(c('^out', '^in'), grepl, type))]]
    } else {
      dat <- obj1[[which(sapply(names(obj1), startsWith, strsplit(type, '')[[1]][1]))]]
    }
    dat <- rbind(dat, cbind(setNames(dat[, 2:1], names(dat)[1:2]), dat[, -(1:2)]))
    dat$sig <- ifelse(dat$sig, 'sig', 'nonsig')
    d2 <- setNames(data.frame(unique(dat[, 1]), unique(dat[, 1]), 0, 0, 'same'), names(dat))
    dat <- dat2 <- rbind(dat, d2)
    colorValues <- c(same = 'white', nonsig = 'lightgray', sig = 'black')
    colnames(dat)[1:2] <- paste0('id', 1:2)
    id0 <- ifelse(type == 'edges', 'edge', 'node')
    dat$id1 <- factor(dat$id1, levels = dat0[[id0]][order(dat0[[order]])])
    dat$id2 <- factor(dat$id2, levels = dat0[[id0]][order(dat0[[order]])])
    dat$type <- factor(paste0(ifelse(pp == 'pairwise', 'Pairwise (', 'Interactions ('), capitalize(type), ')'))
    if(color & type %in% c('edges', 'pairwise', 'interactions')){
      # Adding threshold = TRUE has some issues when type = 'interactions'
      if(net == 'ggm'){net <- 'temporal'}
      nn <- switch(2 - pairwise, net(fit0, n = net), netInts(fit0))
      if(lags2 & !pairwise){nn <- t(nn)}
      nn <- rownames(nn)
      obj0 <- switch(2 - pairwise, fit0, t(netInts(fit0, threshold = threshold, avg = !lags2)))
      graph <- plotNet(obj0, which.net = net, threshold = threshold, plot = FALSE)
      edgelist <- data.frame(from = nn[graph$Edgelist$from], to = nn[graph$Edgelist$to],
                             col = graph$graphAttributes$Edges$color, stringsAsFactors = FALSE)
      if(all(grepl(':', edgelist$from))){
        n1 <- gsub(':.*', '', edgelist$from)
        n2 <- gsub(':.*', '', edgelist$to)
        mm <- unique(gsub('.*:', '', edgelist$from))
        edgelist$id <- paste0(paste0(n1, ifelse(lags2, '-->', '--'), n2), '|', mm)
      } else {
        if(lags){
          edgelist$from <- gsub('[.]y$', '', edgelist$from)
          edgelist$to <- gsub('[.]y$', '', edgelist$to)
        }
        edgelist$id <- paste0(edgelist$from, ifelse(lags2, '-->', '--'), edgelist$to)
      }
      if(any(is.na(edgelist$col))){edgelist$col[is.na(edgelist$col)] <- 'white'}
      dat$sig <- ifelse(dat$sig == 'same', as.character(dat$id1), dat$sig)
      colorValues <- c(setNames(edgelist$col, edgelist$id), colorValues)
    }
    if(is.null(labels)){labels <- ifelse(nrow(dat) <= 50, TRUE, ifelse(type == 'edges', 'numbers', FALSE))}
    if(is.character(labels) & length(labels) == 1){
      k <- unique(dat0$node1)
      kk <- unique(dat0$edge)
      if(!pairwise & interactions){
        mname <- unique(gsub('^.*.:', '', k))
        k <- gsub(paste0(':', mname), '', k)
        kk <- gsub(paste0('[|]', mname), '', kk)
      }
      for(i in seq_along(k)){kk <- gsub(k[i], i, kk)}
      if(!pairwise & interactions){kk <- paste0(kk, '|M')}
      names(kk) <- unique(dat0$edge)
      kk <- kk[match(levels(dat$id1), names(kk))]
      levels(dat$id1) <- levels(dat$id2) <- unname(kk)
      k1 <- dat0$edge; k2 <- dat0$node1; k3 <- dat0$node2
      if(!pairwise & interactions){
        k1 <- gsub(paste0('[|]', mname), '', k1)
        k2 <- gsub(paste0(':', mname), '', k2)
        k3 <- gsub(paste0(':', mname), '', k3)
      }
      for(i in seq_along(k)){
        k1 <- gsub(k[i], i, k1)
        k2 <- gsub(k[i], i, k2)
        k3 <- gsub(k[i], i, k3)
      }
      if(!pairwise & interactions){
        k1 <- paste0(k1, '|M')
        k2 <- paste0(k2, ':M')
        k3 <- paste0(k3, ':M')
      }
      dat0$edge <- k1
      dat0$node1 <- k2
      dat0$node2 <- k3
    }
    if(nonzero & type == 'edges' & any(dat0$sample == 0)){
      zeros <- dat0$edge[dat0$sample == 0]
      dat <- subset(dat, !id1 %in% zeros & !id2 %in% zeros)
      dat$id1 <- factor(dat$id1); dat$id2 <- factor(dat$id2)
    }
    p <- ggplot(dat, aes(x = id1, y = id2, fill = sig)) +
      geom_tile(colour = 'white') + xlab('') + ylab('') +
      scale_fill_manual(values = colorValues) + theme(legend.position = 'none') + facet_grid(~ type)
    p <- p + theme_grey(base_size = 9) +
      theme(legend.position = 'none', axis.text.x = element_text(
        size = 7.2, angle = 270, hjust = 0, colour = 'grey50'))
    if(identical(labels, FALSE)){p <- p + theme(axis.text.y = element_blank(), axis.text.x = element_blank())}
    if(text){
      n3 <- as.character(dat[!dat$sig %in% c('sig', 'nonsig'), 'id1'])
      dat3 <- data.frame(id1 = as.character(dat[!dat$sig %in% c('sig', 'nonsig'), 'id1']), value = NA)
      dat3[match(dat0[[id0]], dat3$id1), 'value'] <- format(signif(dat0$sample, 2), scientific = FALSE)
      dat3$id1 <- dat3$id2 <- factor(dat3$id1)
      dat3$sig <- dat[!dat$sig %in% c('sig', 'nonsig'), 'sig']
      p <- p + geom_text(data = dat3, aes(label = value))
    }
  } else {
    nets <- c('mean', 'sample')
    if(!is.null(true)){
      getType <- switch(type, edges = switch(net, temporal = c, function(x){x[lower.tri(x)]}),
                        strength = , outStrength = function(x){diag(x) <- 0; colSums(abs(x))},
                        inStrength = function(x){diag(x) <- 0; rowSums(abs(x))}, ExpectedInfluence = ,
                        outExpectedInfluence = function(x){diag(x) <- 0; colSums(x)},
                        inExpectedInfluence = function(x){diag(x) <- 0; rowSums(x)})
      dat$true <- numeric(nrow(dat))
      dat[dat$type == 'Pairwise', 'true'] <- getType(net(true, n = switch(net, ggm = 'between', net)))
      if('Interactions' %in% dat$type){dat[dat$type == 'Interactions', 'true'] <- getType(netInts(true))}
      nets <- c(nets, 'true')
    }
    if(is.null(true) & order %in% c('true', 6)){order <- TRUE}
    if(is.logical(order)){order <- ifelse(order, 'mean', 'id')}
    if(is.numeric(order)){order <- switch(order, 'id', 'sample', 'mean', 'v1', 'v2', 'true')}
    order <- match.arg(order, c('id', 'sample', 'mean', 'v1', 'v2', 'true'))
    order <- which(order == c('id', 'sample', 'mean', 'v1', 'v2', 'true'))
    id <- ifelse(type == 'edges', 'edge', 'node')
    lci <- ifelse(cis == 'quantile', 'boot_lower', 'sample_lower')
    uci <- ifelse(cis == 'quantile', 'boot_upper', 'sample_upper')
    if(ci != .95 & cis == 'quantile' & type == 'edges'){
      dat$boot_lower <- unname(apply(x$boots$boot_edges, 1, quantile, probs = (1 - ci)/2))
      dat$boot_upper <- unname(apply(x$boots$boot_edges, 1, quantile, probs = 1 - (1 - ci)/2))
    }
    v2 <- all(unique(dat$type) == c('Pairwise', 'Interactions'))
    if(!v2 & order > 3 & order != 6){order <- 1}
    if(order > 3 & order != 6){
      dat <- dat[order(dat$type, decreasing = TRUE), ]
      ord <- c('sample', 'boot_mean')[order - 3]
      o1 <- order(dat[dat$type == 'Pairwise', ord])
      o2 <- order(dat[dat$type == 'Interactions', ord]) + max(o1)
      dat <- dat[c(o1, o2), ]
      order <- 1
    } else if(order == 6){
      order <- 4
    }
    v2 <- TRUE
    type2 <- switch(2 - v2, paste0(rep(dat$type, length(nets)), ' (', capitalize(type), ')'), type)
    dat2 <- cbind.data.frame(
      id = rep(dat[[id]], length(nets)), value = unname(unlist(dat[, gsub('mean', 'boot_mean', nets)])),
      lower = rep(dat[[lci]], length(nets)), upper = rep(dat[[uci]], length(nets)),
      net = rep(nets, each = nrow(dat)), type = factor(capitalize(type2))
    )
    colnames(dat2)[1] <- id
    dat2[[id]] <- factor(dat2[[id]], levels = dat[[id]][switch(
      order, 1:nrow(dat), order(dat$sample), order(dat$boot_mean), order(dat$true))])
    colnames(dat2)[1] <- 'id'
    if(is.null(labels)){labels <- isTRUE(nrow(dat) <= 50)}
    scaleSize <- c('mean' = .7, 'sample' = 1, 'true' = .9)[seq_along(nets)]
    scaleAlpha <- c('mean' = .5, 'sample' = 1, 'true' = .5)[seq_along(nets)]
    legLabs <- c('Bootstrap mean', 'Sample', 'True network')[seq_along(nets)]
    legCols <- c('black', 'darkred', 'blue')[seq_along(nets)]
    if(!isTRUE(runonce) & !(pairwise & interactions)){
      call <- replace(as.list(match.call())[-1], c('type', 'pairwise', 'interactions', 'plot'),
                      list(type = runonce, pairwise = pairwise, interactions = interactions, plot = FALSE))
      dat2 <- rbind.data.frame(dat2, do.call(plotBoot, call))
    }
    p <- ggplot(dat2, aes(x = id, y = value, group = net, colour = net, fill = net))
    p <- p + geom_point(aes(alpha = net), show.legend = c(alpha = FALSE))
    if(!errbars){
      p <- p + geom_line(aes(size = net, alpha = net), show.legend = c(colour = FALSE, alpha = FALSE, size = FALSE))
      p <- p + geom_ribbon(aes(ymin = lower, ymax = upper, colour = NULL, fill = NULL), alpha = .1)
    } else {
      p <- p + geom_errorbar(aes(ymin = lower, ymax = upper, colour = NULL, fill = NULL), width = .2)
    }
    p <- p + facet_grid(~ type) + coord_flip() + theme_bw() + xlab('') + ylab('')
    p <- p + guides(colour = guide_legend(title = title), fill = 'none') + theme(legend.position = "top")
    p <- p + scale_size_manual(values = scaleSize) + scale_alpha_manual(values = scaleAlpha)
    p <- p + scale_color_manual('', values = legCols, labels = legLabs)
    if(!isTRUE(labels)){p <- p + theme(axis.text.y = element_blank())}
    if(!identical(vline, FALSE)){
      p <- p + geom_hline(yintercept = ifelse(!is.numeric(vline), 0, vline), linetype = 2, alpha = .35)
    }
    if(!identical(text, FALSE) & type == 'edges'){
      # Would be cool to make a warning if no thresholds are used
      if(isTRUE(text)){text <- 1.2}
      rnames <- rownames(boots)
      dat2$prop0 <- rep(NA, nrow(dat2))
      dat2 <- dat2[dat2$net == 'mean', ]
      dat2[match(rnames, as.character(dat2$id)), 'prop0'] <- rowMeans(data.frame((boots != 0) * 1))
      if(!identical(pp, FALSE) & textPos != 'value'){dat2$value <- textPos}
      p <- p + geom_label(aes(y = value, label = format(round(prop0, 2), nsmall = 2), fill = NULL),
                          data = dat2, cex = text * 2, label.padding = unit(0.1, 'lines'),
                          label.size = 0.1, alpha = .8, colour = 'black')
    }
  }
  if(pp != FALSE){return(p)} else {return(dat2)}
}

#' @rdname plotBoot
#' @export
plot.bootNet <- function(x, type = 'edges', net = 'temporal', plot = 'all', cor = .7,
                         order = 'mean', ci = .95, pairwise = TRUE, interactions = TRUE,
                         labels = NULL, title = NULL, cis = 'quantile', true = NULL,
                         errbars = FALSE, vline = FALSE, threshold = FALSE,
                         difference = FALSE, color = FALSE, text = FALSE,
                         textPos = 'value', multi = NULL, directedDiag = FALSE, ...){
  args <- as.list(match.call())[-1]
  do.call(plotBoot, args)
}

#' Plot the ECDF of p-values from resampling
#'
#' Plots the empirical cumulative distribution function of the p-values related
#' to iterated resampling via bootstrapping or multi-sample splitting.
#'
#' See Meinshausen, Meier, & Buhlmann (2009) for details.
#'
#' @param x Output from \code{\link{resample}}, given that \code{sampMethod =
#'   "bootstrap"} or \code{sampMethod = "split"}.
#' @param outcome Character string or numeric value (in terms of columns in the
#'   dataset) to indicate which outcome to plot the p-value distribution for.
#' @param predictor Character string or numeric value (in terms of columns in
#'   the dataset) to indicate which predictor to plot the p-value distribution
#'   for.
#' @param title If \code{TRUE}, then a default title will be given according to
#'   the outcome and predictor that are specified. If \code{FALSE}, then no
#'   title will be plotted. A custom title may also be supplied by the user.
#' @param alpha The false discovery rate. Defaults to .05
#'
#' @return Returns a plot based on the relationship between a particular outcome
#'   and predictor.
#' @export
#' @references Meinshausen, N., Meier, L., & Buhlmann, P. (2009). P-values for
#'   high-dimensional regression. Journal of the American Statistical
#'   Association. 104, 1671-1681.
#'
#' @seealso \code{\link{resample}}
#'
#' @examples
#' \donttest{
#' x <- resample(ggmDat, sampMethod = "bootstrap")
#' plot(x, what = 'pvals')
#' plot(x, 'pvals', outcome = 'V2', predictor = 'V1')
#' }
plotPvals <- function(x, outcome = 1, predictor = 1, title = TRUE, alpha = .05){
  stopifnot("adjCIs" %in% names(x))
  pvals <- lapply(x$samples$coefs, '[[', "P")
  if(is.character(outcome)){outcome <- which(names(pvals) == outcome)}
  if(is.character(predictor)){predictor <- which(colnames(pvals[[outcome]]) == predictor)}
  if(isTRUE(title)){
    title <- paste(colnames(pvals[[outcome]])[predictor], '-->', names(pvals)[outcome])
  } else if(identical(title, FALSE)){
    title <- ''
  }
  ps <- pvals[[outcome]][, predictor]
  n <- length(ps)
  gammas <- seq(ceiling(alpha * n)/n, 1 - 1/n, by = 1/n)
  h <- hist(ps, breaks = n, plot = FALSE)
  h$density <- h$density/sum(h$density)
  fdr_line <- function(p, gammas, alpha = .05){pmax(.05, (1 - (log(min(gammas)))/alpha) * p)}
  plot(h, freq = FALSE, main = title, ylim = c(0, 1), xlab = "Adjusted P-values")
  lines(sort(ps), fdr_line(sort(ps), gammas, alpha), lty = 2, lwd = 2)
  lines(sort(ps), seq(.05, 1, length = n), lwd = 2)
}

#' Plot stability selection paths for a given outcome
#'
#' Creates a plot to show the stability path for a particular variable in terms
#' of how frequently it was chosen in stability selection.
#'
#' See Meinshausen & Buhlmann (2010) for details on stability selection. Cannot
#' be used when the criterion for stability selection was set as
#' cross-validation.
#'
#' @param x Output of \code{\link{resample}} where \code{sampMethod =
#'   "stability"}.
#' @param outcome Character string or numeric value (in terms of columns in the
#'   dataset) to indicate which outcome to plot the stability path for.
#' @param s Character string or numeric value. This indicates which stability
#'   path to return a plot for. Either the first sample split \code{"split1"},
#'   the second sample split \code{"split2"}, or the path for simultaneous
#'   selection \code{"simult"}, which is the default.
#' @param thresh The selection threshold, which is represented as a horizontal
#'   red line on the plot. Defaults to .5
#' @param typeLegend Logical. If \code{FALSE}, linetype legend is removed. Only
#'   applies if there is a moderator in the model.
#'
#' @return Plot of the stability path associated with a given outcome.
#' @export
#' @references Meinshausen, N., & Buhlmann, P. (2010). Stability selection.
#'   Journal of the Royal Statistical Society: Series B (Statistical
#'   Methodology). 72, 417-423
#'
#' @seealso \code{\link{resample}}
#'
#' @examples
#' \donttest{
#' x <- resample(ggmDat, sampMethod = "stability")
#' plot(x, what = "stability")
#' plot(x, 'stability', outcome = 'V3')
#' }
plotStability <- function(x, outcome = 1, s = c('simult', 'split1', 'split2'),
                          thresh = .5, typeLegend = TRUE){
  if(x$call$criterion == "CV"){stop("Not possible when criterion == CV")}
  objcall <- x$call
  mnull <- is.null(objcall$moderators)
  x <- x$stability
  if(!'simult' %in% names(x)){
    stop('No stability paths available for object. Try re-running resample with a different seed.')
  }
  if(is.character(s)){s <- switch(match.arg(s), simult = 3, split1 = 1, split2 = 2)}
  if(is.character(outcome)){outcome <- which(names(x[[s]]) == outcome)}
  node <- names(x[[s]])[outcome]
  p <- length(x[[s]]) - ifelse(s == 3, 0, ifelse(!objcall$exogenous, 2, 1))
  stopifnot(outcome <= p + 1)
  x <- x[[s]][[outcome]]
  k <- ncol(x)
  n <- nrow(x)
  x <- data.frame(stack(x), lambda = rep(1:nrow(x), ncol(x)))
  lambda <- x$lambda # removing NOTE
  values <- x$values # removing NOTE
  ind <- x$ind # removing NOTE
  if(!mnull){
    x$Type <- factor(rep(c('Main Effect', 'Interaction'), c(p, k - p) * n),
                     levels = c('Main Effect', 'Interaction'))
  }
  g <- ggplot(x, aes(x = lambda, y = values, color = ind)) +
    xlab(expression(lambda)) + ylab('Selection Probability') + ggtitle(node) +
    labs(color = 'Predictor') + theme_classic()
  if(!mnull){
    g <- g + geom_line(aes(linetype = Type))
    if(!isTRUE(typeLegend)){g <- g + guides(linetype = 'none')}
  } else {
    g <- g + geom_line()
  }
  if(!is.null(thresh)){
    if(!identical(as.numeric(thresh), 0)){
      g <- g + geom_hline(yintercept = thresh, alpha = .3, lty = 4)
    }
  }
  return(g)
}

#' Plot model coefficients with confidence intervals
#'
#' Return a plot or dataframe showing the point estimates from each model, along
#' with confidence intervals based on the estimated standard errors.
#'
#' This is differentiated from the output of \code{\link{bootNet}} and
#' \code{\link{plotBoot}} in that the confidence intervals are computed directly
#' from model parameters rather than estimated from bootstrapping.
#'
#' @param fit Output from \code{\link{fitNetwork}}, \code{\link{bootNet}}, or
#'   \code{\link{resample}}. Can also be the \code{fixedNets} or
#'   \code{betweenNet} elements of the \code{\link{mlGVAR}} output.
#' @param true An adjacency matrix containing the true parameter values, if
#'   known. This can be used in conjunction with a simulated network, in that
#'   the user can supply the true network and plot those values against the
#'   estimated values.
#' @param alpha Alpha level that is used to compute confidence intervals.
#' @param plot Logical. If \code{FALSE}, a dataframe containing all of the
#'   confidence interval data will be returned.
#' @param col Character string. Color of the points associated with the
#'   \code{true} values.
#' @param flip Logical. If \code{FALSE}, the facets will be turned 90 degrees.
#' @param data Supply the original dataset if not already included in the
#'   \code{fit} object.
#' @param select Relevant to the \code{\link{resample}} output. Determines
#'   whether all variables should be plotted, or only those that were selected
#'   according to the resampling or variable selection procedure.
#' @param size Numeric. Size of the point estimates.
#' @param labels If logical, determines whether or not variable labels should be
#'   included. If a character vector, can be used to customize variable labels.
#' @param title Custom plot title.
#' @param vars Defaults to \code{"all"}. Determines which variables should be
#'   plotted.
#'
#' @return Plot displaying estimated model coefficients and confidence
#'   intervals.
#' @export
#'
#' @seealso \code{\link{fitNetwork}, \link{resample}, \link{getFitCIs},
#'   \link{plot.resample}, \link{plotNet}}
#'
#' @examples
#' \donttest{
#' x <- fitNetwork(ggmDat)
#' plot(x, which.net = 'coefs')
#' plotCoefs(x) # This is the same as the above command
#' }
plotCoefs <- function(fit, true = FALSE, alpha = .05, plot = TRUE, col = "blue",
                      flip = TRUE, data = NULL, select = TRUE, size = 1,
                      labels = TRUE, title = NULL, vars = 'all'){
  if(isTRUE(attr(fit, 'mlGVAR')) & 'varMods' %in% names(fit)){
    if(class(true) %in% c('logical', 'character', 'numeric')){
      if(is.numeric(true)){true <- as.logical(true - 1)}
      true <- ifelse(is.logical(true), c('fixed', 'between')[true + 1], match.arg(
        true, c('fixed', 'between', 'temporal', 'contemporaneous', 'ggm')))
    } else stop('Must provide resample output to include true network')
    true <- switch(true, ggm = , between = 'between', 'fixed')
    fit <- fit$varMods[[grep(true, names(fit$varMods))]]
    if(is.null(title) & !identical(title, FALSE)){
      title <- switch(true, fixed = 'Temporal Network', 'Between-subjects network')
    }
    true <- FALSE
  }
  if(is(fit, "systemfit")){
    x <- as.matrix(confint(fit, level = (1 - alpha)))
    x <- cbind(x, coef(fit))
    x <- data.frame(x)
    colnames(x) <- c("lower", "upper", "b")
    x <- x[!grepl("Intercept", rownames(x)),]
    x$predictor <- rownames(x)
    rownames(x) <- 1:nrow(x)
    x$Y <- NA
    for(i in 1:length(fit$eq)){
      x[grepl(paste0("eq", i), x$predictor), "Y"] <- as.character(fit$eq[[i]]$terms[[2]])
    }
    x$Y <- as.factor(x$Y)
    #x$predictor <- as.factor(qdapRegex::rm_between(x$predictor, "eq", "_"))
    x$predictor <- as.factor(sapply(strsplit(x$predictor, '_'), function(z) paste0(z[-1], collapse = '_')))
    x <- data.frame(Y = x$Y, predictor = x$predictor, lower = x$lower, b = x$b, upper = x$upper)
    if(!is.null(true)){
      dat <- tryCatch({eval(fit$call$data)}, error = function(e){
        if(!is.null(data)){data} else {stop("Need to supply data to index names of predictors")}
      })
      if(class(dat) == "list"){dat <- cbind.data.frame(dat$X, dat$Y)} # FULLFIX
      nY <- rep(levels(x$Y), each = (ncol(true) - 1))
      nP <- rep(colnames(dat)[!colnames(dat) %in% levels(x$Y)], length(levels(x$Y)))
      trueDat <- data.frame(nY, nP, B = as.vector(t(true[,-1])))
    }
    dat <- list()
    for(i in 1:length(fit$eq)){dat[[i]] <- x[grepl(paste0("^X", i, ".y$"), x$Y),]}
    dat <- do.call(rbind, dat)
  } else {
    if(all(c("mod0", "mod1se") %in% names(fit))){
      stop("Must index either 'mod0' or 'mod1se' to indicate which to use")
    } else if(any(grepl("CIs|freqs", names(fit)))){
      if(is.null(true) & !any(grepl("freqs", names(fit)))){
        true <- if("adjCIs" %in% names(fit)){fit$fits$fit1} else {fit$fit1}}
      fit <- ifelse("adjCIs" %in% names(fit), list(fit$adjCIs),
                    ifelse("freqs" %in% names(fit), list(fit$freqs),
                           list(fit$fitCIs)))[[1]]
    } else if(any(c('fitobj', 'SURfit') %in% names(fit))){
      fit <- getFitCIs(fit)
    }
    if(is.logical(true)){true <- NULL}
    ynames <- names(fit)
    x <- do.call(rbind, fit)
    if("select" %in% colnames(fit[[1]])){
      x <- cbind.data.frame(Y = rep(ynames, sapply(fit, nrow)), x)
    } else {
      x <- cbind.data.frame(Y = rep(ynames, each = length(levels(x$predictor))), x)
    }
    if(!is.null(true)){
      if(is.list(true)){
        true2 <- unlist(lapply(lapply(true$mods, '[[', "model"), function(z) z[-1, ]))
        if(length(as.vector(true2)) != nrow(x)){
          true3 <- names(true2)
          true4 <- rep(NA, nrow(x))
          xx <- unname(apply(x[, 1:2], 1, paste0, collapse = "."))
          true4[match(true3, xx)] <- true2
          true2 <- true4
        }
        trueDat <- data.frame(nY = x$Y, nP = as.character(x$predictor), B = as.vector(unname(true2)))
      } else if(!is.logical(true)){
        if(ncol(true) == (2 * nrow(true))){true <- true[, -1]}
        trueDat <- data.frame(nY = x$Y, nP = as.character(x$predictor), B = as.vector(t(true)))
      }
    }
    if(any(rowSums(is.na(x)) != 0)){x <- x[rowSums(is.na(x)) == 0, ]}
    x$Y <- factor(x$Y)
    x$predictor <- factor(x$predictor)
    dat <- x
  }
  rownames(dat) <- 1:nrow(dat)
  dat$ord <- paste(dat$Y, dat$predictor, sep = "_")
  if(!is.null(true)){
    trueDat$ord <- paste(trueDat$nY, trueDat$nP, sep = "_")
    dat$B <- trueDat[trueDat$ord %in% dat$ord, "B"]
  }
  if(select != FALSE & "select" %in% names(dat)){
    if(!isTRUE(select)){dat$select <- !(dat$lower <= 0 & dat$upper >= 0)}
    dat <- dat[dat$select, ]
  }
  if(!identical(vars, 'all')){
    vs <- gsub('[.]y$', '', as.character(unique(dat$Y)))
    if(is.numeric(vars)){vars <- vs[vars]}
    dat$Y <- as.character(dat$Y)
    dat <- subset(dat, grepl(paste0('^', vars, collapse = '|'), Y))
    dat$Y <- factor(dat$Y)
  }
  if(plot == TRUE){
    plotCI <- function(dat, xlabs, true = NULL, flip = TRUE,
                       size = 1, labels = TRUE, title = NULL){
      #invisible(suppressMessages(require(ggplot2)))
      if(is.logical(labels)){if(!labels){xlabs <- NULL}} else {xlabs <- labels}
      p <- ggplot(dat, aes(x = reorder(ord, b), y = b, group = Y)) +
        facet_wrap(Y ~., scales = "free") + geom_point(size = size, na.rm = TRUE) +
        geom_errorbar(aes(ymin = lower, ymax = upper), width = .2) +
        geom_hline(aes(yintercept = 0), linetype = "dashed") +
        scale_x_discrete(labels = xlabs) + xlab("Predictor") +
        ylab(expression(hat(beta))) + theme_bw()
      if(!is.null(true)){
        p <- p + geom_point(data = dat, aes(x = reorder(ord, B), y = B, group = Y),
                            col = col, alpha = .85, shape = 8, na.rm = TRUE, size = size)
      }
      if(!is.null(title)){p <- p + ggtitle(title)}
      if(flip == TRUE){return(p + coord_flip())} else {return(p)}
    }
    plotCI(dat, xlabs = function(x){sub("[^_]*_", "", x)}, true = true,
           flip = flip, size = size, labels = labels, title = title)
  } else {
    return(dat)
  }
}

#' Conditional effects plot
#'
#' Creates a plot of the relationships between two variables at different levels
#' of the moderator. Only works for relationships that include an interaction.
#'
#' @param out Output from \code{\link{fitNetwork}} or \code{\link{resample}}.
#'   Can also provide the \code{fixedNets} or \code{betweenNet} element of the
#'   \code{\link{mlGVAR}} output.
#' @param to Outcome variable, specified with character string or numeric value.
#' @param from Predictor variable, specified with character string or numeric
#'   value.
#' @param swap Logical. Serves to switch the arguments for \code{to} and
#'   \code{from}.
#' @param avg Logical. If \code{TRUE} then the average relationship between the
#'   two variables is displayed. Only works for GGMs.
#' @param compare Two values can be supplied to indicate levels of the moderator
#'   to be compared.
#' @param hist Logical. Determines whether to show a histogram of the data
#'   distribution at the bottom of the plot.
#' @param xlab Character string for labeling the x-axis.
#' @param mods This argument will be removed. Model output is automatically
#'   detected based on \code{fit} argument.
#' @param nsims Number of iterations to simulate the posterior distribution.
#' @param xn Numeric value to indicate how many values of the moderator should
#'   be evaluated.
#' @param getCIs Logical. Only applies when \code{avg = TRUE}. If \code{getCIs =
#'   TRUE}, then the confidence intervals for the average difference between the
#'   maximum and minimum of the moderator will be returned.
#' @param discrete Logical. Determines whether to treat the moderator as a
#'   discrete or continuous variable.
#' @param ylab Character string for labeling the y-axis.
#' @param main Character string for labeling the title of the plot.
#' @param midline Logical. Only applies when \code{discrete = TRUE}. Shows a
#'   line at the average level of the outcome.
#'
#' @return A plot of the conditional effects of one variable on another given
#'   different levels of the moderator.
#' @export
#'
#' @seealso \code{\link{fitNetwork}, \link{resample}}
#'
#' @examples
#' \donttest{
#' fit <- fitNetwork(ggmDat, 'M')
#' condPlot(fit, to = 'V5', from = 'V4')
#' condPlot(fit, to = 2, from = 3, avg = TRUE)
#' }
condPlot <- function(out, to, from, swap = FALSE, avg = FALSE, compare = NULL,
                     hist = FALSE, xlab = NULL, mods = NULL, nsims = 500,
                     xn = NULL, getCIs = FALSE, discrete = FALSE,
                     ylab = NULL, main = NULL, midline = TRUE){
  #suppressMessages(invisible(require(ggplot2)))
  if(isTRUE(discrete)){compare <- 0:1}
  if(is(out, 'resample')){out <- out$fit0}
  if("adjMat" %in% names(out)){out <- out$mods0}
  if(any(c("models", "SURnet") %in% names(out))){
    out <- condEffects(out, xn = xn, x = compare)}
  if(is.null(xlab)){xlab <- attributes(out)$moderator}
  if(length(to) == 2){from <- to[2]; to <- to[1]}
  if(class(to) == "numeric"){to <- names(out)[to]}
  if(class(from) == "numeric"){from <- names(out)[from]}
  if(swap){
    toX <- to
    to <- from
    from <- toX
  }
  data <- out[[to]][[from]]
  stopifnot(!is.null(data))
  stopifnot(nrow(data) > 1)
  if(avg & !"SURnet" %in% names(attributes(out))){
    stopifnot(!is.null(out[[from]][[to]]))
    newdat <- data.frame(matrix(NA, ncol = 5, nrow = nrow(data)))
    newdat[, 1] <- data$x
    newdat[, 2] <- rowMeans(cbind(data$y, out[[from]][[to]]$y))
    newdat[, 3] <- rowMeans(cbind(data$se, out[[from]][[to]]$se))
    newdat[, 4] <- rowMeans(cbind(data$lower, out[[from]][[to]]$lower))
    newdat[, 5] <- rowMeans(cbind(data$upper, out[[from]][[to]]$upper))
    colnames(newdat) <- c("x", "y", "se", "lower", "upper")
    data <- newdat
  }
  to2 <- capitalize(to)
  fr2 <- capitalize(from)
  xlab <- capitalize(xlab)
  if(!avg){
    if(is.null(ylab)){ylab <- bquote(.(fr2)%->%.(to2))}
    #if(is.null(ylab)){ylab <- paste0(fr2, " ---> ", to2)}
    if(is.null(main)){mainlab <- paste0("Estimated Effect of ", fr2, " on ", to2, " across levels of ", xlab)}
  } else if(!"SURnet" %in% names(attributes(out))){
    if(is.null(ylab)){ylab <- paste0(fr2, " ~ ", to2)}
    if(is.null(main)){mainlab <- paste0("Mean Relation between ", fr2, " and ", to2, " across levels of ", xlab)}
  }
  if(!is.null(main)){mainlab <- main}
  alpha <- attributes(out)$alpha
  if("mods" %in% names(attributes(out))){mods <- attributes(out)$mods}
  if(!is.null(mods) & !"SURnet" %in% names(attributes(out))){
    ci_diff <- margCIs(mods = mods, alpha = alpha, nsims = nsims, compare = compare)
    if(avg){ci_diff2 <- ci_diff[[from]][to, ]}
    ci_diff <- ci_diff[[to]][from, ]
  } else {
    ci_diff <- NULL
  }
  if(nrow(data) != 2 & discrete == FALSE){
    if(!hist){
      pp <- ggplot(data = data, aes(x = x, y = y)) +
        geom_line(color = "red") +
        geom_ribbon(alpha = .2, aes(ymin = lower, ymax = upper))
      if(min(data$lower) < 0 & max(data$upper) > 0){pp <- pp + geom_hline(yintercept = 0, linetype = 2)}
      pp <- pp + xlab(xlab) + ylab(ylab) + ggtitle(mainlab) + theme_bw()
    } else {
      var2_dt <- ifelse(
        "SURnet" %in% names(attributes(out)),
        list(mods), mods$dat[, ncol(mods$dat)])[[1]]
      yrange <- c(data$upper, data$lower)
      maxdiff <- (max(yrange) - min(yrange))
      breaks_var2 <- nrow(data)
      if(hist == "unique"){breaks_var2 <- length(unique(var2_dt))}
      hist.out <- hist(var2_dt, breaks = breaks_var2, plot = F)
      n.hist <- length(hist.out$mids)
      dist <- hist.out$mids[2] - hist.out$mids[1]
      hist.max <- max(hist.out$counts)
      histX <- data.frame(ymin = rep(min(yrange) - maxdiff/5, n.hist),
                          ymax = hist.out$counts/hist.max * maxdiff/5 + min(yrange) - maxdiff/5,
                          xmin = hist.out$mids - dist/2, xmax = hist.out$mids + dist/2)
      pp <- ggplot()
      pp <- pp + geom_rect(data = histX, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                           colour = "gray50", alpha = 0, size = .5)
      pp <- pp + geom_line(data = data, aes(x = x, y = y), color = "red")
      pp <- pp + geom_ribbon(data = data, aes(x = x, ymin = lower, ymax = upper), alpha = .2)
      if(min(data$lower) < 0 & max(data$upper) > 0){pp <- pp + geom_hline(yintercept = 0, linetype = 2)}
      pp <- pp + xlab(xlab) + ylab(ylab) + ggtitle(mainlab) + theme_bw()
    }
  } else {
    if(isTRUE(discrete)){data$x <- factor(data$x, levels = 0:1)}
    pp <- ggplot(data = data, aes(x = x, y = y)) +
      geom_point() + geom_errorbar(aes(ymin = lower, ymax = upper), width = .05)
    if(midline){pp <- pp + geom_hline(yintercept = mean(data$y), linetype = 3, alpha = .6, color = "red")}
    if(min(data$lower) < 0 & max(data$upper) > 0){pp <- pp + geom_hline(yintercept = 0, linetype = 2, alpha = .6)}
    pp <- pp + scale_x_discrete(limits = data$x)
    pp <- pp + xlab(xlab) + ylab(ylab) + ggtitle(mainlab) + theme_bw()
  }
  if(!is.null(ci_diff)){
    if(!is.null(compare)){
      compare <- round(compare, 2)
      citxt <- paste0(100 - (alpha * 100), "% CI(", compare[2], " - ", compare[1], "): [")
    } else {
      citxt <- paste0(100 - (alpha * 100), "% CI(Max - Min): [")
    }
    if(avg){
      a1 <- paste0(fr2, " ---> ", to2, ":  ")
      b1 <- paste0(to2, " ---> ", fr2, ":  ")
      a2 <- paste0(round(ci_diff[2], 3), ", ", round(ci_diff[3], 3), "]")
      b2 <- paste0(round(ci_diff2[2], 3), ", ", round(ci_diff2[3], 3), "]")
      pp <- pp + labs(subtitle = paste0(a1, citxt, a2, "\n", b1, citxt, b2))
      if(getCIs){
        ci_diff0 <- list(ci_diff, ci_diff2)
        names(ci_diff0) <- c(paste0(fr2, ":", to2), paste0(to2, ":", fr2))
        return(ci_diff0)
      }
    } else {
      pp <- pp + labs(subtitle = paste0(citxt, round(ci_diff[2], 3), ", ", round(ci_diff[3], 3), "]"))
    }
  }
  pp
}

#' Plot confidence intervals for interaction terms
#'
#' Allows one to plot the confidence intervals associated with interaction
#' terms. Provides an easy way to look at whether there are any significant
#' interactions, and if so which interactions are important.
#'
#' The default setting \code{y = "all"} shows all interaction terms associated
#' with the model. But the user can also home-in on specific variables to see
#' what interactions might be relevant. When \code{y = "all"}, the axis labels
#' should be explained. These follow the format of \code{predictor:outcome}. The
#' title reflects the name of the moderator variable. For instance, if a
#' variable named \code{"M"} moderates the relationship between \code{"X"} and
#' \code{"Y"}, where \code{"X"} predicts \code{"Y"}, the title of the plot will
#' list the variable \code{"M"} as the moderator, and the label (shown on the
#' y-axis), will read \code{"X:Y"}. When \code{y != "all"} (that is, a specific
#' value for \code{y} is provided), then the title will still reflect the
#' moderator, but the labels will simply show which predictor interacts with
#' that moderator to predict the outcome.
#'
#' @param out GGM moderated network output from \code{\link{fitNetwork}}, or
#'   output from a moderated between-subjects network fit with
#'   \code{\link{mlGVAR}} (e.g., when \code{bm = TRUE}).
#' @param y Character string. The name of the outcome variable for which to
#'   create the plot. If \code{y = "all"}, then all interaction terms associated
#'   with all outcomes will be plotted.
#' @param nsims The number of simulations to estimate the posterior distribution
#'   of the difference between high and low levels of the confidence interval.
#' @param alpha Alpha level that is used to compute confidence intervals.
#'
#' @return A plot showing the spread of different interactions.
#' @export
#'
#' @seealso \code{\link{fitNetwork}, \link{plotNet}, \link{mlGVAR}}
#'
#' @examples
#' fit <- fitNetwork(ggmDat, 'M')
#' plot(fit, 'ints', y = 'all')
intsPlot <- function(out, y = 'all', nsims = 500, alpha = .05){
  if(is(out, 'SURnet')){stop('Cannot use this function with temporal networks')}
  if(is(out, 'mlGVAR')){out <- out$betweenNet}
  if("adjMat" %in% names(out)){out <- out$mods0}
  if("models" %in% names(out)){out <- margCIs(out, nsims = nsims, alpha = alpha)}
  moderator <- attributes(out)$moderator
  if(class(y) == "character"){
    if(y == "all"){
      y0 <- as.vector(sapply(seq_along(out), function(z)
        paste0(names(out)[z], ":", rownames(out[[z]]))))
      dd <- suppressWarnings(data.frame(do.call(rbind, out)))
      if(any(grepl(":$", y0))){y0 <- y0[-grep(":$", y0)]}
      if(class(y0) == "list"){y0 <- unlist(y0)}
      rownames(dd) <- y0
    } else if(y == "sig"){
      out2 <- lapply(seq_along(out), function(z){
        z2 <- data.frame(out[[z]], sig = as.numeric(!(out[[z]][, 2] < 0 & out[[z]][, 3] > 0)),
                         check.names = FALSE)
        return(z2[z2$sig == 1, -4])
      })
      names(out2) <- names(out)
      if(any(sapply(out2, length) == 0)){
        out2 <- out2[-which(sapply(out2, length) == 0)]
      }
      y0 <- unlist(as.vector(sapply(seq_along(out2), function(z)
        paste0(names(out2)[z], ":", rownames(out2[[z]])))))
      dd <- suppressWarnings(data.frame(do.call(rbind, out2)))
      if(any(grepl(":$", y0))){y0 <- y0[-grep(":$", y0)]}
      if(class(y0) == "list"){y0 <- unlist(y0)}
      rownames(dd) <- y0
    } else {
      y <- which(names(out) == y)
      dd <- data.frame(out[[y]])
      Y <- capitalize(names(out)[y])
    }
  } else {
    dd <- data.frame(out[[y]])
    Y <- capitalize(names(out)[y])
  }
  M <- capitalize(moderator)
  if(y %in% c("sig", "all")){
    Y <- c("Significant interaction effects", "All interactions")[which(c("sig", "all") %in% y)]
    main <- paste0(Y, " across levels of ", M)
  } else {
    main <- paste0("Effects on ", Y, " across levels of ", M)
  }
  colnames(dd)[2:3] <- c("lower", "upper")
  dd$pred <- paste0("_", rownames(dd))
  p <- ggplot(data = dd, aes(x = reorder(pred, b), y = b)) +
    geom_point() + geom_errorbar(aes(ymin = lower, ymax = upper), width = .2) +
    geom_hline(yintercept = 0, linetype = 2) +
    scale_x_discrete(labels = function(x){sub("[^_]*_", "", x)}) +
    ggtitle(label = main) + xlab(Y) + ylab(expression(hat(beta))) +
    theme_bw() + coord_flip()
  p
}

#' Plot results of power simulations
#'
#' Plots the output from the \code{\link{mnetPowerSim}} function.
#'
#' The options of what performance metrics to plot include: \itemize{
#' \item{Sensitivity} \item{Specificity} \item{Correlation} \item{MAE (Mean
#' Absolute Error)} \item{Precision} \item{Accuracy} \item{FDR (False Discovery
#' Rate)} }
#'
#' @param x \code{\link{mnetPowerSim}} output
#' @param by In development. Currently only supports \code{"type"} for creating
#'   different facets for Pairwise and Interaction effects. \code{"network"} for
#'   creating facets based on different networks (e.g., temporal,
#'   contemporaneous). \code{"p"} for creating facets based on the number of
#'   nodes in the network.
#' @param yvar The performance metrics to plot. Options include:
#'   \code{"sensitivity", "specificity", "correlation", "precision", "MAE",
#'   "FDR", "accuracy"}. The option \code{"default"} automatically sets this to
#'   sensitivity, specificity, and correlation.
#' @param yadd Specify additional performance metrics to plot. The final
#'   performance metrics that end up being plotted are simply: \code{c(yvar,
#'   yadd)}. Thus, this argument is only useful as a shortcut for keeping the
#'   default values of \code{yvar}, but adding more metrics to plot.
#' @param hline Numeric value between 0 and 1 for where to plot a horizontal
#'   line of interest. Can set to \code{FALSE} to remove line.
#' @param xlab Character string for the x-axis label.
#' @param title Character string for the title of the plot.
#' @param ... Additional arguments.
#'
#' @return Plots the results of a power simulation according to a variety of
#'   performance metrics.
#' @export
#'
#' @seealso \code{\link{mnetPowerSim}}
#'
#' @examples
#' \donttest{
#' x <- mnetPowerSim(niter = 10, N = c(100, 200))
#' summary(x)
#' plot(x)
#' }
plotPower <- function(x, by = 'type', yvar = 'default', yadd = NULL, hline = .8,
                      xlab = 'Number of cases', title = NULL, ...){
  args <- tryCatch({list(...)}, error = function(e){list()})
  if(is(x, 'list')){x <- x$Results}
  if(identical(yvar, 'default')){yvar <- c('sensitivity', 'specificity', 'correlation')}
  if(any(colnames(x) == 'cor')){colnames(x)[colnames(x) == 'cor'] <- 'correlation'}
  if(!is.null(yadd)){yvar <- c(yvar, yadd)}
  cents <- c('Strength', 'EI')
  xvar0 <- c('N', 'p', 'm2', 'index', 'iter', 'type', 'network')
  yvar0 <- setdiff(c(colnames(x), cents, 'fdr'), xvar0)
  yvar <- match.arg(tolower(yvar), unique(tolower(yvar0)), several.ok = TRUE)
  if('fdr' %in% tolower(yvar)){
    x$FDR <- 1 - x$precision
    yvar[tolower(yvar) == 'fdr'] <- 'FDR'
  }
  if(any(tolower(cents) %in% tolower(yvar))){
    for(i in which(tolower(cents) %in% tolower(yvar))){
      yvar <- c(yvar[-which(tolower(yvar) == tolower(cents[i]))],
                colnames(x)[grepl(cents[i], colnames(x))])
    }
  }
  yvar <- gsub('EI', 'ExpectedInfluence', capitalize(yvar))
  colnames(x) <- gsub('EI', 'ExpectedInfluence', capitalize(colnames(x)))
  if(!is.null(by)){
    by <- match.arg(tolower(by), c(tolower(xvar0), 'pairwise', 'interactions', 'beta', 'kappa', 'pcc'), several.ok = TRUE)
    if('pcc' %in% by){by <- gsub('pcc', 'PCC', by)}
    by <- capitalize(by)
    if(by %in% c('Pairwise', 'Interactions')){
      x <- subset(x, Type == by)
      by <- switch(2 - (length(unique(x$Network)) > 1), 'Network', NULL)
    }
  }
  if(length(unique(x$Type)) != 1 & !identical(by, 'Type')){warning('Interactions and pairwise parameters aggregated; try setting by = "type"')}
  if(length(unique(x$Network)) != 1 & !identical(by, 'Network')){warning('Multiple networks aggregated; try setting by = "network"')}
  if(any(colnames(x) == 'N')){colnames(x)[colnames(x) == 'N'] <- 'nCases'}
  if(length(unique(x$nCases)) == 1){stop('Must have used more than one sample size to plot')}
  #FUN <- bootnet:::plot.netSimulator
  args0 <- setdiff(formalArgs(netsimulator), '...')
  allargs <- formals(netsimulator)[intersect(names(formals(netsimulator)), args0)]
  allargs <- replace(allargs, c('x', 'yvar', 'xlab'), list(x = x, yvar = yvar, xlab = xlab))
  if(!is.null(by)){allargs$yfacet <- by}
  args <- args[intersect(names(args), args0[-1])]
  if(length(args) > 0){allargs <- replace(allargs, names(args), args)}
  g <- do.call(netsimulator, replace(allargs, 'print', FALSE))
  if(!is.null(hline) & !identical(hline, FALSE)){
    g <- g + ggplot2::geom_hline(yintercept = hline, linetype = 2,
                                 colour = 'red', alpha = .3)
  }
  if(!is.null(title)){g <- g + ggplot2::ggtitle(title)}
  return(g)
}

#' @rdname plotPower
#' @export
plot.mnetPower <- function(x, by = 'type', yvar = 'default', yadd = NULL, hline = .8,
                           xlab = 'Number of cases', title = NULL, ...){
  args <- as.list(match.call())[-1]
  do.call(plotPower, args)
}

#' Plot temporal and contemporaneous networks in the same window
#'
#' Designed for easy-to-use plotting with temporal networks. Essentially just a
#' wrapper for running \code{\link{plotNet}} twice---once for a temporal
#' network, and again for a contemporaneous network---and plotting the two
#' networks in the same window. Good for a quick glance at results from a SUR
#' network. Also compatible with \code{\link{mlGVAR}} and \code{\link{lmerVAR}}
#' outputs, although can only plot two networks in the same window.
#' \code{\link{plotNet3}} can be used to plot 3 networks.
#'
#' @param object Output from \code{\link{fitNetwork}}, specifically with a SUR
#'   model.
#' @param whichNets Vector of length 2 indicating which networks to plot.
#'   \code{"beta"} and \code{"temporal"} both refer to the unstandardized
#'   temporal network coefficients, while \code{"PDC"} refers to the
#'   standardized temporal network. \code{"PCC"} and \code{"contemporaneous"}
#'   both refer to the standardized residual covariance matrix (the
#'   contemporaneous network). If the \code{object} is fitted with
#'   \code{\link{mlGVAR}} or \code{\link{lmerVAR}}, then \code{"between"} is
#'   also an option for plotting the between-subjects network.
#' @param whichTemp Which version of the temporal network should be plotted,
#'   either \code{"temporal"} or \code{"PDC"}. This argument is ignored if
#'   \code{whichNets} is not \code{NULL}.
#' @param titles Character vector of length 2 where custom names for the two
#'   network plots can be supplied.
#' @param ... Additional arguments.
#'
#' @return Returns two network plots side by side.
#' @export
#'
#' @seealso \code{\link{fitNetwork}}
#'
#' @examples
#' x <- fitNetwork(gvarDat, lags = TRUE)
#' plotNet2(x)
plotNet2 <- function(object, whichNets = NULL, whichTemp = c("temporal", "PDC"),
                     titles = c("PDC ", "PCC "), ...){
  whichTemp <- match.arg(whichTemp)
  if("lmerVAR" %in% names(attributes(object))){whichTemp <- "temporal"}
  if(!is.null(whichNets)){
    tt <- ifelse(all(tolower(whichNets) == "all"), 3, 2)
    if(tt == 3){whichNets <- c(whichTemp, "contemporaneous", "between")}
    if(tt == 2){
      l <- averageLayout(plotNet(object, which.net = whichNets[1], plot = FALSE),
                         plotNet(object, which.net = whichNets[2], plot = FALSE))
    } else if(tt == 3){
      l <- averageLayout(plotNet(object, which.net = whichNets[1], plot = FALSE),
                         plotNet(object, which.net = whichNets[2], plot = FALSE),
                         plotNet(object, which.net = whichNets[3], plot = FALSE))
      if(length(titles) == 2){titles <- c("Temporal Effects", "Contemporaneous Effects", "Between-Subject Effects")}
    }
  } else {
    tt <- 2
    l <- averageLayout(plotNet(object, which.net = whichTemp, plot = FALSE),
                       plotNet(object, which.net = "contemporaneous", plot = FALSE))
  }
  layout(t(1:tt))
  if(tt == 2){
    if(all(titles == c("PDC ", "PCC "))){
      if("lmerVAR" %in% names(attributes(object)) & all(whichNets == c("temporal", "contemporaneous"))){
        title1 <- "Temporal Effects"
        title2 <- "Contemporaneous Effects"
        titles <- c(title1, title2)
      } else {
        title1 <- ifelse(whichTemp == "temporal", "Temporal Effects", "Partial Directed Correlations")
        title2 <- "Partial Contemporaneous Correlations"
      }
    } else {
      title1 <- titles[1]
      title2 <- titles[2]
    }
  }
  if(!is.null(whichNets)){
    if(length(titles) == 2){if(all(titles == c("PDC ", "PCC "))){titles <- whichNets}}
    for(i in 1:tt){
      plotNet(object, which.net = whichNets[i], layout = l, title = titles[i], ...)
    }
  } else {
    plotNet(object, which.net = whichTemp, layout = l, title = title1, ...)
    plotNet(object, which.net = "contemporaneous", layout = l, title = title2, ...)
  }
}

#' Plot temporal, contemporaneous, and between-subject networks
#'
#' Quick, easy plotting for \code{\link{mlGVAR}} and \code{\link{lmerVAR}}
#' output. Allows one to plot three networks in the same window: temporal,
#' contemporaneous, and between-subject.
#'
#' @param object Output from \code{\link{mlGVAR}} or \code{\link{lmerVAR}}.
#' @param ... Additional arguments.
#' @param nets Character vector of length 3 indicating which networks to plot,
#'   and in what order. Same options as for \code{which.net} in
#'   \code{\link{plotNet}}.
#' @param titles If \code{TRUE}, then uses default titles for temporal,
#'   contemporaneous, and between-subject networks. If \code{FALSE}, then no
#'   titles will be used. Can also be a character vector to provide custom plot
#'   titles.
#' @param l A numeric value to indicate a type of pane layout.
#' @param label Can include a character string for text annotation.
#' @param xpos Numeric, x-axis value for where the text annotation should go.
#'   Between 0 and 1.
#' @param ypos numeric, y-axis value for where the text annotation should go.
#'   Between 0 and 1.
#'
#' @return Returns 3 network plots.
#' @export
#'
#' @seealso \code{\link{mlGVAR}, \link{lmerVAR}}
#'
#' @examples
#' \donttest{
#' x <- mlGVAR(mlgvarDat, 'M')
#' plotNet3(x)
#' }
plotNet3 <- function(object, ..., nets = c('temporal', 'contemporaneous', 'between'),
                     titles = TRUE, l = 3, label = NULL, xpos = 0, ypos = .5){
  args0 <- list(...)
  avlay <- function(..., which.net = 'temporal', threshold = FALSE,
                    collapse = FALSE, net = NULL, l = NULL, args = NULL){
    if(is.null(l)){
      x <- list(...)
      stopifnot(length(x) > 0)
      #if(length(x) == 1){x <- x[[1]]} else if(collapse){x <- appd(x)}
      args0 <- list(x = NA, which.net = which.net,
                    threshold = threshold, plot = FALSE)
      if(!is.null(net)){args0$which.net <- net}
      args <- append(args[setdiff(names(args), names(args0))], args0)
      averageLayout(lapply(x, function(z) do.call(
        plotNet, replace(args, 'x', list(x = z)))))
    } else if(identical(l, 1) | identical(l, 2) | identical(l, 3) | identical(l, 4)){
      layout(switch(l, t(1:2), t(matrix(1:4, 2, 2)),
                    matrix(c(1, 1, 2, 2, 4, 3, 3, 4), nrow = 2, ncol = 4, byrow = T),
                    matrix(c(4, 1, 1, 4, 2, 2, 3, 3), nrow = 2, ncol = 4, byrow = T)))
    }
  }
  avlay(l = l)
  args0$x <- object
  stopifnot(length(nets) == 3)
  mains <- c('Temporal network', 'Partial Contemporaneous Correlations', 'Between-subjects network')
  invisible(lapply(1:3, function(i){
    args <- replace(args0, 'which.net', list(which.net = nets[i]))
    if(isTRUE(titles)){
      args$title <- mains[i]
    } else if(is.character(titles)){
      args$title <- titles[i]
    }
    do.call(plotNet, args)
  }))
  if(!is.null(label)){
    plot.new()
    text(x = xpos, y = ypos, labels = label)
  }
}

##### g_legend
g_legend <- function(a.gplot, silent = TRUE){
  fun <- switch(2 - silent, suppressWarnings, function(x){return(x)})
  tmp <- fun(ggplot_gtable(ggplot_build(a.gplot)))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

##### covNet
covNet <- function(object, mnet = TRUE, threshold = .05){
  if(isTRUE(threshold)){threshold <- .05}
  if(class(mnet) == "character"){
    object$mods0$covariates$Bcovs <- list(object$mods0$Bm)
    names(object$mods0$covariates$Bcovs) <- mnet
    data <- object$mods0$dat
    mnet <- FALSE
  } else {
    data <- if(mnet){object$mods0$dat} else {object$data}
    data <- data.frame(data, object$mods0$covariates$covs)
  }
  if(mnet){if(!"mnet" %in% names(object)){mnet <- FALSE}}
  cs <- length(object$mods0$covariates$Bcovs)
  bb <- if(mnet){object$mnet$adjMat} else {object$adjMat}
  cp <- ncol(bb)
  cadj <- cedges <- matrix(0, cp + cs, cp + cs)
  cd <- matrix(FALSE, cp + cs, cp + cs)
  cadj[1:cp, 1:cp] <- bb
  dimnames(cadj) <- dimnames(cedges) <- dimnames(cd) <- rep(list(colnames(data)), 2)
  if(mnet){cd[1:cp, 1:cp] <- object$mnet$d}
  for(i in 1:cs){
    np <- nrow(object$mods0$covariates$Bcovs[[i]])
    if(threshold != FALSE){
      cadj[cp + i, 1:np] <- ifelse(object$mods0$covariates$Bcovs[[i]][, 4] <= threshold,
                                   object$mods0$covariates$Bcovs[[i]][, 1], 0)
    } else {
      cadj[cp + i, 1:np] <- object$mods0$covariates$Bcovs[[i]][, 1]
    }
  }
  p <- ifelse(mnet, cp - 1, cp)
  if("modEdges" %in% names(object)){
    if(mnet){cedges[1:cp, 1:cp] <- object$mnet$modEdges} else {cedges[1:cp, 1:cp] <- object$modEdges}
    cedges[-c(1:cp), 1:p] <- 1
  }
  cd[-c(1:cp), 1:p] <- TRUE
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
  out <- list(adjMat = cadj, edgeColors = getEdgeColors(cadj), modEdges = cedges, d = cd, data = data)
  if(all(cedges == 0)){out$modEdges <- NULL}
  attributes(out)$mnet <- mnet
  attributes(out)$threshold <- threshold
  attributes(out)$covs <- ifelse(mnet, cs + 1, cs)
  class(out) <- c('list', 'covNet')
  return(out)
}
