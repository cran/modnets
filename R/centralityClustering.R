#' Create table of centrality values or clustering coefficients
#'
#' Mimics the output of the
#' \code{\link[qgraph:centralityTable]{qgraph::centralityTable}} and
#' \code{\link[qgraph:clusteringTable]{qgraph::clusteringTable}} functions. The
#' purpose of revising these function was to make them compatible with outputs
#' from the \code{modnets} package.
#'
#' For \code{\link{centTable}}, centrality values can be computed for the matrix
#' of interactions within a temporal network.
#'
#' @param Wmats Output from one of the primary \code{modnets} functions.
#' @param scale Logical. Determines whether to standardize values within each
#'   measure (i.e., convert to z-scores).
#' @param which.net Only applies to SUR networks, as well as those fit with the
#'   \code{\link{mlGVAR}} function. Character string to indicate which type of
#'   network to compute centrality values for. Options are \code{"temporal"} for
#'   the temporal network, \code{"contemporaneous"} for the contemporaneous
#'   network, \code{"PDC"} for the partial directed correlation network, and
#'   \code{"interactions"} for the temporal interaction network.
#' @param labels Character vector to input the names of the nodes. If left
#'   \code{NULL}, the function defaults to the node names specified by the
#'   model.
#' @param relative Logical. Determines whether to scale values within each
#'   measure relative to the largest value within that measure.
#' @param weighted Logical. If \code{TRUE} then results are converted to an
#'   unweighted network.
#' @param signed Logical. Determines whether to ignore the signs of edges or
#'   not. Primarily affects the output for expected influence statistics.
#'
#' @return A table containing the names of nodes, measures of node centrality,
#'   and their corresponding centrality values or clustering coefficients.
#' @export
#' @name CentralityAndClustering
#'
#' @seealso \code{\link{centAuto}, \link{clustAuto}, \link{centPlot},
#'   \link{clustPlot}, \link{plotCentrality},
#'   \link[qgraph:centralityTable]{qgraph::centralityTable},
#'   \link[qgraph:clusteringTable]{qgraph::clusteringTable}}
#'
#' @examples
#' x <- fitNetwork(gvarDat, 'M', lags = TRUE)
#'
#' clustTable(x)
#' centTable(x, which.net = 'interactions')
centTable <- function(Wmats, scale = TRUE, which.net = "temporal", labels = NULL,
                      relative = FALSE, weighted = TRUE, signed = TRUE){
  if(isTRUE(attr(Wmats, "mlGVAR"))){
    Wmats <- switch(which.net, between = Wmats$betweenNet, Wmats$fixedNets)}
  if("SURnet" %in% c(names(Wmats), names(attributes(Wmats)))){
    if("SURnet" %in% names(Wmats)){Wmats <- Wmats$SURnet}
    if(is.numeric(which.net)){which.net <- c("t", "c", "p", 'i')[which.net]}
    which.net <- match.arg(tolower(which.net), c("temporal", "contemporaneous", "pdc", "interactions"))
    Wmats <- Wmats[[ifelse(which.net == "contemporaneous", "contemporaneous", ifelse(which.net == 'interactions', 'interactions', 'temporal'))]]
    if(which.net == "pdc"){Wmats <- Wmats$PDC}
    if(which.net == 'interactions'){names(Wmats)[1] <- 'adjMat'}
  } else if(startsWith(tolower(which.net), 'i')){stop('Interaction centrality only supported for temporal networks.')}
  if("adjMat" %in% names(Wmats)){Wmats <- t(Wmats$adjMat)}
  if(any(grepl("lag", dimnames(Wmats)))){dimnames(Wmats) <- lapply(dimnames(Wmats), function(z) gsub("[.]lag1.*|[.]y$", "", z))}
  if(!is.list(Wmats)){Wmats <- list(Wmats)}
  if(any(sapply(Wmats, ncol) == 1)){stop("Not supported for single-node graphs")}
  #names(Wmats) <- qgraph:::fixnames(Wmats, "graph ")
  names(Wmats) <- fnames(Wmats, 'graph ')
  centOut <- lapply(Wmats, centAuto, which.net = which.net, weighted = weighted, signed = signed)
  for(g in seq_along(centOut)){
    if(!is(centOut[[g]], "centrality_auto")){
      #names(centOut[[g]]) <- qgraph:::fixnames(centOut[[g]], "type ")
      names(centOut[[g]]) <- fnames(centOut[[g]], 'type ')
      for(t in seq_along(centOut[[g]])){
        if(!is.null(labels)){
          centOut[[g]][[t]][["node.centrality"]][["node"]] <- labels
        } else if(!is.null(rownames(centOut[[g]][[t]][["node.centrality"]]))){
          centOut[[g]][[t]][["node.centrality"]][["node"]] <- rownames(centOut[[g]][[t]][["node.centrality"]])
        } else {
          centOut[[g]][[t]][["node.centrality"]][["node"]] <- paste("Node", seq_len(nrow(centOut[[g]][[t]][["node.centrality"]])))
        }
        centOut[[g]][[t]]$node.centrality$graph <- names(centOut)[g]
        centOut[[g]][[t]]$node.centrality$type <- names(centOut[[g]])[t]
      }
    } else {
      centOut[[g]]$node.centrality$graph <- names(centOut)[g]
      if(!is.null(labels)){
        centOut[[g]][["node.centrality"]][["node"]] <- labels
      } else if(!is.null(rownames(centOut[[g]][["node.centrality"]]))){
        centOut[[g]][["node.centrality"]][["node"]] <- rownames(centOut[[g]][["node.centrality"]])
      } else {
        centOut[[g]][["node.centrality"]][["node"]] <- paste("Node", seq_len(nrow(centOut[[g]][["node.centrality"]])))
      }
    }
  }
  isList <- sapply(centOut, function(x) !"centrality_auto" %in% class(x))
  if(any(isList)){
    for(l in which(isList)){centOut <- c(centOut, centOut[[l]])}
    centOut <- centOut[-which(isList)]
  }
  for(i in seq_along(centOut)){
    if(relative | scale){
      if(relative & scale){warning("Using 'relative' and 'scale' together is not recommended")}
      for(j in which(sapply(centOut[[i]][["node.centrality"]], mode) == "numeric")){
        if(scale){
          #centOut[[i]][["node.centrality"]][, j] <- qgraph:::scale2(centOut[[i]][["node.centrality"]][, j])
          centOut[[i]][["node.centrality"]][, j] <- scaleNA(centOut[[i]][["node.centrality"]][, j])
        }
        if(relative){
          mx <- max(abs(centOut[[i]][["node.centrality"]][, j]), na.rm = TRUE)
          if(mx != 0){centOut[[i]][["node.centrality"]][, j] <- centOut[[i]][["node.centrality"]][, j]/mx}
        }
        attributes(centOut[[i]][["node.centrality"]][, j]) <- NULL
      }
    }
  }
  wideCent <- plyr::rbind.fill(lapply(centOut, "[[", "node.centrality"))
  if(is.null(wideCent$type)){wideCent$type <- NA}
  longCent <- reshape2::melt(wideCent, variable.name = "measure", id.var = c("graph", "type", "node"))
  if(any(is.nan(longCent$value))){warning("NaN detected in centrality measures. Try relative = FALSE")}
  return(longCent)
}

#' @rdname CentralityAndClustering
#' @export
clustTable <- function(Wmats, scale = TRUE, labels = NULL,
                       relative = FALSE, signed = TRUE){
  if("SURnet" %in% c(names(Wmats), names(attributes(Wmats)))){
    if("SURnet" %in% names(Wmats)){Wmats <- Wmats$SURnet}
    Wmats <- Wmats$contemporaneous$adjMat
  } else if("adjMat" %in% names(Wmats)){
    Wmats <- Wmats$adjMat
  }
  if(any(grepl("lag", dimnames(Wmats)))){
    dimnames(Wmats) <- lapply(dimnames(Wmats), function(z) gsub("[.]lag1.*|[.]y$", "", z))
  }
  if(!is.list(Wmats)){Wmats <- list(Wmats)}
  if(any(sapply(Wmats, ncol) == 1)){stop("Not supported for single-node graphs")}
  syms <- sapply(Wmats, isSymmetric)
  if(any(!syms)){
    if(all(!syms)){stop("No symmetrical graphs detected")}
    warning(paste0(sum(!syms), " Nonsymmetrical graph", ifelse(sum(!syms) > 1, "s ", " "), "removed"))
    Wmats <- Wmats[-which(!syms)]
  }
  #names(Wmats) <- qgraph:::fixnames(Wmats, "graph ")
  names(Wmats) <- fnames(Wmats, 'graph ')
  clustOut <- lapply(Wmats, clustAuto)
  for(g in seq_along(clustOut)){
    if(!is(clustOut[[g]], "clustcoef_auto")){
      #names(clustOut[[g]]) <- qgraph:::fixnames(clustOut[[g]], "type ")
      names(clustOut[[g]]) <- fnames(clustOut[[g]], 'type ')
      for(t in seq_along(clustOut[[g]])){
        if(!is.null(labels)){
          clustOut[[g]][[t]][["node"]] <- labels
        } else if(!is.null(rownames(clustOut[[g]][[t]]))){
          clustOut[[g]][[t]][["node"]] <- rownames(clustOut[[g]][[t]])
        } else {
          clustOut[[g]][[t]][["node"]] <- paste("Node", seq_len(nrow(clustOut[[g]][[t]])))
        }
        clustOut[[g]][[t]]$graph <- names(clustOut)[g]
        clustOut[[g]][[t]]$type <- names(clustOut[[g]])[t]
      }
    } else {
      clustOut[[g]]$graph <- names(clustOut)[g]
      if(!is.null(labels)){
        clustOut[[g]][["node"]] <- labels
      } else if(!is.null(rownames(clustOut[[g]]))){
        clustOut[[g]][["node"]] <- rownames(clustOut[[g]])
      } else {
        clustOut[[g]][["node"]] <- paste("Node", seq_len(nrow(clustOut[[g]])))
      }
    }
  }
  isList <- sapply(clustOut, function(x) !"clustcoef_auto" %in% class(x))
  if(any(isList)){
    for(l in which(isList)){clustOut <- c(clustOut, clustOut[[l]])}
    clustOut <- clustOut[-which(isList)]
  }
  for(i in seq_along(clustOut)){
    if(any(grepl("signed_", names(clustOut[[i]])))){
      clustOut[[i]] <- clustOut[[i]][, sapply(clustOut[[i]], mode) != "numeric" | grepl("signed_", names(clustOut[[i]])) == signed]
      names(clustOut[[i]]) <- gsub("signed_", "", names(clustOut[[i]]))
    }
    names(clustOut[[i]]) <- gsub("clust", "", names(clustOut[[i]]))
    if(relative | scale){
      if(relative & scale){warning("Using 'relative' and 'scale' together is not recommended")}
      for(j in which(sapply(clustOut[[i]], mode) == "numeric")){
        if(scale){
          #clustOut[[i]][, j] <- qgraph:::scale2(clustOut[[i]][, j])
          clustOut[[i]][, j] <- scaleNA(clustOut[[i]][, j])
        }
        if(relative){
          mx <- max(abs(clustOut[[i]][, j]), na.rm = TRUE)
          if(mx != 0){clustOut[[i]][, j] <- clustOut[[i]][, j]/mx}
        }
        attributes(clustOut[[i]][, j]) <- NULL
      }
    }
  }
  WideCent <- plyr::rbind.fill(clustOut)
  if(is.null(WideCent$type)){WideCent$type <- NA}
  LongCent <- reshape2::melt(WideCent, variable.name = "measure", id.var = c("graph", "type", "node"))
  return(LongCent)
}

#' Node centrality, clustering coefficients, and shortest path lengths
#'
#' Mimics the \code{\link[qgraph:centrality_auto]{qgraph::centrality_auto}} and
#' \code{\link[qgraph:clustcoef_auto]{qgraph::clustcoef_auto}} functions. The
#' purpose of amending these functions was to make them compatible with outputs
#' from the \code{modnets} package. The main use of these functions is as the
#' engines for the \code{\link{centTable}} and \code{\link{clustTable}}
#' functions.
#'
#' Returns several node centrality statistics, edge-betweenness centrality, and
#' shortest path lengths. Betweenness and Closeness centrality are computed for
#' all types of networks, as well as edge-betweenness values and shortest path
#' lengths. For GGMs, Strength centrality and Expected Influence are also
#' computed. For SUR networks, InStrength, OutStrength, InExpectedInfluence, and
#' OutExpectedInfluence are computed instead.
#'
#' The key distinction between these functions and the
#' \code{\link[qgraph:centrality_auto]{qgraph::centrality_auto}} and
#' \code{\link[qgraph:clustcoef_auto]{qgraph::clustcoef_auto}} functions is that
#' centrality and clustering values can be computed for the matrix of
#' interactions within a temporal network.
#'
#' @param x Output from one of the primary \code{modnets} functions. Can also
#'   supply a list of network models, and the function will be applied to all
#'   models in the list.
#' @param which.net Only applies to SUR networks, as well as those fit with the
#'   \code{\link{mlGVAR}} function. Character string to indicate which type of
#'   network to compute centrality values for. Options are \code{"temporal"} for
#'   the temporal network, \code{"contemporaneous"} for the contemporaneous
#'   network, \code{"PDC"} for the partial directed correlation network, and
#'   \code{"interactions"} for the temporal interaction network.
#' @param weighted Logical. If \code{TRUE} then results are converted to an
#'   unweighted network.
#' @param signed Logical. Determines whether to ignore the signs of edges or
#'   not. Primarily affects the output for expected influence statistics.
#' @param thresholdWS Numeric threshold for the WS values.
#' @param thresholdON Numeric threshold for the Zhang values.
#'
#' @return A list containing node centrality statistics, edge-betweenness
#'   values, and shortest path lengths.
#' @export
#' @name CentClust
#'
#' @seealso \code{\link{centTable}, \link{clustTable}, \link{centPlot},
#'   \link{clustPlot}, \link{plotCentrality},
#'   \link[qgraph:centrality_auto]{qgraph::centrality_auto},
#'   \link[qgraph:clustcoef_auto]{qgraph::clustcoef_auto}}
#'
#' @examples
#' x <- fitNetwork(ggmDat, 'M')
#'
#' clustAuto(x)
#' centAuto(x, 'interactions')
centAuto <- function(x, which.net = "temporal", weighted = TRUE, signed = TRUE){
  if(isTRUE(attr(x, "mlGVAR"))){
    x <- switch(which.net, between = x$betweenNet, x$fixedNets)}
  if("SURnet" %in% c(names(x), names(attributes(x)))){
    if("SURnet" %in% names(x)){x <- x$SURnet}
    if(is.numeric(which.net)){which.net <- c("t", "c", "p", "i")[which.net]}
    which.net <- match.arg(tolower(which.net), c("temporal", "contemporaneous", "pdc", "interactions"))
    x <- x[[ifelse(which.net == "contemporaneous", "contemporaneous", ifelse(which.net == 'interactions', 'interactions', "temporal"))]]
    if(which.net == "pdc"){x <- x$PDC}
    if(which.net == 'interactions'){names(x)[1] <- 'adjMat'}
  }
  if("adjMat" %in% names(x)){x <- t(x$adjMat)}
  if(any(grepl("lag", dimnames(x)))){dimnames(x) <- lapply(dimnames(x), function(z) gsub("[.]lag1.*|[.]y$", "", z))}
  if(is.list(x)){return(lapply(x, centAuto, which.net = which.net, weighted = weighted, signed = signed))}
  if(!weighted){x <- sign(x)}
  if(!signed){x <- abs(x)}
  if(!is.matrix(x)){stop("The input network must be an adjacency or weights matrix")}
  diag(x) <- 0
  directed.gr <- ifelse(isSymmetric.matrix(object = x, tol = 0.000000000001), FALSE, TRUE)
  weighted.gr <- ifelse(all(qgraph::mat2vec(x) %in% c(0, 1)), FALSE, TRUE)
  net_qg <- qgraph::qgraph(x, diag = FALSE, labels = colnames(x), DoNotPlot = TRUE, minimum = 0)
  centr <- qgraph::centrality(net_qg)
  if(directed.gr & !weighted.gr){
    centr1 <- data.frame(cbind(
      Betweenness = centr$Betweenness, Closeness = centr$Closeness,
      InDegree = centr$InDegree, OutDegree = centr$OutDegree,
      OutExpectedInfluence = centr$OutExpectedInfluence,
      InExpectedInfluence = centr$InExpectedInfluence))
  }
  if(directed.gr & weighted.gr){
    centr1 <- data.frame(cbind(
      Betweenness = centr$Betweenness, Closeness = centr$Closeness,
      InStrength = centr$InDegree, OutStrength = centr$OutDegree,
      OutExpectedInfluence = centr$OutExpectedInfluence,
      InExpectedInfluence = centr$InExpectedInfluence))
  }
  if(!directed.gr & !weighted.gr){
    centr1 <- data.frame(cbind(
      Betweenness = centr$Betweenness/2, Closeness = centr$Closeness,
      Degree = centr$OutDegree, ExpectedInfluence = centr$OutExpectedInfluence))
  }
  if(!directed.gr & weighted.gr){
    centr1 <- data.frame(cbind(
      Betweenness = centr$Betweenness/2, Closeness = centr$Closeness,
      Strength = centr$OutDegree, ExpectedInfluence = centr$OutExpectedInfluence))
  }
  row.names(centr1) <- colnames(x)
  log <- capture.output({
    graph <- igraph::graph.adjacency(
      adjmatrix = 1 * (x != 0),
      mode = ifelse(directed.gr, "directed", "undirected"))
    comps <- igraph::components(graph)
    largcomp <- comps$membership == which.max(comps$csize)
  })
  if(sum(largcomp) < ncol(x) & sum(largcomp) > 1){
    x2 <- x[largcomp, largcomp]
    clos <- qgraph::centrality(qgraph::qgraph(
      x2, diag = FALSE, labels = colnames(x)[largcomp],
      DoNotPlot = TRUE, minimum = 0))$Closeness
    centr1$Closeness[largcomp] <- clos
    centr1$Closeness[!largcomp] <- NA
  }
  net_ig_abs <- igraph::graph.adjacency(
    adjmatrix = abs(1/x), mode = ifelse(directed.gr, "directed", "undirected"),
    weighted = ifelse(weighted.gr, list(TRUE), list(NULL))[[1]], diag = FALSE)
  edgebet <- igraph::edge.betweenness(graph = net_ig_abs, directed = directed.gr)
  el <- data.frame(igraph::get.edgelist(graph = net_ig_abs), stringsAsFactors = FALSE)
  edgebet <- merge(el, edgebet, by = 0)
  edgebet$Row.names <- NULL
  names(edgebet) <- c("from", "to", "edgebetweenness")
  edgebet <- edgebet[order(edgebet$edgebetweenness, decreasing = TRUE), ]
  ShortestPathLengths <- centr$ShortestPathLengths
  rownames(ShortestPathLengths) <- colnames(ShortestPathLengths) <- colnames(x)
  Res <- list(node.centrality = centr1, edge.betweenness.centrality = edgebet,
              ShortestPathLengths = ShortestPathLengths)
  class(Res) <- c("list", "centrality_auto")
  return(Res)
}

#' @rdname CentClust
#' @export
clustAuto <- function(x, thresholdWS = 0, thresholdON = 0){
  if("SURnet" %in% c(names(x), names(attributes(x)))){
    if("SURnet" %in% names(x)){x <- x$SURnet}
    x <- x$contemporaneous$adjMat
  } else if("adjMat" %in% names(x)){
    x <- x$adjMat
  }
  if(any(grepl("lag", dimnames(x)))){
    dimnames(x) <- lapply(dimnames(x), function(z) gsub("[.]lag1.*|[.]y$", "", z))
  }
  if(is.list(x)){
    return(lapply(x, clustAuto, thresholdWS = thresholdWS, thresholdON = thresholdWS))
  }
  dim = dim(x)
  if(is.null(dim) || length(dim) != 2){stop("adjacency is not two-dimensional")}
  if(!is.numeric(x)){stop("adjacency is not numeric")}
  if(dim[1] != dim[2]){stop("adjacency is not square")}
  if(max(abs(x - t(x)), na.rm = TRUE) > 0.000000000001){stop("adjacency is not symmetric")}
  if(min(x, na.rm = TRUE) < -1 || max(x, na.rm = TRUE) > 1){x <- x/max(abs(x))}
  weighted.gr <- ifelse(all(abs(x) %in% c(0, 1)), FALSE, TRUE)
  signed.gr <- ifelse(all(x >= 0), FALSE, TRUE)
  net_ig <- igraph::graph.adjacency(
    adjmatrix = abs(x), mode = "undirected",
    weighted = ifelse(weighted.gr, list(TRUE), list(NULL))[[1]], diag = FALSE)
  cb <- igraph::transitivity(net_ig, type = "barrat", isolates = "zero")
  #cw <- qgraph:::clustWS(x, thresholdWS)
  #cz <- qgraph:::clustZhang(x)
  #co <- qgraph:::clustOnnela(x, thresholdON)
  cw <- WS(x, thresholdWS)
  cz <- zhang(x)
  co <- onnela(x, thresholdON)
  if(!signed.gr & !weighted.gr){output <- cbind(clustWS = cw[, 1])}
  if(!signed.gr & weighted.gr){
    output <- cbind(clustWS = cw[, 1], clustZhang = cz[, 1],
                    clustOnnela = co[, 1], clustBarrat = cb)
  }
  if(signed.gr & !weighted.gr){
    output <- cbind(clustWS = cw[, 1], signed_clustWS = cw[, 2])
  }
  if(signed.gr & weighted.gr){
    output <- cbind(clustWS = cw[, 1], signed_clustWS = cw[, 2],
                    clustZhang = cz[, 1], signed_clustZhang = cz[, 2],
                    clustOnnela = co[, 1], signed_clustOnnela = co[, 2],
                    clustBarrat = cb)
  }
  output[is.na(output)] <- 0
  Res <- data.frame(output)
  class(Res) <- c("data.frame", "clustcoef_auto")
  rownames(Res) <- colnames(x)
  return(Res)
}

#' Plots for node centrality values or clustering coefficients
#'
#' Mimics the \code{\link[qgraph:centralityPlot]{qgraph::centralityPlot}} and
#' \code{\link[qgraph:clusteringPlot]{qgraph::clusteringPlot}} functions. The
#' purpose of revising this function was to make it compatible with outputs from
#' the \code{modnets} package.
#'
#' The only utility of the \code{\link{plotCentrality}} function is as an easy
#' way to combine centrality measures and clustering coefficients into a single
#' plot.
#'
#' @param Wmats Output from one of the primary \code{modnets} functions.
#' @param scale If \code{"z-scores"}, then standardized values will be plotted.
#'   If \code{"relative"}, then values will be scaled relative to the largest
#'   value on each measure. \code{"raw"} can be used to plot raw values.
#' @param which.net Only applies to SUR networks, as well as those fit with the
#'   \code{\link{mlGVAR}} function. Character string to indicate which type of
#'   network to compute centrality values for. Options are \code{"temporal"} for
#'   the temporal network, \code{"contemporaneous"} for the contemporaneous
#'   network, \code{"PDC"} for the partial directed correlation network, and
#'   \code{"interactions"} for the temporal interaction network.
#' @param include Character vector of which centrality measures to plot.
#'   \code{"Betweenness"} and \code{"Closeness"} are available for all types of
#'   network. \code{"Strength"} and \code{"ExpectedInfluence"} are only
#'   available for GGMs. And \code{"InStrength", "OutStrength",
#'   "InExpectedInfluence", "OutExpectedInfluence"} are only available for SUR
#'   networks. Defaults to \code{"all"}
#' @param labels Character vector listing the node names. If \code{NULL}, then
#'   the names specified by the model are used.
#' @param orderBy Character string specifying which measure to order values by.
#' @param decreasing Logical. Only relevant if \code{orderBy} is specified.
#'   Determines whether values are organized from highest to lowest, or vice
#'   versa.
#' @param plot Logical. Determines whether to plot the output or not.
#' @param verbose Logical. Determines whether to return a message about the plot
#'   (messages are only shown if values are scaled).
#' @param weighted See \code{\link{centTable}} or \code{\link{clustTable}}.
#' @param signed See \code{\link{centTable}} or \code{\link{clustTable}}.
#' @param centrality Character vector of centrality measures to plot. Defaults
#'   to \code{"all"}.
#' @param clustering Character vector of clustering measures to plot. Defaults
#'   to \code{"Zhang"}.
#'
#' @return A plot of centrality values or clustering coefficients for several
#'   measures.
#' @export
#' @name CentralityAndClusteringPlots
#'
#' @seealso \code{\link{centTable}, \link{clustTable}, \link{centAuto},
#'   \link{clustAuto}, \link[qgraph:centralityPlot]{qgraph::centralityPlot},
#'   \link[qgraph:clusteringPlot]{qgraph::clusteringPlot}}
#'
#' @examples
#' x <- fitNetwork(ggmDat)
#'
#' centPlot(x)
#' clustPlot(x)
#' plotCentrality(x)
centPlot <- function(Wmats, scale = c("z-scores", "raw", "raw0", "relative"),
                     which.net = "temporal", include = "all", labels = NULL,
                     orderBy = NULL, decreasing = FALSE, plot = TRUE,
                     verbose = TRUE, weighted = TRUE, signed = TRUE){
  if(isTRUE(attr(Wmats, "mlGVAR"))){
    Wmats <- switch(which.net, between = Wmats$betweenNet, Wmats$fixedNets)}
  if(is.logical(scale)){scale <- ifelse(scale, "z-scores", "raw")}
  #invisible(suppressMessages(require(ggplot2)))
  measure <- value <- node <- type <- NULL
  scale <- tryCatch({match.arg(scale)}, error = function(e){scale})
  #scale <- match.arg(scale)
  include0 <- c("Degree", "Strength", "OutDegree", "InDegree", "OutStrength",
                "InStrength", "Closeness", "Betweenness", "ExpectedInfluence",
                "OutExpectedInfluence", "InExpectedInfluence")
  if(all(tolower(include) == "all")){
    include <- include0
  } else if(isTRUE(attr(Wmats, "SURnet")) & which.net != "contemporaneous"){
    include0 <- include0[!grepl("Degree|^S|^E", include0)]
    include <- include0[grep(paste(tolower(include), collapse = "|"), tolower(include0))]
  } else if(which.net %in% c("between", "contemporaneous")){
    include0 <- include0[!grepl("Degree|^Out|^In", include0)]
    include <- include0[grep(paste(tolower(include), collapse = "|"), tolower(include0))]
  }
  include <- match.arg(include, c("Degree", "Strength", "OutDegree", "InDegree", "OutStrength",
                                  "InStrength", "Closeness", "Betweenness", "ExpectedInfluence",
                                  "OutExpectedInfluence", "InExpectedInfluence"), several.ok = TRUE)
  if(scale == "z-scores" & verbose & plot){message("Note: z-scores are shown on x-axis.")}
  if(scale == "relative" & verbose & plot){message("Note: relative centrality indices are shown on x-axis.")}
  Long <- centTable(Wmats = Wmats, scale = (scale == "z-scores"), labels = labels,
                    which.net = which.net, relative = (scale == "relative"),
                    weighted = weighted, signed = signed)
  Long <- subset(Long, measure %in% include)
  Long$measure <- factor(Long$measure)
  if(ifelse(is.null(orderBy), FALSE, ifelse(orderBy == "default", TRUE, FALSE))){
    nodeLevels <- unique(gtools::mixedsort(
      as.character(Long$node), decreasing = decreasing))
  } else if(!is.null(orderBy)){
    nodeLevels <- names(sort(tapply(
      Long$value[Long$measure == orderBy],
      Long$node[Long$measure == orderBy], mean),
      decreasing = decreasing))
  } else {
    nodeLevels <- rev(unique(as.character(Long$node)))
  }
  Long$node <- factor(as.character(Long$node), levels = nodeLevels)
  Long <- Long[gtools::mixedorder(Long$node), ]
  if(length(unique(Long$type)) > 1){
    g <- ggplot(Long, aes(x = value, y = node, group = type, colour = type))
  } else {
    g <- ggplot(Long, aes(x = value, y = node, group = type))
  }
  g <- g + geom_path() + xlab("") + ylab("") + geom_point() + theme_bw()
  if(length(unique(Long$graph)) > 1){
    g <- g + facet_grid(graph ~ measure, scales = "free")
  } else {
    g <- g + facet_grid(~measure, scales = "free")
  }
  if(scale == "raw0"){g <- g + xlim(0, NA)}
  if(plot){plot(g)} else {invisible(g)}
}

#' @rdname CentralityAndClusteringPlots
#' @export
clustPlot <- function(Wmats, scale = c("z-scores", "raw", "raw0", "relative"),
                      include = "all", labels = NULL, orderBy = NULL,
                      decreasing = FALSE, plot = TRUE, signed = TRUE,
                      verbose = TRUE){
  if(is.logical(scale)){scale <- ifelse(scale, "z-scores", "raw")}
  #invisible(suppressMessages(require(ggplot2)))
  measure <- value <- node <- type <- NULL
  scale <- match.arg(scale)
  if(scale == "z-scores" & verbose & plot){message("Note: z-scores are shown on x-axis.")}
  if(scale == "relative" & verbose & plot){message("Note: relative centrality indices are shown on x-axis.")}
  Long <- clustTable(Wmats = Wmats, scale = (scale == "z-scores"), labels = labels,
                     relative = (scale == "relative"), signed = signed)
  Long$value[!is.finite(Long$value)] <- 0
  if(all(include == "all")){include <- c("WS", "Zhang", "Onnela", "Barrat")}
  include <- match.arg(include, c("WS", "Zhang", "Onnela", "Barrat"), several.ok = TRUE)
  Long <- subset(Long, measure %in% include)
  Long$measure <- factor(Long$measure)
  if(ifelse(is.null(orderBy), FALSE, ifelse(orderBy == "default", TRUE, FALSE))){
    nodeLevels <- unique(gtools::mixedsort(
      as.character(Long$node), decreasing = decreasing))
  } else if(!is.null(orderBy)){
    nodeLevels <- names(sort(tapply(
      Long$value[Long$measure == orderBy],
      Long$node[Long$measure == orderBy], mean),
      decreasing = decreasing))
  } else {
    nodeLevels <- rev(unique(as.character(Long$node)))
  }
  Long$node <- factor(as.character(Long$node), levels = nodeLevels)
  Long <- Long[gtools::mixedorder(Long$node), ]
  if(length(unique(Long$type)) > 1){
    g <- ggplot(Long, aes(x = value, y = node, group = type, colour = type))
  } else {
    g <- ggplot(Long, aes(x = value, y = node, group = type))
  }
  g <- g + geom_path() + xlab("") + ylab("") + geom_point() + theme_bw()
  if(length(unique(Long$graph)) > 1){
    g <- g + facet_grid(graph ~ measure, scales = "free")
  } else {
    g <- g + facet_grid(~measure, scales = "free")
  }
  if(scale == "raw0"){g <- g + xlim(0, NA)}
  if(plot){plot(g)} else {invisible(g)}
}

#' @rdname CentralityAndClusteringPlots
#' @export
plotCentrality <- function(Wmats, which.net = "temporal", scale = TRUE,
                           labels = NULL, plot = TRUE, centrality = "all",
                           clustering = "Zhang"){
  if(any(c("ggm", "SURnet", "mlGVAR") %in% names(attributes(Wmats)))){Wmats <- list(net1 = Wmats)}
  if(all(sapply(Wmats, function(z) isTRUE(attr(z, "mlGVAR"))))){
    Wmats <- lapply(Wmats, function(z) switch(
      which.net, between = z$betweenNet, z$fixedNets))
  }
  if(any(grepl("ggm", lapply(Wmats, function(z) names(attributes(z)))))){which.net <- "contemporaneous"}
  if(length(unique(lapply(Wmats, checkInclude))) != 1){stop("All networks must be of the same type")}
  if(is.null(names(Wmats))){names(Wmats) <- paste0("net", seq_along(Wmats))}
  which.net <- match.arg(tolower(which.net), c("temporal", "contemporaneous", "pdc"))
  c0 <- c01 <- do.call(rbind, lapply(seq_along(Wmats), function(z){
    cbind.data.frame(
      centTable(Wmats[[z]], scale = scale, which.net = which.net,
                labels = labels), group = names(Wmats)[z])
  }))
  if(all(centrality != "all")){
    include0 <- checkInclude(Wmats[[1]], which.net = which.net)
    include0 <- include0[!grepl(ifelse(
      which.net != "contemporaneous", "Degree|^S|^E",
      "Degree|^Out|^In"), include0)]
    centrality <- include0[grep(paste(tolower(
      centrality), collapse = "|"), tolower(include0))]
    c0 <- c01 <- subset(c0, measure %in% centrality)
  }
  if(which.net == "contemporaneous" & clustering != FALSE){
    c1 <- do.call(rbind, lapply(seq_along(Wmats), function(z){
      z1 <- clustTable(Wmats[[z]], scale = scale, labels = labels)
      z1 <- z1[z1$measure == ifelse(is.logical(clustering), "Zhang", clustering), ]
      z1$measure <- "Clust. coef."
      z1$node <- as.character(z1$node)
      rownames(z1) <- 1:nrow(z1)
      return(cbind.data.frame(z1, group = names(Wmats)[z]))
    }))
    c01 <- rbind(c0, c1)
  }
  c01 <- c01[order(c01$node), ]
  c01 <- c01[order(c01$group), ]
  rownames(c01) <- 1:nrow(c01)
  #c01$node <- stringr::str_sub(c01$node, 1, 6)
  c01$node <- substr(c01$node, 1, 6)
  if(!plot){
    if(which.net == 'contemporaneous' & clustering != FALSE){
      out <- list(centrality = c0, clustering = c1, combined = c01)
    } else {
      out <- c0
    }
    return(out)
    #list2env(list(c0 = c0, c1 = c1, c01 = c01), .GlobalEnv)
  } else {
    #invisible(suppressMessages(require(ggplot2)))
    g1 <- ggplot(c01, aes(x = value, y = node, group = group, color = group, shape = group)) +
      geom_path(alpha = 1, size = 1) + geom_point(size = 2) + xlab("") + ylab("") + theme_bw() +
      facet_grid(. ~ measure, scales = "free") + scale_x_continuous(breaks = c(-1, 0, 1)) +
      theme(axis.line.x = element_line(colour = "black"),
            axis.ticks.x = element_line(colour = "black"),
            axis.ticks.y = element_line(colour = "white", size = 0),
            axis.text.y = element_text(colour = "black"),
            axis.text.x = element_text(angle = 45, colour = "black"))
    g1
  }
}
