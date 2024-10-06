#' @title Plot a dendrogram of the merges.
#' @description Converts results to a `dendrogram` (see [`stats::as.dendrogram()`]) object and plots it.
#' @import graphics
#'
#' @param result The `dendrogram.input` value from [`psa()`].
#' @param plot logical; if TRUE, plot the dendrogram output.
#' @param edge.root logical; if TRUE, draw an edge to the root node.
#' @param edgePar a list of plotting parameters for edges. See [`graphics::segments()`].
#' @param horiz logical; if TRUE, draw the dendrogram horizontally.
#' @param nodeheight a string specifying node height.
#' "step" for heights given by the number of vertices after the merge
#' or "rmse scores" for heights given by the RMSE of the scores of the merge.
#' @param colors a vector of colors of nodes
#' @param ... Passed to [`graphics::plot()`].

#' @export
plotdendrogram2 <- function(result, plot = TRUE, edge.root = TRUE,
                            edgePar = list(p.lty = 0, t.cex = 1, t.adj = c(1, 1)),
                            horiz = TRUE, nodeheight = "step",
                            colors = NULL, ...){
  result = result$dendrogram.input

  startvertices <- result$vertices[[1]]
  if (is.null(rownames(startvertices))){
    rownames(startvertices) <- 1:nrow(startvertices)
  }
  nodes <- lapply(1:nrow(startvertices), function(x){
    aleaf <- as.integer(x)
    attributes(aleaf) <- list(leaf = TRUE, class = "dendrogram",
                              members = as.integer(1), height = 0,
                              label = rownames(startvertices)[[x]])
    aleaf})
  nodecoords <- matrix(NA, ncol = 2, nrow = nrow(result$merges)) #for adding text later
  colnames(nodecoords) <- c("height", "midpoint")
  rownames(nodecoords) <- rownames(result$merges)
  for (i in 2:nrow(result$merges)){
    minfo <- result$merges[i, ]
    v1 <- nodes[[minfo$v1]]
    attributes(v1)$edgetext <- round(minfo$w, 2)
    v2 <- nodes[[minfo$v2]]
    attributes(v2)$edgetext <- round(1 - minfo$w, 2)
    ## Put smaller cluster on the top of the dendrogram
    if(attr(v1, 'members') < attr(v2, 'members')){
      v = v2
      v2 = v1
      v1 = v
    }
    if (nodeheight == "step"){height <- i-1}
    else if (nodeheight == "rmse scores"){height <- sqrt(mean(result$scores[[i]]^2))}
    newnode <- merge(v1, v2,
                     height = height,
                     adjust = "none")
    attributes(newnode)$label <- rownames(result$merges)[[i]]
    nodecoords[i, ] <- c(height, attributes(newnode)$midpoint)
    nodes <- c(list(newnode), nodes[-unlist(minfo[c("v1", "v2")])])
  }
  dend <- nodes[[1]]

  ## Put smaller cluster on the top of the dendrogram
  # dend <- sort_dendrogram(dend)
  ## Color labels
  if(!is.null(colors)){
    dendextend::labels_colors(dend) = colors[length(colors):1]
  }

  if (plot){
    graphics::plot(dend, edge.root = edge.root,
                   edgePar = edgePar,
                   horiz = horiz, axes = FALSE, ...)
    if (nodeheight == "step"){
      axis(side = 2-horiz,
           at=seq(0, nrow(result$merges) - 1),
           labels = as.character(seq(nrow(result$merges), 1, by = -1)))
      axislab <- "Number of Vertices"
      title(xlab = switch(1+horiz, NULL, axislab),
            ylab = switch(1+!horiz, NULL, axislab))
      grid(nx = switch(1+horiz, NA, NULL), ny = switch(1+!horiz, NA, NULL))
    } else {
      nodecoords <- nodecoords[seq(nrow(nodecoords), by = -1, length.out = 4), ]
      text(x = nodecoords[, "height"],
           y = 1 + nodecoords[, "midpoint"],
           labels = rownames(nodecoords))
      axis(side = 2-horiz)
      axislab <- "rmse scores"
      title(xlab = switch(1+horiz, NULL, axislab),
            ylab = switch(1+!horiz, NULL, axislab))
    }
  }

  invisible(dend)
}
