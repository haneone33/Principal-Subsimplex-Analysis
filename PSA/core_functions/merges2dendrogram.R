#' @title Plot a dendrogram of the merges.
#' @description Converts results to a `dendrogram` (see [`stats::as.dendrogram()`]) object and plots it.
#' @param ... Passed to [`stats:::plot.dendrogram()`].
#' Height is the sqrt of the squared sum of scores of a merge.
#' @param result The return value from [`pssa()`], [`psoa()`] or [`recurse_tolowerdimension()`].
#' @param If "step", then merges have 'height' given by the number of vertices after the merge. Otherwise, if "rmse scores", the node height is the RMSE of the scores of the merge.
#' @examples

#' @export
plotdendrogram <- function(result, plot = TRUE, edge.root = TRUE,
                           edgePar = list(p.lty = 0, t.cex = 0.5), horiz = TRUE, nodeheight = "step", ...){
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
  if (plot){
  stats:::plot.dendrogram(dend, edge.root = edge.root,
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
