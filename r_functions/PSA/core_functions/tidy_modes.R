#' @title Prepare modes of variation for ggplotting
#' @param ... Data frames of covariates with rows corresponding exactly to the input data. Column names should not be the same as column names of the input.
#' @param dimensions Supply the number of vertices after the merge corresponding to the modes of variation. By default all dimensions are included.
# eventually could be new registered methods with 'broom' package

tidy_modes <- function(result, ..., dimensions = NULL){
  if (is.null(dimensions)){dimensions <- vapply(result$vertices[-1], nrow, FUN.VALUE = 2)}
  lapply(result$modes[paste0("r=", dimensions)], function(m){cbind(as.data.frame(m), ...)}) |>
    dplyr::bind_rows(.id = "Dimension") |>
    dplyr::mutate(Dimension = factor(Dimension, levels = paste0("r=", dimensions), ordered = TRUE))
}

ggmodes <- function(result, ..., dimensions = NULL,
                    geom_line_args = list(alpha = 0.2)){
  modes <- tidy_modes(result, ..., dimensions = dimensions)
  baseplt <- tidyr::pivot_longer(modes, colnames(result$pts[[1]]), names_to = "Component") |>
    ggplot2::ggplot(ggplot2::aes(x = Component, y = value)) +
    ggplot2::facet_wrap(dplyr::vars(Dimension), scales="free")

  line <- do.call(ggplot2::geom_line, geom_line_args)

  baseplt + line +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(angle = 90)) +
    ggplot2::ggtitle("Modes")
}


tidy_scores <- function(result, ..., dimensions = NULL){
  if (is.null(dimensions)){dimensions <- vapply(result$vertices[-1], nrow, FUN.VALUE = 2)}
  scores <- dplyr::bind_cols(result$scores[paste0("r=", dimensions)], .name_repair = "minimal")
  dplyr::bind_cols(..., scores, .name_repair = "minimal")
}


ggscores <- function(result, ..., dimensions = NULL,
                     ggpairs_args = list(progress = FALSE, upper = list(continuous = "points"), legend = c(2, 1))){
  scores <- tidy_scores(result, ..., dimensions = dimensions)
  do.call(GGally::ggpairs, c(list(data = scores, columns = (ncol(...) + 1):ncol(scores)), ggpairs_args)) +
    ggplot2::ggtitle("Scores")
}

#' @title Plot scores from two merges
#' @param result Output from psoa or pssa
#' @param ... Passed to tidy_scores
#' @param dimensions Passed to tidy_scores. Note dimensions is a misnomer - it is actually the number of vertices.
#' @param mapping Create using `ggplot2::aes()`, do not specify `x` and `y` aesthetics - these are created from `dimensions` (possibly hacky, expoiting the fact that aes() creates a list).
#' @description Scatter plot of scores. Along the border of the plot the component names making up the vertices involved are printed.
#' The top border has the components that make up the positive-direction vertex of the merge specified by `y` in `mapping`.
#' The bottom border has the components that make up the negative-direction vertex of the merge.
#'
#' For PSSA scores depend on the ratio of mass in each of the vertices in the merge, and also on the total mass in both vertices (this could potentially be removed later for ).

ggscores_pairs <- function(result, ..., dimensions, mapping = NULL){
  stopifnot(length(dimensions) == 2)
  dimensiontext <- paste0("r=", dimensions)
  scores <- tidy_scores(result, ..., dimensions = dimensions)
  mdirX <- result$mergedirections[dimensiontext[[1]], ]
  mdirY <- result$mergedirections[dimensiontext[[2]], ]
  Xposcomponents <- colnames(result$mergedirections)[mdirX>0]
  Xnegcomponents <- colnames(result$mergedirections)[mdirX<0]
  Yposcomponents <- colnames(result$mergedirections)[mdirY>0]
  Ynegcomponents <- colnames(result$mergedirections)[mdirY<0]
  Xrange <- range(scores[, dimensiontext[[1]]])
  Yrange <- range(scores[, dimensiontext[[2]]])
  Xtenpcnt <- 0.1 * (Xrange[[2]] - Xrange[[1]])
  Ytenpcnt <- 0.1 * (Yrange[[2]] - Yrange[[1]])
  Xposlabels <- data.frame(label = Xposcomponents,
                           x = Xrange[[2]] + Xtenpcnt,
                           y = seq_c(Yrange[[1]], Yrange[[2]], length.out = length(Xposcomponents)))
  Xneglabels <- data.frame(label = Xnegcomponents,
                           x = Xrange[[1]] - Xtenpcnt,
                           y = seq_c(Yrange[[1]], Yrange[[2]], length.out = length(Xnegcomponents)))

  Yposlabels <- data.frame(label = Yposcomponents,
                           y = Yrange[[2]] + Ytenpcnt,
                           x = seq_c(Xrange[[1]] , Xrange[[2]], length.out = length(Yposcomponents)))
  Yneglabels <- data.frame(label = Ynegcomponents,
                           y = Yrange[[1]] - Ytenpcnt,
                           x = seq_c(Xrange[[1]], Xrange[[2]], length.out = length(Ynegcomponents)))

  labels <- dplyr::bind_rows(Xposlabels, Xneglabels, Yposlabels, Yneglabels)

  basemapping <- aes_all(c("x", "y"))
  basemapping[[1]] <- sym(dimensiontext[[1]])
  basemapping[[2]] <- sym(dimensiontext[[2]])
  ggplot(scores, mapping = basemapping) +
    geom_point(mapping = mapping) +
    geom_text(aes(x = x, y = y, label = label), data = labels)
}
