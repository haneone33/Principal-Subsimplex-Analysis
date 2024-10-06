#' @title Plot the backwards mean
#' @description I *think* this is equivalent to the zero-dimensional result, i.e. `result$vertices["r=1"]]`
#' @import ggplot2
#'
#' @param result Output of psoa() pssa() or similar.
#' @param plot If TRUE, the makes a parallel coordinates plot using `ggplot2`. Otherwise returns the backward mean.
ggbmean <- function(result, plot = TRUE){
  bmean <- result$vertices$`r=1`
  if (plot){
  data.frame(bmean) |>
    tidyr::pivot_longer(tidyr::everything()) |>
    ggplot2::ggplot() +
    ggplot2::geom_line(aes(x = .data$name, y = .data$value, group = 1)) +
    ggplot2::scale_x_discrete(guide = guide_axis(angle = 90)) +
    ggplot2::ggtitle("Backwards Mean")
  } else {
    bmean
  }
}
