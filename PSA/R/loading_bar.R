#' @title Loading Bar Plots for PSA
#' @description Creating a bar plot summarizing a loading vector
#' @import ggplot2
#'
#' @param v a loading vector
#' @param max.k maximum number of elements to display
#'
#' @return a ggplot object
#'
#' @export
loading_bar <- function(v, max.k = 12){

  v = v[abs(v)>1e-8]
  if(length(v) > max.k){
    threshold = sort(abs(v), decreasing = T)[max.k]
    v = v[abs(v) >= threshold]
  }

  v = sort(v, decreasing = F)
  df = data.frame(variable = factor(names(v), levels = names(v)),
                  value = v,
                  vertex = factor(v>0, levels = c(T,F)))
  rownames(df) = NULL

  g = ggplot2::ggplot(data = df, aes(x = .data$variable, y = .data$value, fill = .data$vertex)) +
    geom_bar(stat = 'identity') +
    coord_flip() +
    theme_bw() +
    theme(axis.ticks.y = element_blank()) +
    labs(x = '', y = '') +
    theme(legend.position = 'none') +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(plot.margin = unit(c(2.5,2.5,2.5,2.5), "points"))

  return(g)
}
