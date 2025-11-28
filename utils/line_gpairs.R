line_lowerfun <- function(data, mapping, col.point, col.line, shape, size){
  ggplot(data = data, mapping = mapping) +
    geom_path(col = col.line, linewidth = 0.5) +
    geom_point(col = col.point, stroke = 1, shape = shape, size = size)
}

jitter_diag <- function(data, mapping, col.point, shape, size){
  set.seed(1)
  ggplot(data, mapping) +
    geom_jitter(aes(y = 0.5), height = 0.2,
                col = col.point, stroke = 1, shape = shape, size = size) +
    scale_y_continuous(limits = c(0,1)) +
    theme_void()
}

quantile_diag <- function(data, mapping, col.point, col.line, shape, size){
  ggplot(data, mapping) +
    scale_y_continuous(trans = 'reverse') +
    geom_path(aes(y = Depth), col = col.line, linewidth = 0.5) +
    geom_point(aes(y = Depth), col = col.point, stroke = 1, shape = shape, size = size) +
    theme_void()
}

line_gpairs <- function(data, col.point, col.line, shape, size, columns = 1:5){
  data = as.data.frame(data)
  data$Depth = diatom.df$Depth
  g = ggpairs(data, columns = columns,
              lower = list(continuous = wrap(line_lowerfun,
                                             col.point = col.point,
                                             col.line = col.line,
                                             shape = shape,
                                             size = size)),
              diag = NULL,
              upper = NULL) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 15))
  return(g)
}

#' @title Single Bar Plot for PSA Loadings or Vertices
#' @description Creating a bar plot summarizing a loading vector or a vertex vector
#' @import ggplot2
#'
#' @param v a loading vector
#' @param max.k maximum number of elements to display
#'
#' @return a ggplot object
plot_vertex <- function(v, max.k = 12){

  v = v[abs(v)>1e-10]
  if(length(v) > max.k){
    threshold = sort(abs(v), decreasing = T)[max.k]
    v = v[abs(v) >= threshold]
  }

  v = sort(v, decreasing = F)
  df = data.frame(variable = factor(names(v), levels = names(v)),
                  value = v,
                  vertex = factor(v>0, levels = c(T,F)))
  rownames(df) = NULL

  g = ggplot(data = df, aes(x = .data$variable, y = .data$value, fill = .data$vertex)) +
    geom_bar(stat = 'identity') +
    coord_flip() +
    theme_bw() +
    theme(axis.ticks = element_blank()) +
    labs(x=NULL, y=NULL) +
    theme(legend.position = 'none') +
    theme(plot.title = element_text(hjust = 0.5))

  if(min(v) >= 0){
    g = g + scale_fill_manual(values = c('darkgray','darkgray'))
  }

  return(g)
}

#' @title Loading Plot for PSA
#' @description Creating a bar plot summarizing a loading vectors of PSA
#' @import ggplot2
#'
#' @param psa.res output of `psa()`
#' @param k number of loading vectors to display
#' @inheritParams plot_bar
#'
#' @return a ggplot object
#'
#' @export
plot_loadings <- function(V, k = 4, max.k = 12){

  ls = list()
  for(i in 1:k){
    ls[[i]] = plot_vertex(V[,i], max.k) + ggtitle(colnames(V)[i])
  }
  g = cowplot::plot_grid(plotlist = ls, nrow = 1, align = 'v', axis = 'l')
  return(g)
}
