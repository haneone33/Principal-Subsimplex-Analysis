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
plot_loadings <- function(V, k = base::min(ncol(V),4), max.k = 12){

  ls = list()
  for(i in 1:k){
    ls[[i]] = plot_vertex(V[,i], max.k) + ggtitle(colnames(V)[i])
  }
  g = cowplot::plot_grid(plotlist = ls, nrow = 1, align = 'v', axis = 'l')
  return(g)
}
