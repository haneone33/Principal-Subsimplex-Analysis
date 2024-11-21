point_gpairs <- function(scores, X.label, k = min(ncol(scores)-1,5)){
  scores.df = as.data.frame(scores[,1:k])
  scores.df$label = X.label
  lims = plot.lims(scores[,1:k])
  score.plot = ggpairs(scores.df, aes(col = label), columns = 1:k, upper = NULL,
                       lower = list(continuous = wrap(lowerfun, lims = lims)),
                       diag = list(continuous = wrap(diagfun, lims = lims))) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 15))
  return(score.plot)
}

lowerfun <- function(data, lims, mapping){
  ggplot(data = data, mapping = mapping) +
    geom_point(shape = 1, size = 2) +
    scale_x_continuous(limits = lims, breaks = breaks_pretty(n = 2)) +
    scale_y_continuous(limits = lims, breaks = breaks_pretty(n = 2))
}  

diagfun <- function(data, lims, mapping){
  ggplot(data = data, mapping = mapping) +
    geom_density(aes(fill = label), alpha = 0.5) +
    scale_x_continuous(limits = lims, breaks = breaks_pretty(n = 2))
}