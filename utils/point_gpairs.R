point_gpairs <- function(scores, X.label, k = min(ncol(scores)-1,5)){
  scores.df = as.data.frame(scores[,1:k])
  scores.df$label = X.label
  lims = get_lims(scores[,1:k])
  score.plot = ggpairs(scores.df, mapping = aes(col = label, shape = label, fill = label), columns = 1:k, upper = NULL,
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
  p <- ggplot(data = data, mapping = mapping) +
    scale_x_continuous(limits = lims, breaks = breaks_pretty(n = 2))

  var_name = ggplot2::as_label(mapping$x)
  group_var = ggplot2::as_label(mapping$fill)
  groups = unique(data[[group_var]])

  for(g in groups){
    data_g <- data[data[[group_var]] == g, , drop = FALSE]

    min_bw = (lims[2] - lims[1]) * 0.02
    p <- p + geom_density(data = data_g,
                          bw = base::max(bw.nrd0(data_g[[var_name]]), min_bw),
                          alpha = 0.5)
  }
  p
}
