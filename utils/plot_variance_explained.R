plot_variance_explained <- function(rss, k = ncol(rss)){
  rss = rss[,1:k]/rowSums(rss)
  rss.cum = t(apply(rss, 1, cumsum))

  method = factor(rownames(rss), levels = rownames(rss))
  rss = as.data.frame(rss)
  rss.cum = as.data.frame(rss.cum)
  rss$Method = method
  rss.cum$Method = method

  g = ggplot(data = melt(rss, id.vars = 'Method'),
             aes(x = variable, y = value, group = Method, shape = Method, linetype = Method)) +
    theme_bw() +
    geom_line(col = 'black') +
    geom_point(size = 2) +
    geom_line(data = melt(rss.cum, id.vars = 'Method'), col = 'black') +
    geom_point(data = melt(rss.cum, id.vars = 'Method'), size = 2) +
    scale_y_continuous(limits = c(0,1)) +
    labs(x = 'Rank', y = 'Proportion of variance explained')

  return(g)
}
