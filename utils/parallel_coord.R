parallel_coord <- function(modes, X.label, var.names = colnames(modes)){

  modes = as.matrix(modes)
  n = nrow(modes)
  p = ncol(modes)
  mat = matrix(0, nrow = n, ncol = 2*p)
  mat[,(1:p)*2-1] = modes
  mat[,(1:p)*2] = modes

  if(is.null(var.names)){
    var.names = paste0('V',1:p)
  }

  df = as.data.frame(mat)
  df$label = X.label
  df$id = 1:nrow(df)

  df.melt = melt(df, id.vars = c('label','id'))
  df.melt$variable = as.numeric(as.factor(df.melt$variable))

  g <- ggplot(df.melt, aes(x=variable, y=value, group = id, col = label)) +
    geom_line(linewidth = 0.25) +
    theme_bw() +
    ylab(NULL) +
    scale_x_continuous(name = NULL, breaks = 2*(1:p)-0.5,
                       labels = var.names,
                       minor_breaks = 1:(ncol(modes)*2)) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_line(color = 'gray'),
          legend.position = 'none',
          plot.title = element_text(hjust = 0.5, size = 20))

  return(g)
}
