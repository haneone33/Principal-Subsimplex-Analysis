score_plot_matrix <- function(scores, X.label, k = min(ncol(scores)-1,5),
                              plot.diag = (k!=2)){
  scores.df = as.data.frame(scores[,1:k])
  scores.df$label = X.label
  lims = plot.lims(scores[,1:k])
  score.plot = ggpairs(scores.df, aes(col = label), columns = 1:k, upper = NULL,
                       lower = list(continuous = wrap(lowerfun, lims = lims)),
                       diag = list(continuous = wrap(diagfun, lims = lims))) +
    theme_bw()
  if(!plot.diag){
    score.plot = gpairs_lower(score.plot)
  }
  return(score.plot)
}


lowerfun <- function(data, lims, mapping){
  ggplot(data = data, mapping = mapping) +
    geom_point(shape = 1, size = 2) +
    scale_x_continuous(limits = lims) +
    scale_y_continuous(limits = lims)
}  

diagfun <- function(data, lims, mapping){
  ggplot(data = data, mapping = mapping) +
    geom_density(aes(fill = label), alpha = 0.5) +
    xlim(lims)
}

compare_score_plot_matrix <- function(X.res,
                                      X.label = as.factor(rep(1,nrow(X.res$X))),
                                      k = min(ncol(X.res$X)-1,5),
                                      plot.diag = (k!=2),
                                      ternary = (ncol(X.res$X)==3),
                                      legend.plot = NULL,
                                      main = 'Score Plot Matrix'){
  if(ternary){
    g.data <- ggtern_colored(X.res$X, X.res$label) +
      ggtitle('Data') +
      theme(plot.title = element_text(hjust = 0.5, size = 15))
    g.data <- ggplot_gtable(ggtern::ggplot_build(g.data))
  }else if(!is.null(legend.plot)){
    g.data = legend.plot
  }
  else{
    g.data <- ggplot() + theme_void()
    g.data <- ggplot_gtable(ggtern::ggplot_build(g.data))
  }

  g.psas = score_plot_matrix(X.res$psas$scores, X.label, k = k, plot.diag = plot.diag) +
            ggtitle('PSA-S') +
            theme(plot.title = element_text(hjust = 0.5, size = 15))
  g.psao = score_plot_matrix(X.res$psao$scores, X.label, k = k, plot.diag = plot.diag) +
            ggtitle('PSA-O') +
            theme(plot.title = element_text(hjust = 0.5, size = 15))
  g.pca = score_plot_matrix(X.res$pca$scores, X.label, k = k, plot.diag = plot.diag) +
            ggtitle('PCA') +
            theme(plot.title = element_text(hjust = 0.5, size = 15))
  g.apca = score_plot_matrix(X.res$apca$scores, X.label, k = k, plot.diag = plot.diag) +
            ggtitle('Log-ratio PCA') +
            theme(plot.title = element_text(hjust = 0.5, size = 15))
  g.power = score_plot_matrix(X.res$power_pca$scores, X.label, k = k, plot.diag = plot.diag) +
              ggtitle('Power Transform PCA') +
              theme(plot.title = element_text(hjust = 0.5, size = 15))
  
  g = plot_grid(g.data,
                ggmatrix_gtable(g.psas),
                ggmatrix_gtable(g.psao),
                ggmatrix_gtable(g.pca), 
                ggmatrix_gtable(g.power),
                ggmatrix_gtable(g.apca), nrow = 2)
  
  if(!is.null(main)){
    g.title = ggdraw() + draw_label(main, size = 20)
    g = plot_grid(g.title, g, ncol=1, rel_heights=c(0.1, 1))
  }
  
  return(g)
}
