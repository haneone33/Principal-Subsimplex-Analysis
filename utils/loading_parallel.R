loading_parallel <- function(v, var.names = NULL){
  
  v[v>0] = v[v>0]/sum(v[v>0])
  v[v<0] = -v[v<0]/sum(v[v<0])
  
  v1 = pmin(v, 0)
  v2 = pmax(v, 0)
  df = rbind(rep(0, length(v)), v2, v1)
  
  parallel_coord(df, X.label = as.factor(c(1,2,3)), var.names) +
    scale_color_manual(values = c('gray',gg_color_hue(2))) +
    geom_line(linewidth = 1) +
    scale_y_continuous(limits = c(-1,1)) +
    theme(title = element_text(hjust = 0.5, size = 10),
          axis.text.x = element_text(size = 15))
}
