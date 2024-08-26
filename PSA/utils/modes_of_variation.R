get_loading <- function(dend.input){
  max_r = nrow(dend.input$merges)
  modes = vector(mode = 'list', length = max_r)
  names(modes) = rownames(dend.input$merges)
  
  for(r in (max_r-1):1){
    v1 = dend.input$vertices[[paste0('r=',r+1)]][dend.input$merges[paste0('r=',r),'v1'],]
    v2 = dend.input$vertices[[paste0('r=',r+1)]][dend.input$merges[paste0('r=',r),'v2'],]
    modes[[paste0('r=',r)]] = list(v1 = v1, v2 = v2)
  }
  
  return(modes)
}

legend_mode_vertex <- function(){
  g = parallel_coord(matrix(c(1,-1),ncol = 1),
                     X.label = factor(c('V2','-V1'), levels = c('V2','-V1'))) +
    geom_line(size = 1) +
    scale_color_manual(values = gg_color_hue(2)) +
    scale_color_discrete(name = 'Vertex') +
    theme(legend.position = 'right',
          legend.text = element_text(size = 20),
          legend.title =element_text(size = 20))
  plot_grid(get_legend(g))
}

loading_bar <- function(v, max.k = 20){
  
  v[v>0] = v[v>0]/sum(v[v>0])
  v[v<0] = -v[v<0]/sum(v[v<0])
  v = v[abs(v)>0]
  
  if(length(v) > max.k){
    threshold = sort(abs(v), decreasing = T)[max.k]
    v = v[abs(v) >= threshold]
  }
  
  v = sort(v, decreasing = T)
  df = data.frame(variable = factor(names(v), levels = names(v)),
                  value = v,
                  vertex = factor(v>0, levels = c(T,F)))
  rownames(df) = NULL
  
  g = ggplot(data = df, aes(x = variable, y = value, fill = vertex)) +
    geom_bar(stat = 'identity') +
    coord_flip() +
    theme_bw() +
    scale_x_discrete(limits = rev(levels(df$variable))) +
    theme(axis.ticks.y = element_blank()) +
    labs(x = '', y = '') +
    theme(legend.position = 'none') +
    theme(plot.title = element_text(hjust = 0.5, size = 14)) +
    theme(plot.margin = unit(c(0.5,0.5,0.5,-0.5), "cm"))
  
  return(g)
}

fix_label_ratio <- function(g, label.ratio = 1){
  g1 = plot_grid(get_y_axis(g))
  g2 = g + theme(axis.text.y = element_blank(),
                 axis.ticks.y = element_line())
  g3 = plot_grid(g1, g2, nrow = 1, rel_widths = c(label.ratio, 1),
                align = 'h', axis = 'bt')
  return(g3)
}

loading_parallel <- function(v){
  
  v[v>0] = v[v>0]/sum(v[v>0])
  v[v<0] = -v[v<0]/sum(v[v<0])
  
  v1 = pmin(v, 0)
  v2 = pmax(v, 0)
  df = rbind(rep(0, length(v)), v2, v1)
  
  parallel_coord(df, X.label = as.factor(c(1,2,3))) +
    scale_color_manual(values = c('gray',gg_color_hue(2))) +
    geom_line(linewidth = 1) +
    scale_y_continuous(limits = c(-1,1)) +
    theme(title = element_text(hjust = 0.5, size = 10),
          axis.text.x = element_text(size = 15))
}
