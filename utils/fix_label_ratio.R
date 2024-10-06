fix_label_ratio <- function(g, label.ratio = 1){
  g1 = plot_grid(get_y_axis(g))
  g2 = g + theme(axis.text.y = element_blank(),
                 axis.ticks.y = element_line(),
                 plot.margin = unit(c(15,15,15,-15), 'points'))
  g3 = plot_grid(g1, g2, nrow = 1, rel_widths = c(label.ratio, 1),
                 align = 'h', axis = 'bt')
  return(g3)
}
