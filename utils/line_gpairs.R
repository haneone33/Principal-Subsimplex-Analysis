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
