line_lowerfun <- function(data, mapping, col.point, col.line, shape, size, limits){
  ggplot(data = data, mapping = mapping) +
    geom_path(col = col.line, linewidth = 0.5) +
    geom_point(col = col.point, stroke = 1, shape = shape, size = size) +
    scale_x_continuous(limits = limits) +
    scale_y_continuous(limits = limits)
}

jitter_diag <- function(data, mapping, col.point, shape, size){
  set.seed(1)
  ggplot(data, mapping) +
    geom_jitter(aes(y = 0.5), height = 0.2,
                col = col.point, stroke = 1, shape = shape, size = size) +
    scale_y_continuous(limits = c(0,1)) +
    theme_void()
}

quantile_diag <- function(data, mapping, col.point, col.line, shape, size, limits){
  ggplot(data, mapping) +
    scale_y_continuous(trans = 'reverse') +
    geom_path(aes(y = Depth), col = col.line, linewidth = 0.5) +
    geom_point(aes(y = Depth), col = col.point, stroke = 1, shape = shape, size = size) +
    theme_void() +
    scale_x_continuous(limits = limits)
}

line_gpairs <- function(data, col.point, col.line, shape, size, limits, columns = 1:4){
  data = as.data.frame(data)
  data$Depth = diatom.df$Depth
  g = ggpairs(data, columns = columns,
              lower = list(continuous = wrap(line_lowerfun,
                                             col.point = col.point,
                                             col.line = col.line,
                                             shape = shape,
                                             size = size,
                                             limits = limits)),
              diag = list(continuous = wrap(quantile_diag,
                                            col.point = col.point,
                                            col.line = col.line,
                                            shape = shape,
                                            size = size,
                                            limits = limits)),
              upper = NULL) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 15))
}
