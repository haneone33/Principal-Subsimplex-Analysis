line_lowerfun <- function(data, mapping, col.point, col.line){
  ggplot(data = data, mapping = mapping) +
    geom_path(col = col.line, linewidth = 0.5) +
    geom_point(col = col.point, shape = 1, stroke = 1) 
}

jitter_diag <- function(data, mapping, col.point){
  set.seed(1)
  ggplot(data, mapping) +
    geom_jitter(aes(y = 0.5), height = 0.2, col = col.point, shape = 1, stroke = 1) +
    scale_y_continuous(limits = c(0,1)) +
    theme_void()
}

quantile_diag <- function(data, mapping, col.point, col.line){
  ggplot(data, mapping) +
    scale_y_continuous(trans = 'reverse') +
    geom_point(aes(y = Depth), col = col.point, shape = 1, stroke = 1) +
    geom_path(aes(y = Depth), col = col.line, linewidth = 0.5) +
    theme_void()
}

line_gpairs <- function(data, col.point, col.line, columns = 1:4){
  data = as.data.frame(data)
  warning('palo.df updated needed')
  data$Depth = diatom.df$Depth
  g = ggpairs(data, columns = columns,
              lower = list(continuous = wrap(line_lowerfun,
                                             col.point = col.point,
                                             col.line = col.line)),
              diag = list(continuous = wrap(quantile_diag,
                                            col.point = col.point,
                                            col.line = col.line)),
              upper = NULL) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 15))
}