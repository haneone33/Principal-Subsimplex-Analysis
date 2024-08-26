source('PSA_init.R')

library(xlsx)
library(viridis)
library(ggnewscale)

data.path = 'Diatom/Data/'
image.path = 'Diatom/Figures/'

################################################################################
## diatom preparation

## read data
diatom.df = read.csv(paste0(data.path, 'diatom.csv'))
diatom.X = diatom.df[,-1]
diatom.X = diatom.X[,names(sort(apply(diatom.X, 2, var), decreasing = T))] # sort by variance

## PSA application
# diatom.res = compare_analysis(diatom.X) # takes approximately an hour
# saveRDS(diatom.res, paste0(data.path, 'diatom_res.rds'))
diatom.res = readRDS(paste0(data.path, 'diatom_res.rds'))

################################################################################
## descriptive figures

## legend
g = ggplot(diatom.df, aes(x = Depth, y = Depth, col = Depth)) +
  geom_point() +
  theme_bw() +
  scale_color_viridis(name = 'Depth', trans = 'reverse', direction = -1,
                      guide = guide_colorbar(order = 1))
g = g + new_scale_color() +
  geom_point(data = g$data[36,], aes(color = '67.03'), shape = 1, stroke = 1) +
  scale_color_manual(values = c('67.03' = 'red'), name = 'Outlier')

diatom.legend = plot_grid(get_legend(g))

ggsave('1. diatom_legend.png', diatom.legend, path = image.path, width = 3, height = 6)

## color vector

diatom.col = ggplot_build(g)$data[[1]]$colour_ggnewscale_1
diatom.col.outlier = diatom.col
diatom.col.outlier[36] = 'red'

## parallel coordinate plot
g = parallel_coord(diatom.X, diatom.df$Depth) +
  scale_color_viridis(name = 'Depth', trans = 'reverse', direction = -1,
                      guide = guide_colorbar(order = 1)) +
  labs(x = '', y = 'Proportion') +
  theme(legend.position = 'right') +
  ggtitle('Diatom') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
g = g + new_scale_color() +
  geom_line(data = g$data[g$data$id == 36,], aes(color = '67.03'), linewidth = 0.6) +
  scale_color_manual(values = c('67.03' = 'red'), name = 'Outlier')
ggsave('1. Parallel_plot_raw.png', g, path = image.path, width = 12, height = 6)

## depth distribution
n = length(diatom.df$Depth)
diatom.depth.df = data.frame(Depth = diatom.df$Depth,
                           y = rev((1:n)/n-1/2/n))
g = ggplot(diatom.depth.df, aes(x = Depth, y = y, col = Depth)) +
  scale_x_continuous(trans = 'reverse') +
  geom_point(shape = 1) +
  labs(y='') +
  theme_bw() +
  scale_color_viridis(name = 'Depth', trans = 'reverse', direction = -1) +
  ggtitle('Diatom') +
  theme(plot.title = element_text(hjust = 0.5, size = 12))
g = g + new_scale_color() +
  geom_point(data = g$data[36,], aes(color = '67.03'), shape = 1) +
  scale_color_manual(values = c('67.03' = 'red'), name = 'Outlier')

ggsave('1. Depth_distribution.png', g, path = image.path, width = 8, height = 3)

## percent nonzero
df = data.frame(variable = colnames(diatom.X),
                percent.nonzero = colSums(diatom.X > 0)/nrow(diatom.X))
df$variable = factor(df$variable, levels = df$variable)
g = ggplot(df, aes(x = variable, y = percent.nonzero)) +
  geom_col() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(x = '', y = 'Proportion Nonzero') +
  geom_hline(yintercept = c(0.25,0.5,0.75), col = 'red') +
  ggtitle('Diatom') +
  theme(plot.title = element_text(hjust = 0.5, size = 20))

ggsave('1. Percent_nonzero.png', g, path = image.path, width = 8, height = 6)

################################################################################
## Compare 5 methods

## score scatter plot matrices
g.psas = line_gpairs(diatom.res$psas$scores, diatom.col.outlier, diatom.col) + ggtitle('PSA-S')
g.psao = line_gpairs(diatom.res$psao$scores, diatom.col.outlier, diatom.col) + ggtitle('PSA-O')
g.pca = line_gpairs(diatom.res$pca$scores, diatom.col.outlier, diatom.col) + ggtitle('PCA')
g.power_pca = line_gpairs(diatom.res$power_pca$scores, diatom.col.outlier, diatom.col) + ggtitle('Power Transform PCA')
g.apca = line_gpairs(diatom.res$apca$scores, diatom.col.outlier, diatom.col) + ggtitle('Log-ratio PCA')
g = plot_grid(diatom.legend, 
              ggmatrix_gtable(g.psas), 
              ggmatrix_gtable(g.psao),
              ggmatrix_gtable(g.pca), 
              ggmatrix_gtable(g.power_pca),
              ggmatrix_gtable(g.apca), nrow = 2)

ggsave('2. Score_plot_matrix.png', g, path = image.path, width = 18, height = 12)
ggsave('2-1. Score_psas.png', g.psas, path = image.path, width = 9, height = 9)
ggsave('2-2. Score_psao.png', g.psao, path = image.path, width = 9, height = 9)
ggsave('2-3. Score_pca.png', g.pca, path = image.path, width = 9, height = 9)
ggsave('2-4. Score_power_pca.png', g.power_pca, path = image.path, width = 9, height = 9)
ggsave('2-5. Score_apca.png', g.apca, path = image.path, width = 9, height = 9)

g = plot_grid(diatom.legend,
              getPlot(g.psas, 2, 1) + ggtitle('PSA-S'),
              getPlot(g.psao, 2, 1) + ggtitle('PSA-O'),
              getPlot(g.pca, 2, 1) + ggtitle('PCA'),
              getPlot(g.power_pca, 2, 1) + ggtitle('Power Transform PCA'),
              getPlot(g.apca, 2, 1) + ggtitle('Log-ratio PCA'),
              nrow = 2)
ggsave('2. Score_1vs2.png', g, path = image.path, width = 9, height = 6)

g = plot_grid(ggdraw() + draw_label('PSA-O'),
              plot_grid(getPlot(g.psao, 1, 1),
                        getPlot(g.psao, 2, 2),
                        getPlot(g.psao, 3, 3),
                        getPlot(g.psao, 4, 4),
                        diatom.legend,
                        nrow = 1, rel_widths = c(2,2,2,2,1)),
              ncol = 1, rel_heights = c(0.2, 3.3))
ggsave('2. Score_psao_diag.png', g, path = image.path, width = 18*0.8*9/8, height = 3.5)

## loading bar plot for PSA

g1 = loading_bar(diatom.res$psas$loadings[,1], max.k = 12) + ggtitle('Comp.1')
g2 = loading_bar(diatom.res$psas$loadings[,2], max.k = 12) + ggtitle('Comp.2')
g3 = loading_bar(diatom.res$psas$loadings[,3], max.k = 12) + ggtitle('Comp.3')
g4 = loading_bar(diatom.res$psas$loadings[,4], max.k = 12) + ggtitle('Comp.4')
g = plot_grid(ggdraw() + draw_label('Diatom PSA-S'),
              plot_grid(fix_label_ratio(g1, label.ratio = 1.5),
                        fix_label_ratio(g2, label.ratio = 1.5),
                        fix_label_ratio(g3, label.ratio = 1.5),
                        fix_label_ratio(g4, label.ratio = 1.5), nrow = 1),
              ncol = 1, rel_heights = c(1,9))
ggsave('3. loading_bar_psas.png', g, path = image.path, width = 16, height = 4)

g1 = loading_bar(diatom.res$psao$loadings[,1], max.k = 12) + ggtitle('Comp.1')
g2 = loading_bar(diatom.res$psao$loadings[,2], max.k = 12) + ggtitle('Comp.2')
g3 = loading_bar(diatom.res$psao$loadings[,3], max.k = 12) + ggtitle('Comp.3')
g4 = loading_bar(diatom.res$psao$loadings[,4], max.k = 12) + ggtitle('Comp.4')
g = plot_grid(ggdraw() + draw_label('Diatom PSA-O'),
              plot_grid(fix_label_ratio(g1, label.ratio = 1.5),
                        fix_label_ratio(g2, label.ratio = 1.5),
                        fix_label_ratio(g3, label.ratio = 1.5),
                        fix_label_ratio(g4, label.ratio = 1.5), nrow = 1),
              ncol = 1, rel_heights = c(1,9))
ggsave('3. loading_bar_psao.png', g, path = image.path, width = 16, height = 4)

## loading bar plot for 5 methods
diatom.loading <- function(loadings){
  g = plot_grid(fix_label_ratio(loading_bar(loadings[,1], max.k = 12), label.ratio = 1.5),
                fix_label_ratio(loading_bar(loadings[,2], max.k = 12), label.ratio = 1.5),
                fix_label_ratio(loading_bar(loadings[,3], max.k = 12), label.ratio = 1.5),
                fix_label_ratio(loading_bar(loadings[,4], max.k = 12), label.ratio = 1.5), nrow = 1)
  return(g)
}

title.col = plot_grid(ggdraw(), ggdraw() + draw_label('Comp.1'),
                      ggdraw(), ggdraw() + draw_label('Comp.2'),
                      ggdraw(), ggdraw() + draw_label('Comp.3'),
                      ggdraw(), ggdraw() + draw_label('Comp.4'),
                      nrow = 1, rel_widths = c(1.5, 1, 1.5, 1, 1.5, 1, 1.5, 1))
title.row = plot_grid(ggdraw() + draw_label('PSA-S'),
                      ggdraw() + draw_label('PSA-O'),
                      ggdraw() + draw_label('PCA'),
                      ggdraw() + draw_label('Power Transform PCA'),
                      ggdraw() + draw_label('Log-ratio PCA'), ncol = 1)
g.main = plot_grid(diatom.loading(diatom.res$psas$loadings),
                   diatom.loading(diatom.res$psao$loadings),
                   diatom.loading(diatom.res$pca$loadings),
                   diatom.loading(diatom.res$power_pca$loadings),
                   diatom.loading(diatom.res$apca$loadings), ncol = 1)

g = plot_grid(ggdraw(), title.col,
              title.row, g.main, nrow = 2,
              rel_heights = c(0.2, 9), rel_widths = c(0.8,4))

ggsave('3. loading_bar_all.png', g, path = image.path, width = 16, height = 16)

## loading parallel coordinate plot

g.axis = loading_parallel(diatom.res$psas$loadings[,1]) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.ticks.x = element_blank())

g = plot_grid(ggdraw() + draw_label('Diatom PSA-S', size = 20),
              plot_grid(plot_grid(ggdraw() + draw_label('Comp.1'), ggdraw() + draw_label('Comp.2'),
                                  ggdraw() + draw_label('Comp.3'), ggdraw() + draw_label('Comp.4'),
                                  ggdraw(),
                                  ncol = 1, rel_heights = c(1,1,1,1,3)),
                        plot_grid(loading_parallel(diatom.res$psas$loadings[,1]) +
                                    theme(axis.text.x = element_blank()),
                                  loading_parallel(diatom.res$psas$loadings[,2]) +
                                    theme(axis.text.x = element_blank()),
                                  loading_parallel(diatom.res$psas$loadings[,3]) +
                                    theme(axis.text.x = element_blank()),
                                  loading_parallel(diatom.res$psas$loadings[,4]) +
                                    theme(axis.text.x = element_blank()),
                                  ggdraw(get_x_axis(g.axis)),
                                  ncol = 1, rel_heights = c(1,1,1,1,3), align = 'v', axis = 'lr'),
                        nrow = 1, rel_widths = c(1,9)),
              ncol = 1, rel_heights = c(1,9))

ggsave('4. loading_parallel_psas.png', g, path = image.path, width = 12, height = 12)

g = plot_grid(ggdraw() + draw_label('Diatom PSA-O', size = 20),
              plot_grid(plot_grid(ggdraw() + draw_label('Comp.1'), ggdraw() + draw_label('Comp.2'),
                                  ggdraw() + draw_label('Comp.3'), ggdraw() + draw_label('Comp.4'),
                                  ggdraw(),
                                  ncol = 1, rel_heights = c(1,1,1,1,3)),
                        plot_grid(loading_parallel(diatom.res$psao$loadings[,1]) +
                                    theme(axis.text.x = element_blank()),
                                  loading_parallel(diatom.res$psao$loadings[,2]) +
                                    theme(axis.text.x = element_blank()),
                                  loading_parallel(diatom.res$psao$loadings[,3]) +
                                    theme(axis.text.x = element_blank()),
                                  loading_parallel(diatom.res$psao$loadings[,4]) +
                                    theme(axis.text.x = element_blank()),
                                  ggdraw(get_x_axis(g.axis)),
                                  ncol = 1, rel_heights = c(1,1,1,1,3), align = 'v', axis = 'lr'),
                        nrow = 1, rel_widths = c(1,9)),
              ncol = 1, rel_heights = c(1,9))

ggsave('4. loading_parallel_psao.png', g, path = image.path, width = 12, height = 12)

## ternary plot
df.psas = as.data.frame(diatom.res$psas$pts$`r=3`)
colnames(df.psas) = c('V1','V2','V3')
df.psas$Depth = diatom.df$Depth
psas.tern = ggtern(df.psas, aes(x=V1, y=V2, z=V3,
                                col = sort(diatom.df$Depth, index.return = T)[[2]])) +
  geom_point(shape = 1, stroke = 1) +
  geom_path(linewidth = 0.5) +
  theme_bw() +
  theme(legend.position = 'none') +
  scale_color_viridis() +
  ggtitle('PSA-S') +
  theme(plot.title = element_text(hjust = 0.5, size = 20)) +
  geom_point(data = data.frame(df.psas)[36,], col = 'red', shape = 1, stroke = 1)
ggsave('3-1. ternary_psas.png', psas.tern, path = image.path, width = 6, height = 6)

df.psao = as.data.frame(diatom.res$psao$pts$`r=3`)
colnames(df.psao) = c('V1','V2','V3')
df.psao$Depth = diatom.df$Depth
psao.tern = ggtern(df.psao, aes(x=V1, y=V2, z=V3,
                                col = sort(diatom.df$Depth, index.return = T)[[2]])) +
  geom_point(shape = 1, stroke = 1) +
  geom_path(linewidth = 0.5) +
  theme_bw() +
  theme(legend.position = 'none') +
  scale_color_viridis() +
  ggtitle('PSA-O') +
  theme(plot.title = element_text(hjust = 0.5, size = 20)) +
  geom_point(data = data.frame(df.psao)[36,], col = 'red', shape = 1, stroke = 1)
ggsave('3-2. ternary_psao.png', psao.tern, path = image.path, width = 6, height = 6)

## ternary vertices

g1 = loading_bar(diatom.res$psas$vertices$`r=3`[1,], max.k = 12) +
  ggtitle('V1') +
  scale_fill_manual(values = 'darkgray')
g2 = loading_bar(diatom.res$psas$vertices$`r=3`[2,], max.k = 12) +
  ggtitle('V2') +
  scale_fill_manual(values = 'darkgray')
g3 = loading_bar(diatom.res$psas$vertices$`r=3`[3,], max.k = 12) +
  ggtitle('V3') +
  scale_fill_manual(values = 'darkgray')
g = plot_grid(ggdraw() + draw_label('Diatom PSA-S Vertices'),
              plot_grid(fix_label_ratio(g1, label.ratio = 1.5),
                        fix_label_ratio(g2, label.ratio = 1.5),
                        fix_label_ratio(g3, label.ratio = 1.5),
                        nrow = 1),
              ncol = 1, rel_heights = c(1,9))
ggsave('3-3. ternary_vertex_psas.png', g, path = image.path, width = 12, height = 4)          

g1 = loading_bar(diatom.res$psao$vertices$`r=3`[1,], max.k = 12) +
  ggtitle('V1') +
  scale_fill_manual(values = 'darkgray')
g2 = loading_bar(diatom.res$psao$vertices$`r=3`[2,], max.k = 12) +
  ggtitle('V2') +
  scale_fill_manual(values = 'darkgray')
g3 = loading_bar(diatom.res$psao$vertices$`r=3`[3,], max.k = 12) +
  ggtitle('V3') +
  scale_fill_manual(values = 'darkgray')
g = plot_grid(ggdraw() + draw_label('Diatom PSA-O Vertices'),
              plot_grid(fix_label_ratio(g1, label.ratio = 1.5),
                        fix_label_ratio(g2, label.ratio = 1.5),
                        fix_label_ratio(g3, label.ratio = 1.5),
                        nrow = 1),
              ncol = 1, rel_heights = c(1,9))
ggsave('3-4. ternary_vertex_psao.png', g, path = image.path, width = 12, height = 4)     

