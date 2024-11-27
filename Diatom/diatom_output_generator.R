devtools::install_github('haneone33/Principal-Subsimplex-Analysis', subdir = 'PSA')
library(PSA)

require(compositions)
require(ggplot2)
require(ggtern)
require(GGally)
require(cowplot) # plot_grid
require(reshape) # parallel coordinate plot
library(viridis)
library(ggnewscale)

invisible(lapply(list.files('utils', pattern = '.R', full.names = T), source))

data.path = 'Diatom/Data/'
image.path = 'Diatom/Figures/'
dir.create(data.path, showWarnings = F)
dir.create(image.path, showWarnings = F)

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

diatom.info = read.csv(paste0(data.path,'diatom codes.csv'))
diatom.info$class[diatom.info$class == ''] = 'others'
diatom.info$class = factor(diatom.info$class, levels = c('warm water','open ocean','sea ice','others'))
diatom.info$color = c('hotpink','green3','dodgerblue','black')[diatom.info$class]
rownames(diatom.info) = colnames(diatom.df)[-1]

################################################################################
## descriptive figures

## parallel coordinate plot
g = parallel_coord(diatom.X, diatom.df$Depth, diatom.info[as.character(colnames(diatom.X)),'code']) +
  scale_color_viridis(name = 'Depth', trans = 'reverse', direction = -1,
                      guide = guide_colorbar(order = 1)) +
  labs(x = '', y = 'Proportion') +
  theme(legend.position = 'right', legend.key.size = unit(10, 'points'), legend.title = element_text(size = 10)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.key.size = unit(12, 'points'),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.spacing.y = unit(1,'points')) +
  theme(axis.text.x = element_text(color = diatom.info[as.character(colnames(diatom.X)),'color']))
g = g + new_scale_color() +
  geom_line(data = g$data[g$data$id %in% c(37:39),], aes(color = as.factor(label)), linewidth = 0.6) +
  scale_color_manual(values = c('69.74'='orangered','69.84'='darkorange1','70.04'='orange'),
                     name = 'KM3', guide = guide_legend(order = 2))
g = g + new_scale_color() +
  geom_line(data = g$data[g$data$id == 36,], aes(color = as.factor(label)), linewidth = 0.6) +
  scale_color_manual(values = c('67.03'='red'), name = 'Post-KM3', guide = guide_legend(order = 3))

diatom.legend.line = g
ggsave('1. Parallel_plot_raw.jpeg', g, path = image.path, width = 12, height = 6)

## depth distribution
n = length(diatom.df$Depth)
diatom.depth.df = data.frame(Depth = diatom.df$Depth,
                             y = rev((1:n)/n-1/2/n))
g = ggplot(diatom.depth.df, aes(x = Depth, y = y, col = Depth)) +
  scale_x_continuous(trans = 'reverse') +
  geom_point(shape = 1.5, size = 2, stroke = 1) +
  labs(y='') +
  theme_bw() +
  scale_color_viridis(name = 'Depth', trans = 'reverse', direction = -1,
                      guide = guide_colorbar(order = 1))
g = g + new_scale_color() +
  geom_point(data = g$data[37:39,], aes(color = as.factor(Depth)), shape = 16, size = 3, stroke = 1) +
  scale_color_manual(values = c('69.74'='orangered','69.84'='darkorange1','70.04'='orange'),
                     name = 'KM3', guide = guide_legend(order = 2))
g = g + new_scale_color() +
  geom_point(data = g$data[36,], aes(color = as.factor(Depth)), shape = 16, size = 3, stroke = 1) +
  scale_color_manual(values = c('67.03'='red'), name = 'Post-KM3', guide = guide_legend(order = 3))
g = g + theme(legend.key.size = unit(8, 'points'),
              legend.text = element_text(size = 6),
              legend.title = element_text(size = 8),
              legend.spacing.y = unit(1,'points'),
              plot.margin = unit(rep(15,4),'points'),
              axis.title.x = element_text(size = 8))

diatom.legend.point = g
ggsave('1. Depth_distribution.jpeg', g, path = image.path, width = 8, height = 3)

## color, shape, size vectors
diatom.col = ggplot_build(diatom.legend.point)$data[[1]]$colour_ggnewscale_1
diatom.col.outlier = diatom.col
diatom.col.outlier[36:39] = c('red','orangered','darkorange1','orange')

diatom.shape = c(rep(1,35),rep(16,4),rep(1,32))
diatom.size = c(rep(1,35),rep(2,4),rep(1,32))

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
  theme(axis.text.x = element_text(color = diatom.info[as.character(colnames(diatom.X)),'color'])) +
  scale_x_discrete(labels = diatom.info[as.character(colnames(diatom.X)),'code'])

ggsave('1. Percent_nonzero.jpeg', g, path = image.path, width = 8, height = 4)

################################################################################
## Compare 5 methods

## score scatter plot matrices
matrix.legend = plot_grid(get_legend(diatom.legend.point + theme(legend.key.size = unit(16, 'points'),
                                                                 legend.text = element_text(size = 14),
                                                                 legend.title = element_text(size = 16),
                                                                 legend.spacing.y = unit(1,'points'))))

g.psas = line_gpairs(diatom.res$psas$scores, col.point = diatom.col.outlier, col.line = diatom.col,
                     shape = diatom.shape, size = diatom.size)
g.psao = line_gpairs(diatom.res$psao$scores, col.point = diatom.col.outlier, col.line = diatom.col,
                     shape = diatom.shape, size = diatom.size)
g.pca = line_gpairs(diatom.res$pca$scores, col.point = diatom.col.outlier, col.line = diatom.col,
                    shape = diatom.shape, size = diatom.size)
g.power_pca = line_gpairs(diatom.res$power_pca$scores, col.point = diatom.col.outlier, col.line = diatom.col,
                          shape = diatom.shape, size = diatom.size)
g.apca = line_gpairs(diatom.res$apca$scores, col.point = diatom.col.outlier, col.line = diatom.col,
                     shape = diatom.shape, size = diatom.size)

g = plot_grid(matrix.legend,
              ggmatrix_gtable(g.psas), 
              ggmatrix_gtable(g.psao),
              ggmatrix_gtable(g.pca), 
              ggmatrix_gtable(g.power_pca),
              ggmatrix_gtable(g.apca), nrow = 2)

ggsave('2. Score_plot_matrix.jpeg', g, path = image.path, width = 18, height = 12)


g = plot_grid(get_legend(diatom.legend.point + theme(legend.key.size = unit(14, 'points'),
                                                     legend.text = element_text(size = 12),
                                                     legend.title = element_text(size = 14),
                                                     legend.spacing.y = unit(1,'points'))),
              getPlot(g.psas, 2, 1) + ggtitle('PSA-S'),
              getPlot(g.psao, 2, 1) + ggtitle('PSA-O'),
              getPlot(g.pca, 2, 1) + ggtitle('PCA'),
              getPlot(g.power_pca, 2, 1) + ggtitle('Power Transform PCA'),
              getPlot(g.apca, 2, 1) + ggtitle('Log-ratio PCA'),
              nrow = 2)
ggsave('2. Score_1vs2.jpeg', g, path = image.path, width = 10.5, height = 7)

g = plot_grid(getPlot(g.psao, 1, 1),
              getPlot(g.psao, 2, 2),
              getPlot(g.psao, 3, 3),
              getPlot(g.psao, 4, 4),
              get_legend(diatom.legend.point + theme(legend.key.size = unit(14, 'points'),
                                                     legend.text = element_text(size = 12),
                                                     legend.title = element_text(size = 14),
                                                     legend.spacing.y = unit(1,'points'))),
              nrow = 1, rel_widths = c(2,2,2,2,1))
ggsave('2. Score_psao_diag.jpeg', g, path = image.path, width = 14.5, height = 3)


g.psas = plot_grid(ggmatrix_gtable(g.psas), matrix.legend, nrow = 1, rel_widths = c(5,1))
g.psao = plot_grid(ggmatrix_gtable(g.psao), matrix.legend, nrow = 1, rel_widths = c(5,1))
g.pca = plot_grid(ggmatrix_gtable(g.pca), matrix.legend, nrow = 1, rel_widths = c(5,1))
g.power_pca = plot_grid(ggmatrix_gtable(g.power_pca), matrix.legend, nrow = 1, rel_widths = c(5,1))
g.apca = plot_grid(ggmatrix_gtable(g.apca), matrix.legend, nrow = 1, rel_widths = c(5,1))

ggsave('2-1. Score_psas.jpeg', g.psas, path = image.path, width = 9*1.2, height = 9)
ggsave('2-2. Score_psao.jpeg', g.psao, path = image.path, width = 9*1.2, height = 9)
ggsave('2-3. Score_pca.jpeg', g.pca, path = image.path, width = 9*1.2, height = 9)
ggsave('2-4. Score_power_pca.jpeg', g.power_pca, path = image.path, width = 9*1.2, height = 9)
ggsave('2-5. Score_apca.jpeg', g.apca, path = image.path, width = 9*1.2, height = 9)



## loading bar plot for PSA

g1 = loading_bar(diatom.res$psas$loadings[,1], max.k = 12) + ggtitle('Comp.1')
g2 = loading_bar(diatom.res$psas$loadings[,2], max.k = 12) + ggtitle('Comp.2')
g3 = loading_bar(diatom.res$psas$loadings[,3], max.k = 12) + ggtitle('Comp.3')
g4 = loading_bar(diatom.res$psas$loadings[,4], max.k = 12) + ggtitle('Comp.4')
g1 = g1 + scale_x_discrete(labels = diatom.info[as.character(g1$data$variable),'code']) +
  theme(axis.text.y = element_text(color = diatom.info[as.character(g1$data$variable),'color'], size = 12))
g2 = g2 + scale_x_discrete(labels = diatom.info[as.character(g2$data$variable),'code']) +
  theme(axis.text.y = element_text(color = diatom.info[as.character(g2$data$variable),'color'], size = 12))
g3 = g3 + scale_x_discrete(labels = diatom.info[as.character(g3$data$variable),'code']) +
  theme(axis.text.y = element_text(color = diatom.info[as.character(g3$data$variable),'color'], size = 12))
g4 = g4 + scale_x_discrete(labels = diatom.info[as.character(g4$data$variable),'code']) +
  theme(axis.text.y = element_text(color = diatom.info[as.character(g4$data$variable),'color'], size = 12))
g = plot_grid(fix_label_ratio(g1, label.ratio = 1),
              fix_label_ratio(g2, label.ratio = 1),
              fix_label_ratio(g3, label.ratio = 1),
              fix_label_ratio(g4, label.ratio = 1), nrow = 1)
ggsave('3. loading_bar_psas.jpeg', g, path = image.path, width = 16, height = 4)

g1 = loading_bar(diatom.res$psao$loadings[,1], max.k = 12) + ggtitle('Comp.1')
g2 = loading_bar(diatom.res$psao$loadings[,2], max.k = 12) + ggtitle('Comp.2')
g3 = loading_bar(diatom.res$psao$loadings[,3], max.k = 12) + ggtitle('Comp.3')
g4 = loading_bar(diatom.res$psao$loadings[,4], max.k = 12) + ggtitle('Comp.4')
g1 = g1 + scale_x_discrete(labels = diatom.info[as.character(g1$data$variable),'code']) +
  theme(axis.text.y = element_text(color = diatom.info[as.character(g1$data$variable),'color'], size = 12))
g2 = g2 + scale_x_discrete(labels = diatom.info[as.character(g2$data$variable),'code']) +
  theme(axis.text.y = element_text(color = diatom.info[as.character(g2$data$variable),'color'], size = 12))
g3 = g3 + scale_x_discrete(labels = diatom.info[as.character(g3$data$variable),'code']) +
  theme(axis.text.y = element_text(color = diatom.info[as.character(g3$data$variable),'color'], size = 12))
g4 = g4 + scale_x_discrete(labels = diatom.info[as.character(g4$data$variable),'code']) +
  theme(axis.text.y = element_text(color = diatom.info[as.character(g4$data$variable),'color'], size = 12))
g = plot_grid(fix_label_ratio(g1, label.ratio = 1),
              fix_label_ratio(g2, label.ratio = 1),
              fix_label_ratio(g3, label.ratio = 1),
              fix_label_ratio(g4, label.ratio = 1), nrow = 1)
ggsave('3. loading_bar_psao.jpeg', g, path = image.path, width = 16, height = 4)

## loading bar plot for 5 methods
diatom.loading <- function(loadings){
  g1 = loading_bar(loadings[,1], max.k = 12)
  g2 = loading_bar(loadings[,2], max.k = 12)
  g3 = loading_bar(loadings[,3], max.k = 12)
  g4 = loading_bar(loadings[,4], max.k = 12)
  g1 = g1 + scale_x_discrete(labels = diatom.info[as.character(g1$data$variable),'code']) +
    theme(axis.text.y = element_text(color = diatom.info[as.character(g1$data$variable),'color'], size = 12))
  g2 = g2 + scale_x_discrete(labels = diatom.info[as.character(g2$data$variable),'code']) +
    theme(axis.text.y = element_text(color = diatom.info[as.character(g2$data$variable),'color'], size = 12))
  g3 = g3 + scale_x_discrete(labels = diatom.info[as.character(g3$data$variable),'code']) +
    theme(axis.text.y = element_text(color = diatom.info[as.character(g3$data$variable),'color'], size = 12))
  g4 = g4 + scale_x_discrete(labels = diatom.info[as.character(g4$data$variable),'code']) +
    theme(axis.text.y = element_text(color = diatom.info[as.character(g4$data$variable),'color'], size = 12))
  g = plot_grid(fix_label_ratio(g1, label.ratio = 1.5),
                fix_label_ratio(g2, label.ratio = 1.5),
                fix_label_ratio(g3, label.ratio = 1.5),
                fix_label_ratio(g4, label.ratio = 1.5), nrow = 1)
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
ggsave('3. loading_bar_all.jpeg', g, path = image.path, width = 16, height = 16)

## loading parallel coordinate plot

g.axis = loading_parallel(diatom.res$psas$loadings[,1],
                          diatom.info[as.character(colnames(diatom.X)),'code']) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
        axis.ticks.x = element_blank()) +
  theme(axis.text.x = element_text(color = diatom.info[as.character(colnames(diatom.X)),'color']))

g = plot_grid(plot_grid(ggdraw() + draw_label('Comp.1'), ggdraw() + draw_label('Comp.2'),
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
                        ncol = 1, rel_heights = c(1,1,1,1,2), align = 'v', axis = 'lr'),
              nrow = 1, rel_widths = c(1,9))

ggsave('4. loading_parallel_psas.jpeg', g, path = image.path, width = 12, height = 7)

g = plot_grid(plot_grid(ggdraw() + draw_label('Comp.1'), ggdraw() + draw_label('Comp.2'),
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
                        ncol = 1, rel_heights = c(1,1,1,1,2), align = 'v', axis = 'lr'),
              nrow = 1, rel_widths = c(1,9))

ggsave('4. loading_parallel_psao.jpeg', g, path = image.path, width = 12, height = 7)

## ternary plot
df.psas = as.data.frame(diatom.res$psas$pts$`r=3`)
colnames(df.psas) = c('V1','V2','V3')
df.psas$Depth = diatom.df$Depth
psas.tern = ggtern(df.psas, aes(x=V1, y=V3, z=V2)) +
  geom_path(col = diatom.col, linewidth = 0.5) +
  geom_point(col = diatom.col.outlier, shape = diatom.shape, size = diatom.size*4/3, stroke = 1) +
  theme_bw() +
  theme(legend.position = 'none') +
  scale_color_viridis() +
  geom_point(data = data.frame(df.psas)[36:39,], col = c('red','orangered','darkorange1','orange'), shape = 1, stroke = 1)

psas.tern = plot_grid(ggtern::ggplot_gtable(ggtern::ggplot_build(psas.tern)),
                      get_legend(diatom.legend.point + theme(legend.key.size = unit(16, 'points'),
                                                             legend.text = element_text(size = 14),
                                                             legend.title = element_text(size = 16))),
                      rel_widths = c(3,1))
ggsave('3-1. ternary_psas.jpeg', psas.tern, path = image.path, width = 8, height = 5)

df.psao = as.data.frame(diatom.res$psao$pts$`r=3`)
colnames(df.psao) = c('V1','V2','V3')
df.psao$Depth = diatom.df$Depth
psao.tern = ggtern(df.psao, aes(x=V1, y=V3, z=V2)) +
  geom_path(col = diatom.col, linewidth = 0.5) +
  geom_point(col = diatom.col.outlier, shape = diatom.shape, size = diatom.size*4/3, stroke = 1) +
  theme_bw() +
  theme(legend.position = 'none') +
  scale_color_viridis() +
  geom_point(data = data.frame(df.psao)[36:39,], col = c('red','orangered','darkorange1','orange'), shape = 1, stroke = 1)

psao.tern = plot_grid(ggtern::ggplot_gtable(ggtern::ggplot_build(psao.tern)),
                      get_legend(diatom.legend.point + theme(legend.key.size = unit(16, 'points'),
                                                             legend.text = element_text(size = 14),
                                                             legend.title = element_text(size = 16))),
                      rel_widths = c(3,1))
ggsave('3-2. ternary_psao.jpeg', psao.tern, path = image.path, width = 8, height = 5)

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
g1 = g1 + scale_x_discrete(labels = diatom.info[as.character(g1$data$variable),'code']) +
  theme(axis.text.y = element_text(color = diatom.info[as.character(g1$data$variable),'color'], size = 11))
g2 = g2 + scale_x_discrete(labels = diatom.info[as.character(g2$data$variable),'code']) +
  theme(axis.text.y = element_text(color = diatom.info[as.character(g2$data$variable),'color'], size = 11))
g3 = g3 + scale_x_discrete(labels = diatom.info[as.character(g3$data$variable),'code']) +
  theme(axis.text.y = element_text(color = diatom.info[as.character(g3$data$variable),'color'], size = 11))

g = plot_grid(fix_label_ratio(g1, label.ratio = 1.2),
              fix_label_ratio(g2, label.ratio = 1.2),
              fix_label_ratio(g3, label.ratio = 1.2),
              nrow = 1)
ggsave('3-3. ternary_vertex_psas.jpeg', g, path = image.path, width = 12, height = 4)          


g1 = loading_bar(diatom.res$psao$vertices$`r=3`[1,], max.k = 12) +
  ggtitle('V1') +
  scale_fill_manual(values = 'darkgray')
g2 = loading_bar(diatom.res$psao$vertices$`r=3`[2,], max.k = 12) +
  ggtitle('V2') +
  scale_fill_manual(values = 'darkgray')
g3 = loading_bar(diatom.res$psao$vertices$`r=3`[3,], max.k = 12) +
  ggtitle('V3') +
  scale_fill_manual(values = 'darkgray')
g1 = g1 + scale_x_discrete(labels = diatom.info[as.character(g1$data$variable),'code']) +
  theme(axis.text.y = element_text(color = diatom.info[as.character(g1$data$variable),'color'], size = 11))
g2 = g2 + scale_x_discrete(labels = diatom.info[as.character(g2$data$variable),'code']) +
  theme(axis.text.y = element_text(color = diatom.info[as.character(g2$data$variable),'color'], size = 11))
g3 = g3 + scale_x_discrete(labels = diatom.info[as.character(g3$data$variable),'code']) +
  theme(axis.text.y = element_text(color = diatom.info[as.character(g3$data$variable),'color'], size = 11))

g = plot_grid(fix_label_ratio(g1, label.ratio = 1.2),
              fix_label_ratio(g2, label.ratio = 1.2),
              fix_label_ratio(g3, label.ratio = 1.2),
              nrow = 1)
ggsave('3-4. ternary_vertex_psao.jpeg', g, path = image.path, width = 12, height = 4)     

