require(devtools)
devtools::install_github('haneone33/Principal-Subsimplex-Analysis', subdir = 'PSA')
library(PSA)

require(compositions)
require(ggplot2)
require(ggtern)
require(GGally)
require(cowplot) # plot_grid
require(reshape) # parallel coordinate plot
library(viridis)
library(scales)
library(ggnewscale)

invisible(lapply(list.files('utils', pattern = '.R', full.names = T), source))

data.path = 'Diatom/Data/'
image.path = 'Diatom/Figures_revised/'
dir.create(data.path, showWarnings = T)
dir.create(image.path, showWarnings = T)

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

## heatmap
diatom.X.log = log10(replace_zero(diatom.X))
diatom.ticks = as.numeric(round(quantile(1:nrow(diatom.df))))

heatmap.df = as.data.frame(diatom.X.log) %>%
  cbind(data.frame(Depth = diatom.df$Depth)) %>%
  mutate(Depth_idx = rank(-Depth)) %>%
  select(-Depth) %>%
  melt(id.vars = 'Depth_idx') %>%
  mutate(variable = factor(variable, levels = colnames(diatom.X)))

g2 = ggplot(heatmap.df, aes(x = variable, y = Depth_idx, fill = value)) +
  theme_minimal() +
  geom_tile() +
  scale_fill_gradientn(colors = c('white','black','black'),
                       values = rescale(c(max(diatom.X.log), -1.962437, min(diatom.X.log))),
                       name = 'log10(Proportion)') +
  scale_x_discrete(labels = diatom.info[as.character(colnames(diatom.X)),'code']) +
  scale_y_continuous(breaks = diatom.ticks, expand = c(0,0,0,0),
                     position = 'right') +
  labs(x = NULL, y = 'Depth (Index)') +
  theme(axis.text.x = element_text(color = diatom.info[as.character(colnames(diatom.X)),'color'],
                                   angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.title.align = 0.5)

g1 = data.frame(variable = 'Depth',
                Depth = diatom.df$Depth) %>%
  ggplot(aes(x = variable, y = rank(Depth), fill = Depth)) +
  theme_minimal() +
  geom_tile() +
  scale_fill_gradientn(colors = c('blue','magenta','red','orange','green'),
                       values = rescale(c(56.53,65.04,67.03,69.74,83.99), to = c(0,1)),
                       guide = guide_colorbar(direction = "vertical", reverse = TRUE)) +
  scale_y_continuous(trans = 'reverse',
                     labels = diatom.df$Depth[diatom.ticks], breaks = diatom.ticks,
                     expand = c(0,0,0,0)) +
  labs(x = NULL, y = 'Depth') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.position = 'none')

g = plot_grid(g1, g2, nrow = 1, rel_widths = c(1,10), align = 'h', axis = 'tb')
ggsave('1. diatom_heatmap.jpeg', g, path = image.path, width = 12, height = 4)


## mean proportion versus log scale variance
g = data.frame(mean_prop = colMeans(diatom.X),
               sd_log10 = apply(diatom.X.log, 2, sd),
               label = diatom.info$class) %>%
  ggplot(aes(x = mean_prop, y = sd_log10, col = label, shape = label)) +
  theme_bw() +
  geom_point() +
  scale_color_manual(values = c('hotpink','green3','dodgerblue','black'),
                     name = 'Type') +
  scale_shape_manual(values = c(16,16,16,1),
                     name = 'Type') +
  labs(x = 'Mean proportion', y = 'Standard deviation (log10)')
ggsave('1. diatom_mean_proportion_versus_variance.jpeg', g, path = image.path, width = 6, height = 3)

## parallel coordinate plot
g = parallel_coord(diatom.X, diatom.df$Depth, diatom.info[as.character(colnames(diatom.X)),'code']) +
  scale_color_gradientn(colors = c('blue','magenta','red','orange','green'),
                        values = rescale(c(56.53,65.04,67.03,69.74,83.99), to = c(0,1)),
                        name = 'Depth',
                        guide = guide_colorbar(direction = "vertical", reverse = TRUE)) +
  labs(x = '', y = 'Proportion') +
  theme(legend.position = 'right', legend.key.size = unit(10, 'points'), legend.title = element_text(size = 10)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.key.size = unit(12, 'points'),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.spacing.y = unit(1,'points')) +
  theme(axis.text.x = element_text(color = diatom.info[as.character(colnames(diatom.X)),'color']))

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
  scale_color_gradientn(colors = c('blue','magenta','red','orange','green'),
                        values = rescale(c(56.53,65.04,67.03,69.74,83.99), to = c(0,1)),
                        name = 'Depth',
                        guide = guide_colorbar(direction = "vertical", reverse = TRUE)) +
  theme(legend.key.size = unit(8, 'points'),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 8),
        legend.spacing.y = unit(1,'points'),
        plot.margin = unit(rep(15,4),'points'),
        axis.title.x = element_text(size = 8))

diatom.legend.point = g
ggsave('1. Depth_distribution.jpeg', g, path = image.path, width = 8, height = 3)

## color, shape, size vectors
diatom.col = unique(ggplot_build(diatom.legend.point)$data[[1]]$colour)

# diatom.shape = rep(1,71)
# diatom.size = rep(1,71)
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

g.psas = line_gpairs(cbind(diatom.res$psas$scores[,1:4], data.frame(Depth = diatom.df$Depth)),
                     col.point = diatom.col, col.line = diatom.col,
                     shape = diatom.shape, size = diatom.size)
g.psao = line_gpairs(cbind(diatom.res$psao$scores[,1:4], data.frame(Depth = diatom.df$Depth)),
                     col.point = diatom.col, col.line = diatom.col,
                     shape = diatom.shape, size = diatom.size)
g.pca = line_gpairs(cbind(diatom.res$pca$scores[,1:4], data.frame(Depth = diatom.df$Depth)),
                    col.point = diatom.col, col.line = diatom.col,
                    shape = diatom.shape, size = diatom.size)
g.power_pca = line_gpairs(cbind(diatom.res$power_pca$scores[,1:4], data.frame(Depth = diatom.df$Depth)),
                          col.point = diatom.col, col.line = diatom.col,
                          shape = diatom.shape, size = diatom.size)
g.apca = line_gpairs(cbind(diatom.res$apca$scores[,1:4], data.frame(Depth = diatom.df$Depth)),
                     col.point = diatom.col, col.line = diatom.col,
                     shape = diatom.shape, size = diatom.size)


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

g = plot_grid(getPlot(g.psao, 5, 1),
              getPlot(g.psao, 5, 2),
              getPlot(g.psao, 5, 3),
              getPlot(g.psao, 5, 4),
              get_legend(diatom.legend.point + theme(legend.key.size = unit(14, 'points'),
                                                     legend.text = element_text(size = 12),
                                                     legend.title = element_text(size = 14),
                                                     legend.spacing.y = unit(1,'points'))),
              nrow = 1, rel_widths = c(2,2,2,2,1))
ggsave('2. Score_psao_diag.jpeg', g, path = image.path, width = 14.5, height = 3)



add_gap <- function(g){
  gg = ggmatrix_gtable(gpairs_lower(g))
  gg$heights[[15]] = unit(0.4, 'null')
  gg
}

g.psas = add_gap(g.psas)
g.psao = add_gap(g.psao)
g.pca = add_gap(g.pca)
g.power_pca = add_gap(g.power_pca)
g.apca = add_gap(g.apca)


g = plot_grid(matrix.legend,
              g.psas, g.psao, g.pca, g.power_pca, g.apca, nrow = 2)
ggsave('2. Score_plot_matrix.jpeg', g, path = image.path, width = 18, height = 12)


g.psas = plot_grid(g.psas, matrix.legend, nrow = 1, rel_widths = c(5,1))
g.psao = plot_grid(g.psao, matrix.legend, nrow = 1, rel_widths = c(5,1))
g.pca = plot_grid(g.pca, matrix.legend, nrow = 1, rel_widths = c(5,1))
g.power_pca = plot_grid(g.power_pca, matrix.legend, nrow = 1, rel_widths = c(5,1))
g.apca = plot_grid(g.apca, matrix.legend, nrow = 1, rel_widths = c(5,1))

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
  geom_point(col = diatom.col, shape = diatom.shape, size = diatom.size*4/3, stroke = 1) +
  theme_bw() +
  theme(legend.position = 'none')

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
  geom_point(col = diatom.col, shape = diatom.shape, size = diatom.size*4/3, stroke = 1) +
  theme_bw() +
  theme(legend.position = 'none')

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

rss = rbind(diatom.res$psas$rss,
            diatom.res$psao$rss,
            diatom.res$pca$rss,
            diatom.res$power_pca$rss,
            diatom.res$apca$rss)
rownames(rss) = c('PSA-S','PSA-O','PCA','power transform','log-ratio')
g = plot_variance_explained(rss, 6)
ggsave('variance_explained.jpeg', g, path = image.path, width = 8, height = 4)
