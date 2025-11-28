require(devtools)
devtools::install_github('haneone33/psacomp')
library(psacomp)

require(compositions)
require(ggplot2)
require(ggtern)
require(GGally)
require(cowplot) # plot_grid
require(reshape) # parallel coordinate plot
library(viridis)
library(scales)
library(ggnewscale)
library(dplyr)

invisible(lapply(list.files('utils', pattern = '.R', full.names = T), source))

data.path = 'Diatom/Data_updated/'
image.path = 'Diatom/Figures_updated/'
dir.create(data.path, showWarnings = T)
dir.create(image.path, showWarnings = T)

################################################################################
## diatom preparation

## read data
diatom.df = read.csv(paste0(data.path, 'diatom.csv'))
diatom.info = read.csv(paste0(data.path,'diatom codes.csv'))

diatom.X = as.matrix(diatom.df[,-1])
colnames(diatom.X) = diatom.info$code
diatom.X = diatom.X[,names(sort(colMeans(diatom.X), decreasing = T))] # sort by mean proportion

rownames(diatom.info) = diatom.info$code
diatom.info$class[diatom.info$class == ''] = 'others'
diatom.info$class = factor(diatom.info$class, levels = c('warm water','open ocean','sea ice','others'))
diatom.info$color = c('hotpink','green3','dodgerblue','black')[diatom.info$class]

## PCA application
diatom.res = list(X = diatom.X, info = diatom.info)

# system.time({diatom.res$psas = psa('s', diatom.X)})
# saveRDS(diatom.res$psas, paste0(data.path, 'diatom_psas.rds'))
# system.time({diatom.res$psao = psa('o', diatom.X)})
# saveRDS(diatom.res$psao, paste0(data.path, 'diatom_psao.rds'))
diatom.res$psas = readRDS(paste0(data.path, 'diatom_psas.rds'))
diatom.res$psao = readRDS(paste0(data.path, 'diatom_psao.rds'))

diatom.res$pca = comp_pca(diatom.X)
diatom.res$power_pca = comp_power_pca(diatom.X, 0.5)
diatom.res$apca = comp_apca(diatom.X)

## flip direction of loadings and scores in accordance with PSA
diatom.res$pca$loadings[,c(1,3,4)] = -diatom.res$pca$loadings[,c(1,3,4)]
diatom.res$pca$scores[,c(1,3,4)] = -diatom.res$pca$scores[,c(1,3,4)]
diatom.res$power_pca$loadings[,c(1,3,4)] = -diatom.res$power_pca$loadings[,c(1,3,4)]
diatom.res$power_pca$scores[,c(1,3,4)] = -diatom.res$power_pca$scores[,c(1,3,4)]
diatom.res$apca$loadings[,c(2,4)] = -diatom.res$apca$loadings[,c(2,4)]
diatom.res$apca$scores[,c(2,4)] = -diatom.res$apca$scores[,c(2,4)]

################################################################################
## descriptive figures

## heatmap
diatom.X.log = log10(replace_zero(diatom.X))

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
                       values = rescale(c(max(diatom.X.log), max(apply(diatom.X.log, 2, min)), min(diatom.X.log))),
                       name = 'log10(Proportion)') +
  labs(y = 'Depth (Index)') +
  scale_y_continuous(expand = c(0,0,0,0)) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_x_discrete(labels = diatom.info[as.character(colnames(diatom.X)),'code']) +
  theme(axis.text.x = element_text(color = diatom.info[as.character(colnames(diatom.X)),'color'],
                                   angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.title = element_text(hjust=0.5))

g1 = data.frame(variable = 'Depth',
                Depth = diatom.df$Depth) %>%
  ggplot(aes(x = variable, y = rank(Depth), fill = Depth)) +
  theme_minimal() +
  geom_tile() +
  scale_fill_gradientn(colors = c('blue','magenta','red','orange','green'),
                       values = rescale(c(56.53,65.04,67.03,69.74,83.99), to = c(0,1)),
                       guide = guide_colorbar(direction = "vertical", reverse = TRUE)) +
  scale_y_continuous(trans = 'reverse',
                     labels = diatom.df$Depth[c(1,18,36,54,71)], breaks = c(1,18,36,54,71),
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

## score 1 vs score 2
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

## PSA-O scores
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

## first four scores and depth
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
plot_loadings_diatom <- function(V, k = 4, max.k = 12){
  ls = list()
  for(i in 1:k){
    gg = plot_vertex(V[,i], max.k) +
      ggtitle(colnames(V)[i]) +
      theme(axis.text.y = element_text(size=12))
    ls[[i]] = gg + theme(axis.text.y = element_text(color = diatom.info[as.character(gg$data$variable),'color']))
  }
  g = cowplot::plot_grid(plotlist = ls, nrow = 1, align = 'v', axis = 'l')
  return(g)
}

g = plot_loadings_diatom(diatom.res$psas$loadings)
ggsave('3. loading_bar_psas.jpeg', g, path = image.path, width = 16, height = 4)
g = plot_loadings_diatom(diatom.res$psao$loadings)
ggsave('3. loading_bar_psao.jpeg', g, path = image.path, width = 16, height = 4)

## loading bar plot for 5 methods

g = plot_grid(ggdraw() + draw_label('PSA-S'), plot_loadings_diatom(diatom.res$psas$loadings),
              ggdraw() + draw_label('PSA-O'), plot_loadings_diatom(diatom.res$psao$loadings),
              ggdraw() + draw_label('PCA'), plot_loadings_diatom(diatom.res$pca$loadings),
              ggdraw() + draw_label('Power Transform PCA'), plot_loadings_diatom(diatom.res$power_pca$loadings),
              ggdraw() + draw_label('Log-ratio PCA'), plot_loadings_diatom(diatom.res$apca$loadings),
          ncol=2, align='v', axis='lrtb', rel_widths = c(1,6))
ggsave('3. loading_bar_all.jpeg', g, path = image.path, width = 16, height = 16)


## ternary plot
psas.tern = as.data.frame(diatom.res$psas$Xhat_reduced$`r=2`) %>%
  mutate(Depth = diatom.df$Depth) %>%
  ggtern(aes(x=V1, y=V3, z=V2)) +
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

psao.tern = as.data.frame(diatom.res$psao$Xhat_reduced$`r=2`) %>%
  mutate(Depth = diatom.df$Depth) %>%
  ggtern(aes(x=V1, y=V3, z=V2)) +
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
g = plot_loadings_diatom(diatom.res$psas$Vhat$`r=2`, 3)
ggsave('3-3. ternary_vertex_psas.jpeg', g, path = image.path, width = 12, height = 4)
g = plot_loadings_diatom(diatom.res$psao$Vhat$`r=2`, 3)
ggsave('3-4. ternary_vertex_psao.jpeg', g, path = image.path, width = 12, height = 4)

## RSS
RSS = rbind(c(diatom.res$psas$RSS,0),
            c(diatom.res$psao$RSS,0),
            diatom.res$pca$RSS,
            diatom.res$power_pca$RSS,
            diatom.res$apca$RSS)
rownames(RSS) = c('PSA-S','PSA-O','PCA','Power transform','Log-ratio')
g = plot_variance_explained(RSS, 6)
ggsave('variance_explained.jpeg', g, path = image.path, width = 8, height = 4)

