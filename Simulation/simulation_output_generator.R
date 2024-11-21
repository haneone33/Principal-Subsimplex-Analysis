devtools::install_github('haneone33/Principal-Nested-Simplices', subdir = 'PSA')
library(PSA)

require(compositions)
require(ggplot2)
require(ggtern)
require(GGally)
require(scales)
require(cowplot) # plot_grid
require(reshape) # parallel coordinate plot

invisible(lapply(list.files('utils', pattern = '.R', full.names = T), source))
source('Simulation/simulation_data_generating_functions.R')

data.path = 'Simulation/Data/'
image.path = 'Simulation/Figures/'

################################################################################
## Example 1 data generation
centers = list(c(0.05, 0.05, 0.9),
               c(0.05, 0.9, 0.05),
               c(0.9, 0.05, 0.05),
               c(0.25, 0.7, 0.05))
ns = c(5, 10, 10, 5)
sigma2 = 0.04^2

set.seed(1)
X1 = rbind(make_cluster(centers[[1]], ns[1], sigma2),
           make_cluster(centers[[2]], ns[2], sigma2),
           make_cluster(centers[[3]], ns[3], sigma2),
           make_cluster(centers[[4]], ns[4], sigma2))
ex.label = as.factor(c(rep(1,ns[1]), rep(2,ns[2]), rep(3,ns[3]), rep(4,ns[4])))

df = as.data.frame(cbind(ex.label, X1))
colnames(df) = c('label','V1','V2','V3')
write.csv(df, paste0(data.path, 'ex1.csv'), row.names = F)

## PSA application
ex1.res = compare_analysis(X1)

################################################################################
## Example 1 figures

## ternary plot with rank 1 approximation
g.data = ternary_pc(ex1.res, ex.label, 'data')
g.psas = ternary_pc(ex1.res, ex.label, 'psas')
g.psao = ternary_pc(ex1.res, ex.label, 'psao')
g.pca = ternary_pc(ex1.res, ex.label, 'pca')
g.power_pca = ternary_pc(ex1.res, ex.label, 'power')
g.apca = ternary_pc(ex1.res, ex.label, 'logratio')
g = plot_grid(g.data, g.psas, g.psao, g.pca, g.power_pca, g.apca, nrow = 2)

ggsave('ex1_ternary_pc.jpeg', g, path = image.path, width = 12, height = 8)

## score scatter plot matrices
X.df = as.data.frame(ex1.res$X)
X.df$label = ex.label
g.data <- ggtern(X.df, aes(x = V1, y = V3, z = V2, col = label)) +
  geom_point() +
  theme_classic() +
  ggtitle('Data') +
  theme(plot.title = element_text(hjust = 0.5, size = 15)) +
  theme(legend.position = 'none')
g.data <- ggplot_gtable(ggtern::ggplot_build(g.data))
  
g.psas = point_gpairs(ex1.res$psas$scores, ex.label) +
  ggtitle('PSA-S') +
  theme(plot.title = element_text(hjust = 0.5, size = 15))
g.psao = point_gpairs(ex1.res$psao$scores, ex.label) +
  ggtitle('PSA-O') +
  theme(plot.title = element_text(hjust = 0.5, size = 15))
g.pca = point_gpairs(ex1.res$pca$scores, ex.label) +
  ggtitle('PCA') +
  theme(plot.title = element_text(hjust = 0.5, size = 15))
g.power = point_gpairs(ex1.res$power_pca$scores, ex.label) +
  ggtitle('Power Transform PCA') +
  theme(plot.title = element_text(hjust = 0.5, size = 15))
g.apca = point_gpairs(ex1.res$apca$scores, ex.label) +
  ggtitle('Log-ratio PCA') +
  theme(plot.title = element_text(hjust = 0.5, size = 15))

g = plot_grid(g.data,
              ggmatrix_gtable(g.psas),
              ggmatrix_gtable(g.psao),
              ggmatrix_gtable(g.pca), 
              ggmatrix_gtable(g.power),
              ggmatrix_gtable(g.apca), nrow = 2)

ggsave('ex1_scores.jpeg', g, path = image.path, width = 9, height = 6)

## loading plots
title.row = plot_grid(ggdraw() + draw_label('Comp.1'),
                      ggdraw() + draw_label('Comp.2'), ncol = 1)
title.col = plot_grid(ggdraw(), ggdraw() + draw_label('PSA-S'),
                      ggdraw(), ggdraw() + draw_label('PSA-O'),
                      ggdraw(), ggdraw() + draw_label('PCA'),
                      ggdraw(), ggdraw() + draw_label('Power Transform PCA'),
                      ggdraw(), ggdraw() + draw_label('Log-Ratio PCA'),
                      nrow = 1, rel_widths = rep(c(0.3,1),5))
g.main = plot_grid(loading_bar(ex1.res$psas$loadings[,1]) + theme(plot.margin = unit(c(5, 5, -10, 5), 'pt')),
                   loading_bar(ex1.res$psas$loadings[,2]) + theme(plot.margin = unit(c(5, 5, -10, 5), 'pt')),
                   loading_bar(ex1.res$psao$loadings[,1]) + theme(plot.margin = unit(c(5, 5, -10, 5), 'pt')),
                   loading_bar(ex1.res$psao$loadings[,2]) + theme(plot.margin = unit(c(5, 5, -10, 5), 'pt')),
                   loading_bar(ex1.res$pca$loadings[,1]) + theme(plot.margin = unit(c(5, 5, -10, 5), 'pt')),
                   loading_bar(ex1.res$pca$loadings[,2]) + theme(plot.margin = unit(c(5, 5, -10, 5), 'pt')),
                   loading_bar(ex1.res$power_pca$loadings[,1]) + theme(plot.margin = unit(c(5, 5, -10, 5), 'pt')),
                   loading_bar(ex1.res$power_pca$loadings[,2]) + theme(plot.margin = unit(c(5, 5, -10, 5), 'pt')),
                   loading_bar(ex1.res$apca$loadings[,1]) + theme(plot.margin = unit(c(5, 5, -10, 5), 'pt')),
                   loading_bar(ex1.res$apca$loadings[,2]) + theme(plot.margin = unit(c(5, 5, -10, 5), 'pt')),
                   nrow = 2, byrow = F)

g = plot_grid(ggdraw(), title.col,
              title.row, g.main, nrow = 2,
              rel_heights = c(0.3, 2), rel_widths = c(0.4,5))

ggsave('ex1_loading_bar.jpeg', g, path = image.path, width = 12, height = 2.3)

################################################################################
## Example 2

## data generation
set.seed(1)
X2 = cbind(X1, matrix(rnorm(30*3, mean = 0, sd = sqrt(sigma2)), ncol = 3))
X2 = to_simplex(X2)

df = as.data.frame(cbind(ex.label, X2))
colnames(df) = c('label', 'V1','V2','V3','V4','V5','V6')
write.csv(df, paste0(data.path, 'ex2.csv'), row.names = F)

## PSA application
ex2.res = compare_analysis(X2)

## score scatter plot matrices
g.data <- ggplot() + theme_void()

g.psas = point_gpairs(ex2.res$psas$scores, ex.label, k = 3) +
  ggtitle('PSA-S') +
  theme(plot.title = element_text(hjust = 0.5, size = 15))
g.psao = point_gpairs(ex2.res$psao$scores, ex.label, k = 3) +
  ggtitle('PSA-O') +
  theme(plot.title = element_text(hjust = 0.5, size = 15))
g.pca = point_gpairs(ex2.res$pca$scores, ex.label, k = 3) +
  ggtitle('PCA') +
  theme(plot.title = element_text(hjust = 0.5, size = 15))
g.power = point_gpairs(ex2.res$power_pca$scores, ex.label, k = 3) +
  ggtitle('Power Transform PCA') +
  theme(plot.title = element_text(hjust = 0.5, size = 15))
g.apca = point_gpairs(ex2.res$apca$scores, ex.label, k = 3) +
  ggtitle('Log-ratio PCA') +
  theme(plot.title = element_text(hjust = 0.5, size = 15))

g = plot_grid(g.data,
              ggmatrix_gtable(g.psas),
              ggmatrix_gtable(g.psao),
              ggmatrix_gtable(g.pca), 
              ggmatrix_gtable(g.power),
              ggmatrix_gtable(g.apca), nrow = 2)

ggsave('ex2_scores.jpeg', g, path = image.path, width = 9, height = 6)

## loading plots
title.row = plot_grid(ggdraw() + draw_label('Comp.1'),
                      ggdraw() + draw_label('Comp.2'),
                      ggdraw() + draw_label('Comp.3'),
                      ggdraw() + draw_label('Comp.4'), ncol = 1)
title.col = plot_grid(ggdraw(), ggdraw() + draw_label('PSA-S'),
                      ggdraw(), ggdraw() + draw_label('PSA-O'),
                      ggdraw(), ggdraw() + draw_label('PCA'),
                      ggdraw(), ggdraw() + draw_label('Power Transform PCA'),
                      ggdraw(), ggdraw() + draw_label('Log-Ratio PCA'),
                      nrow = 1, rel_widths = rep(c(0.2,1),5))
g.main = plot_grid(loading_bar(ex2.res$psas$loadings[,1]) + theme(plot.margin = unit(c(5, 5, -10, 5), 'pt')),
                   loading_bar(ex2.res$psas$loadings[,2]) + theme(plot.margin = unit(c(5, 5, -10, 5), 'pt')),
                   loading_bar(ex2.res$psas$loadings[,3]) + theme(plot.margin = unit(c(5, 5, -10, 5), 'pt')),
                   loading_bar(ex2.res$psas$loadings[,4]) + theme(plot.margin = unit(c(5, 5, -10, 5), 'pt')),
                   loading_bar(ex2.res$psao$loadings[,1]) + theme(plot.margin = unit(c(5, 5, -10, 5), 'pt')),
                   loading_bar(ex2.res$psao$loadings[,2]) + theme(plot.margin = unit(c(5, 5, -10, 5), 'pt')),
                   loading_bar(ex2.res$psao$loadings[,3]) + theme(plot.margin = unit(c(5, 5, -10, 5), 'pt')),
                   loading_bar(ex2.res$psao$loadings[,4]) + theme(plot.margin = unit(c(5, 5, -10, 5), 'pt')),
                   loading_bar(ex2.res$pca$loadings[,1]) + theme(plot.margin = unit(c(5, 5, -10, 5), 'pt')),
                   loading_bar(ex2.res$pca$loadings[,2]) + theme(plot.margin = unit(c(5, 5, -10, 5), 'pt')),
                   loading_bar(ex2.res$pca$loadings[,3]) + theme(plot.margin = unit(c(5, 5, -10, 5), 'pt')),
                   loading_bar(ex2.res$pca$loadings[,4]) + theme(plot.margin = unit(c(5, 5, -10, 5), 'pt')),
                   loading_bar(ex2.res$power_pca$loadings[,1]) + theme(plot.margin = unit(c(5, 5, -10, 5), 'pt')),
                   loading_bar(ex2.res$power_pca$loadings[,2]) + theme(plot.margin = unit(c(5, 5, -10, 5), 'pt')),
                   loading_bar(ex2.res$power_pca$loadings[,3]) + theme(plot.margin = unit(c(5, 5, -10, 5), 'pt')),
                   loading_bar(ex2.res$power_pca$loadings[,4]) + theme(plot.margin = unit(c(5, 5, -10, 5), 'pt')),
                   loading_bar(ex2.res$apca$loadings[,1]) + theme(plot.margin = unit(c(5, 5, -10, 5), 'pt')),
                   loading_bar(ex2.res$apca$loadings[,2]) + theme(plot.margin = unit(c(5, 5, -10, 5), 'pt')),
                   loading_bar(ex2.res$apca$loadings[,3]) + theme(plot.margin = unit(c(5, 5, -10, 5), 'pt')),
                   loading_bar(ex2.res$apca$loadings[,4]) + theme(plot.margin = unit(c(5, 5, -10, 5), 'pt')),
                   nrow = 4, byrow = F)

g = plot_grid(ggdraw(), title.col,
              title.row, g.main, nrow = 2,
              rel_heights = c(0.3, 4), rel_widths = c(0.4,5))

ggsave('ex2_loading_bar.jpeg', g, path = image.path, width = 12, height = 4.3)

################################################################################
## Example 1 modes of variation

X.data = to_2d(X1)
X.data$label = ex.label

## data
g0 = empty_tern(c('V1','V2','V3')) +
  geom_point(data = X.data, aes(col = label), shape = 16) +
  ggtitle('Data Distribution') +
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        legend.position = 'none')

## PSAS
v1 = ex1.res$psas$vertices$`r=2`[2,]
v2 = ex1.res$psas$vertices$`r=2`[1,]
mu = ex1.res$psas$vertices$`r=1`
vdiff = c(1,0,-1)
u1 = mu - mu[1]/vdiff[1]*vdiff
u2 = mu - mu[3]/vdiff[3]*vdiff
X.center = as.data.frame(ex1.res$psas$vertices$`r=1`)
mode1 = as.data.frame(rbind(v1, v2))
mode2 = as.data.frame(rbind(u1, u2))
g1 = empty_tern(c('V1','V2','V3')) +
  geom_point(data = X.data, aes(col = label), shape = 16) +
  geom_path(data = to_2d(mode1), col = 'red') +
  geom_path(data = to_2d(mode2), col = 'blue') +
  geom_point(data = to_2d(X.center)) +
  ggtitle('PSA-S Modes') +
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        legend.position = 'none')

for(t in c(1:9)/10){
  x = c(0, t, 1-t)
  y = c(1-t, t, 0)
  df = to_2d(rbind(x,y))
  g1 = g1 + geom_path(data = df, col = 'darkgray')
}

g1 = g1 + annotate(geom = 'text', x = -0.7, y = 0.7,
                   label = TeX(r'($\hat{V}^{(1)}_1$)')) +
  annotate(geom = 'text', x = 1.1, y = 0.1,
           label = TeX(r'($=\hat{V}^{(1)}_2$)')) +
  annotate(geom = 'text', x = 0, y = 0.5,
           label = TeX(r'($\hat{V}^{(0)}_1$)'))

## PSAO
v1 = ex1.res$psao$vertices$`r=2`[2,]
v2 = ex1.res$psao$vertices$`r=2`[1,]
nvec = c(-v2[3], 0, v2[1])
mu = ex1.res$psao$vertices$`r=1`
u1 = nvec - nvec[1]/mu[1]*mu
u2 = nvec - nvec[3]/mu[3]*mu
X.center = as.data.frame(ex1.res$psao$vertices$`r=1`)
mode1 = as.data.frame(rbind(v1, v2))
mode2 = as.data.frame(rbind(u1, u2))
g2 = empty_tern(c('V1','V2','V3')) +
  geom_point(data = X.data, aes(col = label), shape = 16) +
  geom_path(data = to_2d(mode1), col = 'red') +
  geom_path(data = to_2d(mode2), col = 'blue') +
  geom_point(data = to_2d(X.center)) +
  ggtitle('PSA-O Modes') +
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        legend.position = 'none')

for(t in c(1:9)/10){
  x = c(0, t, 1-t)
  x = x/sqrt(sum(x^2))
  v = ex1.res$psao$vertices$`r=2`[1,]
  names(v) = NULL
  v = c(-v[3],0,v[1])
  z = (x-sum(x*v)*v)/sin(acos(sum(x*v)))
  x = x/sum(x)
  z = z/sum(z)
  y = x - (z-x)/(z-x)[3]*x[3]
  df = to_2d(rbind(x,y))
  g2 = g2 + geom_path(data = df, col = 'darkgray')
}

g2 = g2 + annotate(geom = 'text', x = -0.7, y = 0.7,
                   label = TeX(r'($\hat{V}^{(1)}_1$)')) +
  annotate(geom = 'text', x = 1.1, y = 0.1,
           label = TeX(r'($=\hat{V}^{(1)}_2$)')) +
  annotate(geom = 'text', x = 0, y = 0.5,
           label = TeX(r'($\hat{V}^{(0)}_1$)'))

g = plot_grid(g0, g1, g2, nrow = 1)
ggsave('ex1_PSA_modes.jpeg', g, path = image.path, width = 13*0.75, height = 4*0.8)
