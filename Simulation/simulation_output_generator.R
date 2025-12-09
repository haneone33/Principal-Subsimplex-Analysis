source('PSA_setup.R')

data.path = 'Simulation/Data/'
image.path = 'Simulation/Figures/'
dir.create(data.path, showWarnings = F)
dir.create(image.path, showWarnings = F)

################################################################################
## Example 1
################################################################################
## data generation
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
colnames(X1) = paste0('V',1:3)
ex.label = as.factor(c(rep(1,ns[1]), rep(2,ns[2]), rep(3,ns[3]), rep(4,ns[4])))

ex1.df = as.data.frame(X1)
colnames(ex1.df) = c('V1','V2','V3')
ex1.df$label = ex.label
write.csv(ex1.df, paste0(data.path, 'ex1.csv'), row.names = F)

## PSA application

ex1.psas = psa(X1, 's')
ex1.psao = psa(X1, 'o')
ex1.pca = comp_pca(X1)
ex1.power_pca = comp_power_pca(X1, 0.5)
ex1.apca = comp_apca(X1)

ex1.pca = flip_loading(ex1.pca, c(1))
ex1.power_pca = flip_loading(ex1.power_pca, c(1))
ex1.apca = flip_loading(ex1.apca, c(1,2))

################################################################################
## Example 1 figures
################################################################################
## ternary plot with rank 1 approximation
g.data = ternary_pc(X1, ex.label, type = 'data')
g.psas = ternary_pc(X1, ex.label, type = 'psas', X.res = ex1.psas)
g.psao = ternary_pc(X1, ex.label, type = 'psao', X.res = ex1.psao)
g.pca = ternary_pc(X1, ex.label, type = 'pca', X.res = ex1.pca)
g.power_pca = ternary_pc(X1, ex.label, type = 'power', X.res = ex1.power_pca)
g.apca = ternary_pc(X1, ex.label, type = 'logratio', X.res = ex1.apca)
g = plot_grid(g.data, g.psas, g.psao, g.pca, g.power_pca, g.apca, nrow = 2)

ggsave('ex1_ternary_pc.jpeg', g, path = image.path, width = 12, height = 8)

## score scatter plot matrices

g.data <- ggtern(ex1.df, aes(x = V1, y = V3, z = V2, col = label)) +
  geom_point() +
  theme_bw() +
  ggtitle('Data') +
  theme(plot.title = element_text(hjust = 0.5, size = 15)) +
  theme(legend.position = 'none') +
  labs(T = 'V3', L = 'V1', R = 'V2', x = NULL, y = NULL)
g.data <- ggplot_gtable(ggtern::ggplot_build(g.data))

g.psas = point_gpairs(ex1.psas$scores, ex.label, 2) +
  ggtitle('PSA-S') +
  theme(plot.title = element_text(hjust = 0.5, size = 15))
g.psao = point_gpairs(ex1.psao$scores, ex.label, 2) +
  ggtitle('PSA-O') +
  theme(plot.title = element_text(hjust = 0.5, size = 15))
g.pca = point_gpairs(ex1.pca$scores, ex.label, 2) +
  ggtitle('PCA') +
  theme(plot.title = element_text(hjust = 0.5, size = 15))
g.power = point_gpairs(ex1.power_pca$scores, ex.label, 2) +
  ggtitle('Power Transform PCA') +
  theme(plot.title = element_text(hjust = 0.5, size = 15))
g.apca = point_gpairs(ex1.apca$scores, ex.label, 2) +
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
g.main = plot_grid(plot_vertex(ex1.psas$loadings[,1]) + theme(plot.margin = unit(c(5, 5, -10, 5), 'pt')),
                   plot_vertex(ex1.psas$loadings[,2]) + theme(plot.margin = unit(c(5, 5, -10, 5), 'pt')),
                   plot_vertex(ex1.psao$loadings[,1]) + theme(plot.margin = unit(c(5, 5, -10, 5), 'pt')),
                   plot_vertex(ex1.psao$loadings[,2]) + theme(plot.margin = unit(c(5, 5, -10, 5), 'pt')),
                   plot_vertex(ex1.pca$loadings[,1]) + theme(plot.margin = unit(c(5, 5, -10, 5), 'pt')),
                   plot_vertex(ex1.pca$loadings[,2]) + theme(plot.margin = unit(c(5, 5, -10, 5), 'pt')),
                   plot_vertex(ex1.power_pca$loadings[,1]) + theme(plot.margin = unit(c(5, 5, -10, 5), 'pt')),
                   plot_vertex(ex1.power_pca$loadings[,2]) + theme(plot.margin = unit(c(5, 5, -10, 5), 'pt')),
                   plot_vertex(ex1.apca$loadings[,1]) + theme(plot.margin = unit(c(5, 5, -10, 5), 'pt')),
                   plot_vertex(ex1.apca$loadings[,2]) + theme(plot.margin = unit(c(5, 5, -10, 5), 'pt')),
                   nrow = 2, byrow = F)

g = plot_grid(ggdraw(), title.col,
              title.row, g.main, nrow = 2,
              rel_heights = c(0.3, 2), rel_widths = c(0.4,5))

ggsave('ex1_loading_bar.jpeg', g, path = image.path, width = 12, height = 2.3)


ex1.rss = rbind(c(ex1.psas$RSS, 'PC3'=0),
                c(ex1.psao$RSS, 'PC3'=0),
                ex1.pca$RSS,
                ex1.power_pca$RSS,
                ex1.apca$RSS)
rownames(ex1.rss) = c('PSA-S','PSA-O','PCA','Power transform','Log-ratio')
g = plot_variance_explained(ex1.rss)
ggsave('ex1_variance_explained.jpeg', g, path = image.path, width = 8, height = 3)

################################################################################
## Example 2
################################################################################
## data generation
set.seed(1)
X2 = cbind(X1, matrix(rnorm(30*3, mean = 0, sd = sqrt(sigma2)), ncol = 3))
colnames(X2) = paste0('V',1:6)
X2 = to_simplex(X2)

ex2.df = as.data.frame(X2)
colnames(ex2.df) = paste0('V',1:6)
ex2.df$label = ex.label
write.csv(ex2.df, paste0(data.path, 'ex2.csv'), row.names = F)

## PSA application
ex2.psas = psa(X2, 's')
ex2.psao = psa(X2 ,'o')
ex2.pca = comp_pca(X2)
ex2.power_pca = comp_power_pca(X2, 0.5)
ex2.apca = comp_apca(X2)

ex2.pca = flip_loading(ex2.pca, c(1,3))
ex2.power_pca = flip_loading(ex2.power_pca, c(1,3))
ex2.apca = flip_loading(ex2.apca, c(1,2,3))

## score scatter plot matrices
g.data <- ggplot() + theme_void()

g.psas = point_gpairs(ex2.psas$scores, ex.label, k = 3) +
  ggtitle('PSA-S') +
  theme(plot.title = element_text(hjust = 0.5, size = 15))
g.psao = point_gpairs(ex2.psao$scores, ex.label, k = 3) +
  ggtitle('PSA-O') +
  theme(plot.title = element_text(hjust = 0.5, size = 15))
g.pca = point_gpairs(ex2.pca$scores, ex.label, k = 3) +
  ggtitle('PCA') +
  theme(plot.title = element_text(hjust = 0.5, size = 15))
g.power = point_gpairs(ex2.power_pca$scores, ex.label, k = 3) +
  ggtitle('Power Transform PCA') +
  theme(plot.title = element_text(hjust = 0.5, size = 15))
g.apca = point_gpairs(ex2.apca$scores, ex.label, k = 3) +
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
g.main = plot_grid(plot_vertex(ex2.psas$loadings[,1]) + theme(plot.margin = unit(c(5, 5, -10, 5), 'pt')),
                   plot_vertex(ex2.psas$loadings[,2]) + theme(plot.margin = unit(c(5, 5, -10, 5), 'pt')),
                   plot_vertex(ex2.psas$loadings[,3]) + theme(plot.margin = unit(c(5, 5, -10, 5), 'pt')),
                   plot_vertex(ex2.psas$loadings[,4]) + theme(plot.margin = unit(c(5, 5, -10, 5), 'pt')),
                   plot_vertex(ex2.psao$loadings[,1]) + theme(plot.margin = unit(c(5, 5, -10, 5), 'pt')),
                   plot_vertex(ex2.psao$loadings[,2]) + theme(plot.margin = unit(c(5, 5, -10, 5), 'pt')),
                   plot_vertex(ex2.psao$loadings[,3]) + theme(plot.margin = unit(c(5, 5, -10, 5), 'pt')),
                   plot_vertex(ex2.psao$loadings[,4]) + theme(plot.margin = unit(c(5, 5, -10, 5), 'pt')),
                   plot_vertex(ex2.pca$loadings[,1]) + theme(plot.margin = unit(c(5, 5, -10, 5), 'pt')),
                   plot_vertex(ex2.pca$loadings[,2]) + theme(plot.margin = unit(c(5, 5, -10, 5), 'pt')),
                   plot_vertex(ex2.pca$loadings[,3]) + theme(plot.margin = unit(c(5, 5, -10, 5), 'pt')),
                   plot_vertex(ex2.pca$loadings[,4]) + theme(plot.margin = unit(c(5, 5, -10, 5), 'pt')),
                   plot_vertex(ex2.power_pca$loadings[,1]) + theme(plot.margin = unit(c(5, 5, -10, 5), 'pt')),
                   plot_vertex(ex2.power_pca$loadings[,2]) + theme(plot.margin = unit(c(5, 5, -10, 5), 'pt')),
                   plot_vertex(ex2.power_pca$loadings[,3]) + theme(plot.margin = unit(c(5, 5, -10, 5), 'pt')),
                   plot_vertex(ex2.power_pca$loadings[,4]) + theme(plot.margin = unit(c(5, 5, -10, 5), 'pt')),
                   plot_vertex(ex2.apca$loadings[,1]) + theme(plot.margin = unit(c(5, 5, -10, 5), 'pt')),
                   plot_vertex(ex2.apca$loadings[,2]) + theme(plot.margin = unit(c(5, 5, -10, 5), 'pt')),
                   plot_vertex(ex2.apca$loadings[,3]) + theme(plot.margin = unit(c(5, 5, -10, 5), 'pt')),
                   plot_vertex(ex2.apca$loadings[,4]) + theme(plot.margin = unit(c(5, 5, -10, 5), 'pt')),
                   nrow = 4, byrow = F)

g = plot_grid(ggdraw(), title.col,
              title.row, g.main, nrow = 2,
              rel_heights = c(0.3, 4), rel_widths = c(0.4,5))

ggsave('ex2_loading_bar.jpeg', g, path = image.path, width = 12, height = 4.3)


ex2.rss = rbind(c(ex2.psas$RSS, 'PC6'=0),
                c(ex2.psao$RSS, 'PC6'=0),
                ex2.pca$RSS,
                ex2.power_pca$RSS,
                ex2.apca$RSS)
rownames(ex2.rss) = c('PSA-S','PSA-O','PCA','Power transform','Log-ratio')
g = plot_variance_explained(ex2.rss)
ggsave('ex2_variance_explained.jpeg', g, path = image.path, width = 8, height = 3)

################################################################################
## Example 1 modes of variation
################################################################################
X.data = to_2d(X1)
X.data$label = ex.label

## data
g0 = empty_tern(c('V1','V2','V3')) +
  geom_point(data = X.data, aes(col = label), shape = 16) +
  ggtitle('Data Distribution') +
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        legend.position = 'none')

## PSAS
v1 = ex1.psas$Vhat$`r=1`[,2]
v2 = ex1.psas$Vhat$`r=1`[,1]
mu = ex1.psas$backwards_mean
vdiff = c(1,0,-1)
u1 = mu - mu[1]/vdiff[1]*vdiff
u2 = mu - mu[3]/vdiff[3]*vdiff
X.center = t(as.data.frame(ex1.psas$backwards_mean))
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
  annotate(geom = 'text', x = 1.1, y = 0.2,
           label = TeX(r'($=\hat{V}^{(1)}_2$)')) +
  annotate(geom = 'text', x = 0, y = 0.5,
           label = TeX(r'($\hat{V}^{(0)}_1$)'))

## PSAO
v1 = ex1.psao$Vhat$`r=1`[,2]
v2 = ex1.psao$Vhat$`r=1`[,1]
nvec = c(-v2[3], 0, v2[1])
mu = ex1.psao$backwards_mean
u1 = nvec - nvec[1]/mu[1]*mu
u2 = nvec - nvec[3]/mu[3]*mu
X.center = t(as.data.frame(ex1.psao$backwards_mean))
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
  v = ex1.psao$backwards_mean
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
  annotate(geom = 'text', x = 1.1, y = 0.2,
           label = TeX(r'($=\hat{V}^{(1)}_2$)')) +
  annotate(geom = 'text', x = 0, y = 0.5,
           label = TeX(r'($\hat{V}^{(0)}_1$)'))

g = plot_grid(g0, g1, g2, nrow = 1)
ggsave('ex1_PSA_modes.jpeg', g, path = image.path, width = 13*0.75, height = 4*0.8)

