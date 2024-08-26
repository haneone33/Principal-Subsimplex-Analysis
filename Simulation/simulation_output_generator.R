source('PSA_init.R')
source('Simulation/simulation_data_generator.R')

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
Ex.label = as.factor(c(rep(1,ns[1]), rep(2,ns[2]), rep(3,ns[3]), rep(4,ns[4])))

df = as.data.frame(cbind(Ex.label, X1))
colnames(df) = c('label','V1','V2','V3')
write.csv(df, paste0(data.path, 'Ex1.csv'), row.names = F)

## PSA application
Ex1.res = compare_analysis(X1)

################################################################################
## Example 1 figures

## ternary plot with rank 1 approximation
g = compare_ternary(Ex1.res, Ex.label, main = 'Example 1')
ggsave('ternary_pc.png', g, path = image.path, width = 12, height = 8)

## score scatter plot matrices
g = compare_score_plot_matrix(Ex1.res, Ex.label, main = 'Example 1',
                              ternary = F, legend.plot = NULL, plot.diag = T)
ggsave('scores.png', g, path = image.path, width = 9, height = 7)

## PSA modes of variation

## loading plots
g = plot_grid(ggdraw() + draw_label('Example 1 PSA-S'),
              plot_grid(loading_bar(Ex1.res$psas$loadings[,1]) + ggtitle('Comp.1'),
                        loading_bar(Ex1.res$psas$loadings[,2]) + ggtitle('Comp.2'), nrow = 1),
              ncol = 1, rel_heights = c(1,8))
ggsave('loading_bar_psas.png', g, path = image.path, width = 9, height = 3)

g = plot_grid(ggdraw() + draw_label('Example 1 PSA-O'),
              plot_grid(loading_bar(Ex1.res$psao$loadings[,1]) + ggtitle('Comp.1'),
                        loading_bar(Ex1.res$psao$loadings[,2]) + ggtitle('Comp.2'), nrow = 1),
              ncol = 1, rel_heights = c(1,8))
ggsave('loading_bar_psao.png', g, path = image.path, width = 9, height = 3)

## dendrogram
png(paste0(image.path,'PSA-S_dendro.png'), width = 1200, height = 1200, res = 400)
plotdendrogram2(Ex1.res$psas)
title('PSA-S')
dev.off()
png(paste0(image.path,'PSA-O_dendro.png'), width = 1200, height = 1200, res = 400)
plotdendrogram2(Ex1.res$psao)
title('PSA-O')
dev.off()

################################################################################
## Example 2

## data generation
set.seed(1)
X2 = cbind(X1, matrix(rnorm(30*3, mean = 0, sd = sqrt(sigma2)), ncol = 3))
X2 = to_simplex(X2)

df = as.data.frame(cbind(Ex.label, X2))
colnames(df) = c('label', 'V1','V2','V3','V4','V5','V6')
write.csv(df, paste0(data.path, 'Ex2.csv'), row.names = F)

## PSA application
Ex2.res = compare_analysis(X2)

## score scatter plot matrices
g = compare_score_plot_matrix(Ex2.res, Ex.label, main = 'Example 2', k = 3,
                              ternary = F, legend.plot = NULL, plot.diag = T)
ggsave('scores_5dim.png', g, path = image.path, width = 9, height = 7)

################################################################################
## Example 1 modes of variation

X.data = to_2d(X1)
X.data$label = Ex.label

## data
g0 = empty_tern(c('V1','V2','V3')) +
  geom_point(data = X.data, aes(col = label), shape = 16) +
  ggtitle('Data Distribution') +
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        legend.position = 'none')

## PSAS
v1 = Ex1.res$psas$vertices$`r=2`[2,]
v2 = Ex1.res$psas$vertices$`r=2`[1,]
mu = Ex1.res$psas$vertices$`r=1`
vdiff = c(1,0,-1)
u1 = mu - mu[1]/vdiff[1]*vdiff
u2 = mu - mu[3]/vdiff[3]*vdiff
X.center = as.data.frame(Ex1.res$psas$vertices$`r=1`)
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
v1 = Ex1.res$psao$vertices$`r=2`[2,]
v2 = Ex1.res$psao$vertices$`r=2`[1,]
nvec = c(-v2[3], 0, v2[1])
mu = Ex1.res$psao$vertices$`r=1`
u1 = nvec - nvec[1]/mu[1]*mu
u2 = nvec - nvec[3]/mu[3]*mu
X.center = as.data.frame(Ex1.res$psao$vertices$`r=1`)
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
  v = Ex1.res$psao$vertices$`r=2`[1,]
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
ggsave('PSA modes.png', g, path = image.path, width = 13, height = 4)
