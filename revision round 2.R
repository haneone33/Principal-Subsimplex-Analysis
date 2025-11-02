diatom.res = readRDS('../Principal-Subsimplex-Analysis/Diatom/Data/diatom_res.rds')

rss.mat = matrix(NA, 6,6)
rownames(rss.mat) = c('PSAS','PSAO','PCA','PCA_projected','Power_PCA','Logratio')
colnames(rss.mat) = paste0('r=',1:6)

for(r in 1:6){
  r.str = paste0('r=',r+1)
  pca_projected = t(apply(diatom.res$pca$pts.approx[[r.str]], 1, function(x){
    x[x<0] = 0
    x = x/sum(x)
    return(x)
  }))

  power_projected = t(apply(diatom.res$power_pca$pts.approx[[r.str]], 1, function(x){
    x[x<0] = 0
    x = x/sqrt(sum(x**2))
    x = x**2
    return(x)
  }))

  rss.mat[,r] = c(sum(rowSums((diatom.res$psas$pts.approx[[r.str]]-diatom.res$X)**2)),
                  sum(rowSums((diatom.res$psao$pts.approx[[r.str]]-diatom.res$X)**2)),
                  sum(rowSums((diatom.res$pca$pts.approx[[r.str]]-diatom.res$X)**2)),
                  sum(rowSums((pca_projected-diatom.res$X)**2)),
                  sum(rowSums((power_projected-diatom.res$X)**2)),
                  sum(rowSums((diatom.res$apca$pts.approx[[r.str]]-diatom.res$X)**2)))
}


round(rss.mat,4)

library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)


df = rss.mat %>%
  melt() %>%
  setNames(c('Method','Rank','value')) %>%
  mutate(Method = factor(Method, levels = c('PSAS','PSAO','PCA','Power_PCA','Logratio','PCA_projected')))

df %>%
  filter(Method != 'PCA_projected') %>%
  ggplot(aes(x = Rank, y = value, group = Method, shape = Method, linetype = Method)) +
  theme_bw() +
  geom_line() +
  geom_point() +
  geom_line(data = df %>% filter(Method == 'PCA_projected')) +
  geom_point(data = df %>% filter(Method == 'PCA_projected')) +
  scale_linetype_manual(values = c(rep('solid',5), 'dashed')) +
  scale_shape_manual(values = c(16,17,15,3,7,15)) +
  ylab('Residual sums of squares')
ggsave('../Principal-Subsimplex-Analysis/revision round 2/diatom_rss.jpeg', width = 8, height = 4)




range(diatom.res$psas$scores[,1:2])
range(diatom.res$psao$scores[,1:2])
range(diatom.res$pca$scores[,1:2])
range(diatom.res$power_pca$scores[,1:2])
range(diatom.res$apca$scores[,1:2])
g = plot_grid(get_legend(diatom.legend.point + theme(legend.key.size = unit(14, 'points'),
                                                     legend.text = element_text(size = 12),
                                                     legend.title = element_text(size = 14),
                                                     legend.spacing.y = unit(1,'points'))),
              getPlot(g.psas, 2, 1) + ggtitle('PSA-S') +
                scale_x_continuous(limits=c(-0.4,0.63)) + scale_y_continuous(limits=c(-0.4,0.63)),
              getPlot(g.psao, 2, 1) + ggtitle('PSA-O') +
                scale_x_continuous(limits=c(-0.51,0.81)) + scale_y_continuous(limits=c(-0.51,0.81)),
              getPlot(g.pca, 2, 1) + ggtitle('PCA') +
                scale_x_continuous(limits=c(-0.21,0.23)) + scale_y_continuous(limits=c(-0.21,0.23)),
              getPlot(g.power_pca, 2, 1) + ggtitle('Power Transform PCA') +
                scale_x_continuous(limits=c(-0.5,0.55)) + scale_y_continuous(limits=c(-0.5,0.55)),
              getPlot(g.apca, 2, 1) + ggtitle('Log-ratio PCA') +
                scale_x_continuous(limits=c(-7.3,7.3)) + scale_y_continuous(limits=c(-7.3,7.3)),
              nrow = 2, align='hv', axis='tblr')
ggsave('2. Score_1vs2_2.jpeg', g, path = image.path, width = 10.5, height = 7)

draw_biplot <- function(res){
  x = res$scores[,1:2]
  y = res$loadings[,1:2]
  y = y/sqrt(mean(rowSums(y**2))/mean(rowSums(x**2)))/3

  dfx = data.frame(x = x[,1], y = x[,2], Depth = diatom.df$Depth)
  dfy = data.frame(x = 0, y = 0, xend = y[,1], yend = y[,2], size = sqrt(rowSums(y[,1:2]**2)), label=diatom.info[rownames(y),]$code)
  dfy = dfy %>% arrange(desc(size)) %>% head(10)

  ggplot() +
    theme_bw() +
    # geom_hline(yintercept=0, color='gray', linewidth=1) +
    # geom_vline(xintercept=0, color='gray', linewidth=1) +
    geom_point(data=dfx, aes(x=x,y=y,col=Depth)) +
    geom_path(data=dfx, aes(x=x,y=y,col=Depth)) +
    scale_color_gradientn(colors = c('blue','magenta','red','orange','green'),
                          values = rescale(c(56.53,65.04,67.03,69.74,83.99), to = c(0,1)),
                          name = 'Depth',
                          guide = guide_colorbar(direction = "vertical", reverse = TRUE)) +
    geom_segment(data=dfy, aes(x=x,y=y,xend=xend,yend=yend),
                 arrow = arrow(length=unit(5,'points')), linewidth=0.2, color = 'black') +
    geom_text(data=dfy, aes(x=xend, y=yend, label=label), color='black', size=3) +
    labs(x='Comp.1', y='Comp.2') +
    theme(legend.position = 'none',
          plot.title=element_text(hjust=0.5))
}

g = plot_grid(get_legend(diatom.legend.point + theme(legend.key.size = unit(14, 'points'),
                                                     legend.text = element_text(size = 12),
                                                     legend.title = element_text(size = 14),
                                                     legend.spacing.y = unit(1,'points'))),
              draw_biplot(diatom.res$psas) + ggtitle('PSA-S') +
                scale_x_continuous(limits=c(-0.4,0.63)) + scale_y_continuous(limits=c(-0.4,0.63)),
              draw_biplot(diatom.res$psao) + ggtitle('PSA-O') +
                scale_x_continuous(limits=c(-0.51,0.81)) + scale_y_continuous(limits=c(-0.51,0.81)),
              draw_biplot(diatom.res$pca) + ggtitle('PCA') +
                scale_x_continuous(limits=c(-0.21,0.23)) + scale_y_continuous(limits=c(-0.21,0.23)),
              draw_biplot(diatom.res$power_pca) + ggtitle('Power Transform PCA') +
                scale_x_continuous(limits=c(-0.5,0.55)) + scale_y_continuous(limits=c(-0.5,0.55)),
              draw_biplot(diatom.res$apca) + ggtitle('Log-ratio PCA') +
                scale_x_continuous(limits=c(-7.3,7.3)) + scale_y_continuous(limits=c(-7.3,7.3)),
              nrow = 2, align='hv', axis='tblr')
ggsave('2. Score_1vs2_3.jpeg', g, path = image.path, width = 10.5, height = 7)

pca_projected = t(apply(diatom.res$pca$pts.approx[['r=2']], 1, function(x){
  x[x<0] = 0
  x = x/sum(x)
  return(x)
}))

g = parallel_coord(diatom.res$pca$pts.approx$`r=2`, diatom.df$Depth, diatom.info[as.character(colnames(diatom.X)),'code']) +
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
  theme(axis.text.x = element_text(color = diatom.info[as.character(colnames(diatom.X)),'color'])) +
  geom_hline(yintercept = 0, col='black', linewidth=1)
ggsave('revision round 2/parallel_coord_pca_r1.jpeg', g, width = 12, height = 6)

g = parallel_coord(pca_projected, diatom.df$Depth, diatom.info[as.character(colnames(diatom.X)),'code']) +
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
  theme(axis.text.x = element_text(color = diatom.info[as.character(colnames(diatom.X)),'color'])) +
  geom_hline(yintercept = 0, col='black', linewidth=1)
ggsave('revision round 2/parallel_coord_pca_projected_r1.jpeg', g, width = 12, height = 6)


g = parallel_coord(diatom.res$psao$pts.approx$`r=2`, diatom.df$Depth, diatom.info[as.character(colnames(diatom.X)),'code']) +
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
ggsave('revision round 2/parallel_coord_psao_r1.jpeg', g, width = 12, height = 6)

g = parallel_coord(diatom.res$psas$pts.approx$`r=2`, diatom.df$Depth, diatom.info[as.character(colnames(diatom.X)),'code']) +
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
ggsave('revision round 2/parallel_coord_psas_r1.jpeg', g, width = 12, height = 6)

g = parallel_coord(diatom.res$pca$pts.approx$`r=2`, diatom.res$pca$scores[,1], diatom.info[as.character(colnames(diatom.X)),'code']) +
  scale_color_gradient2(high = 'red', mid = 'white', low = 'blue',
                        name = 'Depth',
                        guide = guide_colorbar(direction = "vertical", reverse = TRUE)) +
  labs(x = '', y = 'Proportion') +
  theme(legend.position = 'right', legend.key.size = unit(10, 'points'), legend.title = element_text(size = 10)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.key.size = unit(12, 'points'),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.spacing.y = unit(1,'points')) +
  theme(axis.text.x = element_text(color = diatom.info[as.character(colnames(diatom.X)),'color'])) +
  scale_y_continuous(limits=c(-0.1,0.2))
ggsave('revision round 2/parallel_coord_pca_r1_rb.jpeg', g, width = 12, height = 6)

g = parallel_coord(pca_projected, diatom.res$pca$scores[,1], diatom.info[as.character(colnames(diatom.X)),'code']) +
  scale_color_gradient2(high = 'red', mid = 'white', low = 'blue',
                        name = 'Depth',
                        guide = guide_colorbar(direction = "vertical", reverse = TRUE)) +
  labs(x = '', y = 'Proportion') +
  theme(legend.position = 'right', legend.key.size = unit(10, 'points'), legend.title = element_text(size = 10)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.key.size = unit(12, 'points'),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.spacing.y = unit(1,'points')) +
  theme(axis.text.x = element_text(color = diatom.info[as.character(colnames(diatom.X)),'color'])) +
  scale_y_continuous(limits=c(-0.1,0.2))
ggsave('revision round 2/parallel_coord_pca_projected_r1_rb.jpeg', g, width = 12, height = 6)


data.frame(Method = 'True', Depth = diatom.df$Depth, value = diatom.df$Thalassiosira.insigna) %>%
  rbind(data.frame(Method = 'PCA', Depth = diatom.df$Depth, value = diatom.res$pca$pts.approx$`r=2`[,1])) %>%
  rbind(data.frame(Method = 'PCA_projected', Depth = diatom.df$Depth, value = pca_projected[,1])) %>%
  rbind(data.frame(Method = 'PSAO', Depth = diatom.df$Depth, value = diatom.res$psao$pts.approx$`r=2`[,1])) %>%
  mutate(Method = factor(Method, levels=c('True','PCA','PCA_projected','PSAO'))) %>%
  ggplot(aes(x=Depth, y = value, group = Method, col = Method)) +
  theme_bw() +
  geom_line() +
  scale_color_manual(values=c('black','red','green','cyan')) +
  labs(y = 'Rank 1 approximation', title = 'T. insigna')
ggsave('revision round 2/approximation T insigna.jpeg', width = 8, height = 4)

data.frame(Method = 'True', Depth = diatom.df$Depth, value = diatom.df$Fragilariopsis.barronii..) %>%
  rbind(data.frame(Method = 'PCA', Depth = diatom.df$Depth, value = diatom.res$pca$pts.approx$`r=2`[,2])) %>%
  rbind(data.frame(Method = 'PCA_projected', Depth = diatom.df$Depth, value = pca_projected[,2])) %>%
  rbind(data.frame(Method = 'PSAO', Depth = diatom.df$Depth, value = diatom.res$psao$pts.approx$`r=2`[,2])) %>%
  mutate(Method = factor(Method, levels=c('True','PCA','PCA_projected','PSAO'))) %>%
  ggplot(aes(x=Depth, y = value, group = Method, col = Method)) +
  theme_bw() +
  geom_line() +
  scale_color_manual(values=c('black','red','green','cyan')) +
  labs(y = 'Rank 1 approximation', title = 'F. barronii..')
ggsave('revision round 2/approximation F barronii.jpeg', width = 8, height = 4)
