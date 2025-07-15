library(dplyr)
diatom.df %>%
  mutate(Depth = rank(Depth)) %>%
  melt(id.vars = 'Depth') %>%
  mutate(variable = factor(variable, levels = colnames(diatom.X))) %>%
  ggplot(aes(x = variable, y = Depth, fill = value)) +
  theme_bw() +
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_fill_gradientn(colors = c('darkgray','cyan','purple','orange','red')) +
  scale_x_discrete(labels = diatom.info[as.character(colnames(diatom.X)),'code']) +
  theme(axis.text.x = element_text(color = diatom.info[as.character(colnames(diatom.X)),'color']))
ggsave('Diatom/diatom_distribution.jpeg', width = 12, height = 6)

diatom.df %>%
  mutate(Depth = rank(Depth)) %>%
  melt(id.vars = 'Depth') %>%
  mutate(variable = factor(variable, levels = colnames(diatom.X))) %>%
  ggplot(aes(x = variable, y = Depth, fill = value)) +
  theme_bw() +
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_fill_gradient(high = 'white', low = 'black') +
  scale_x_discrete(labels = diatom.info[as.character(colnames(diatom.X)),'code']) +
  theme(axis.text.x = element_text(color = diatom.info[as.character(colnames(diatom.X)),'color']))
ggsave('Diatom/diatom_distribution_bw.jpeg', width = 12, height = 6)

diatom.ticks = as.numeric(round(quantile(1:nrow(diatom.df))))

heatmap.df = as.data.frame(diatom.X.log) %>%
  cbind(data.frame(Depth = diatom.df$Depth)) %>%
  mutate(Depth_idx = rank(-Depth)) %>%
  select(-Depth) %>%
  melt(id.vars = 'Depth_idx') %>%
  mutate(variable = factor(variable, levels = colnames(diatom.X)))

g2 = ggplot(heatmap.df, aes(x = variable, y = Depth_idx, fill = log10(value))) +
  theme_minimal() +
  geom_tile() +
  scale_fill_gradientn(colors = c('white','black','black'),
                       values = rescale(c(max(log10(diatom.X.log)), -1.962437, min(log10(diatom.X.log)))),
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

plot_grid(g1, g2, nrow = 1, rel_widths = c(1,10), align = 'h', axis = 'tb')
ggsave('Diatom/diatom_heatmap_bw_log10.jpeg', width = 12, height = 4)
