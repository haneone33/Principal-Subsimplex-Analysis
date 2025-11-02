diatom.df[,c('Thalassiosira.insigna',
             'Thalassiosira.inura',
             'Fragilariopsis.barronii..',
             'Fragilariopsis.bohatyi.',
             'Chaetoceros.spp...resting.spores.',
             'Dactyliosolen.antarcticus.',
             'Depth')]

diatom.X[,1:9] %>%
  mutate(Depth = diatom.df$Depth) %>%
  melt(id.vars = 'Depth') %>%
  ggplot(aes(x=Depth, y=value)) +
  theme_bw() +
  geom_line(linewidth=1) +
  facet_wrap(vars(variable), ncol=3)
