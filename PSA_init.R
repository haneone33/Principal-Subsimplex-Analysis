tmp = lapply(list.files('PSA', pattern = '.R',
                        full.names = T, recursive = T), source)
tmp = lapply(list.files('Functions', pattern = '.R',
                        full.names = T, recursive = T), source)
rm(tmp)

require(compositions)
require(MASS)
require(ggplot2)
require(ggtern)
require(latex2exp)
require(GGally)
require(grid)
require(gridExtra)
require(cowplot)
require(reshape)
require(dendextend)
