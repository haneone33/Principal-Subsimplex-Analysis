X22 = X2
X22[X22==0] = min(X22[X22>0])/2
X2.apca = princomp(clr(X22))

point_gpairs(X2.apca$scores, ex.label)
ggsave('Simulation/ex2_replace_by_overall_minimum.jpeg', width = 9, height = 6)
point_gpairs(ex2.res$apca$scores, ex.label)
ggsave('Simulation/ex2_replace_by_columnwise_minimum.jpeg', width = 9, height = 6)

