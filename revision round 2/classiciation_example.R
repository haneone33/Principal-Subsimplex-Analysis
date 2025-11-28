tgrid = seq(0,1,by=0.1)
tt = cbind(tgrid, 1-tgrid)

################################################################################

X2 = tt %*% matrix(c(0, 1, 0, 0,
                     0, 0, 1, 0), byrow = T, ncol = 4)
X3 = tt %*% matrix(c(0, 0.9, 0, 0.1,
                     0, 0, 0.9, 0.1), byrow = T, ncol = 4)
X4 = tt %*% matrix(c(0, 0.8, 0, 0.2,
                     0, 0, 0.8, 0.2), byrow = T, ncol = 4)
X5 = tt %*% matrix(c(0, 0.7, 0, 0.3,
                     0, 0, 0.7, 0.3), byrow = T, ncol = 4)
n = 30
X1 = matrix(rep(c(0.2,0.2,0.2,0.4),n), byrow = T, ncol = 4)
X = rbind(X1, X2, X3, X4, X5)
colnames(X) = paste0('e',1:4)
x.label = factor(c(rep(1,n), rep(2,11), rep(3,11), rep(4,11), rep(5,11)))

res = comp_pca(X)
to_simplex(res$Xhat$`r=1`[c(0,11,22,33,44)+n,])

res2 = psa('s', X)
res2$Vhat$`r=2`
res2$Xhat$`r=1`[c(0,11,22,33,44)+n,]


################################################################################
X = rbind(tt %*% matrix(c(0, 1, 0,
                          0, 0.7, 0.3), byrow = T, ncol = 3),
          tt %*% matrix(c(0.1, 0.9, 0,
                          0.1, 0.6, 0.3), byrow = T, ncol = 3),
          tt %*% matrix(c(0.2, 0.8, 0,
                          0.2, 0.5, 0.3), byrow = T, ncol = 3),
          tt %*% matrix(c(0.3, 0.7, 0,
                          0.3, 0.4, 0.3), byrow = T, ncol = 3),
          tt %*% matrix(c(0.4, 0.6, 0,
                          0.4, 0.3, 0.3), byrow = T, ncol = 3),
          tt %*% matrix(c(0.5, 0.5, 0,
                          0.5, 0.2, 0.3), byrow = T, ncol = 3),
          tt %*% matrix(c(0.6, 0.4, 0,
                          0.6, 0.1, 0.3), byrow = T, ncol = 3),
          tt %*% matrix(c(0.7, 0.3, 0,
                          0.7, 0, 0.3), byrow = T, ncol = 3),
          matrix(rep(c(1,0,0),10),byrow = T, ncol = 3))
X.label = factor(c(rep(1:8, each = 11), rep(9,10)))
X.pca = comp_pca(X)
X.psas = psa('s', X)
X.psao = psa('o', X)
X.apca = comp_apca(X)
g = plot_grid(ternary_pc(X, X.label, type = 'pca', X.pca),
              ternary_pc(X, X.label, type = 'logratio', X.apca),
              ternary_pc(X, X.label, type = 'psas', X.psas),
              ternary_pc(X, X.label, type = 'psao', X.psao),
              nrow = 2)
ggsave('revision round 2/classification_example1.jpeg', g, width = 12, height = 12)

X = rbind(tt %*% matrix(c(0, 1, 0,
                          0, 0.7, 0.3), byrow = T, ncol = 3),
          tt %*% matrix(c(1/7, 6/7, 0,
                          0.1, 0.6, 0.3), byrow = T, ncol = 3),
          tt %*% matrix(c(2/7, 5/7, 0,
                          0.2, 0.5, 0.3), byrow = T, ncol = 3),
          tt %*% matrix(c(3/7, 4/7, 0,
                          0.3, 0.4, 0.3), byrow = T, ncol = 3),
          tt %*% matrix(c(4/7, 3/7, 0,
                          0.4, 0.3, 0.3), byrow = T, ncol = 3),
          tt %*% matrix(c(5/7, 2/7, 0,
                          0.5, 0.2, 0.3), byrow = T, ncol = 3),
          tt %*% matrix(c(6/7, 1/7, 0,
                          0.6, 0.1, 0.3), byrow = T, ncol = 3),
          tt %*% matrix(c(1, 0, 0,
                          0.7, 0, 0.3), byrow = T, ncol = 3))
X.label = factor(c(rep(1:8, each = 11)))
X.pca = comp_pca(X)
X.psas = psa('s', X)
X.psao = psa('o', X)
X.apca = comp_apca(X)
g = plot_grid(ternary_pc(X, X.label, type = 'pca', X.pca),
              ternary_pc(X, X.label, type = 'logratio', X.apca),
              ternary_pc(X, X.label, type = 'psas', X.psas),
              ternary_pc(X, X.label, type = 'psao', X.psao),
              nrow = 2)
ggsave('revision round 2/classification_example2.jpeg', g, width = 12, height = 12)

################################################################################

set.seed(1)
X = rbind(make_cluster(c(0.1, 0.7, 0.2), 10, 0.1^2),
          make_cluster(c(0.7, 0.1, 0.1), 10, 0.1^2),
          make_cluster(c(0.4, 0.4, 0.2), 30, 0.2^2))
df = as.data.frame(X) %>%
  mutate(label = factor(ifelse(V1>0.4, '1', '2')))
ggtern(df, aes(x=V1,y=V3,z=V2, col = label)) +
  geom_point()
X.pca = comp_pca(X)
X.psas = psa('s', X)
X.psao = psa('o', X)
X.apca = comp_apca(X)
g = plot_grid(ternary_pc(X, df$label, type = 'pca', X.pca),
              ternary_pc(X, df$label, type = 'logratio', X.apca),
              ternary_pc(X, df$label, type = 'psas', X.psas),
              ternary_pc(X, df$label, type = 'psao', X.psao),
              nrow = 2)
g

################################################################################

set.seed(1)
X = rbind(make_cluster(c(0.1, 0.8, 0.1), 10, 0.1^2),
          make_cluster(c(0.8, 0.1, 0.1), 10, 0.1^2),
          make_cluster(c(0.45, 0.45, 0.1), 30, 0.1^2))
df = as.data.frame(X) %>%
  mutate(label = factor(ifelse(V1>0.45, '1', '2')))
X.pca = comp_pca(X)
X.power_pca = comp_power_pca(X, 0.5)
X.psas = psa('s', X)
X.psao = psa('o', X)
X.apca = comp_apca(X)
g = plot_grid(ternary_pc(X, df$label, type = 'data'),
              ternary_pc(X, df$label, type = 'psas', X.psas) +
                annotate('text',x=0.7,y=1.5,label='error rate: 0.00'),
              ternary_pc(X, df$label, type = 'psao', X.psao) +
                annotate('text',x=0.7,y=1.5,label='error rate: 0.08'),
              ternary_pc(X, df$label, type = 'pca', X.pca) +
                annotate('text',x=0.7,y=1.5,label='error rate: 0.08'),
              ternary_pc(X, df$label, type = 'power', X.power_pca) +
                annotate('text',x=0.7,y=1.5,label='error rate: 0.06'),
              ternary_pc(X, df$label, type = 'logratio', X.apca) +
                annotate('text',x=0.7,y=1.5,label='error rate: 0.20'),
              nrow = 2)
ggsave('revision round 2/classiciation_example.jpeg', g, width = 12, height = 8)

get_error_rate <- function(y_true){
  n = length(y_true)
  error_rate = rep(0., n)
  for(i in 1:n){
    y1 = y_true[1:i]
    y2 = y_true[(i+1):n]
    e = (sum(y1!='1') + sum(y2!='2'))/n
    error_rate[i] = min(e, 1-e)
  }
  return(error_rate)
}


c(data.frame(PC1 = X.psas$scores[,1],
             label = df$label) %>%
    arrange(PC1) %>%
    pull(label) %>%
    get_error_rate() %>%
    min(na.rm = T),
  data.frame(PC1 = X.psao$scores[,1],
             label = df$label) %>%
    arrange(PC1) %>%
    pull(label) %>%
    get_error_rate() %>%
    min(na.rm = T),
  data.frame(PC1 = X.pca$scores[,1],
             label = df$label) %>%
    arrange(PC1) %>%
    pull(label) %>%
    get_error_rate() %>%
    min(na.rm = T),
  data.frame(PC1 = X.power_pca$scores[,1],
             label = df$label) %>%
    arrange(PC1) %>%
    pull(label) %>%
    get_error_rate() %>%
    min(na.rm = T),
  data.frame(PC1 = X.apca$scores[,1],
             label = df$label) %>%
    arrange(PC1) %>%
    pull(label) %>%
    get_error_rate() %>%
    min(na.rm = T))

