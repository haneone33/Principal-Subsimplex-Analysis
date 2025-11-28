data(Hydrochem)
hc = Hydrochem
head(hc)
table(hc$River)
hc.label = hc$River
levels(hc.label) = list('Anoia' = 'Anoia', 'Cardener' = 'Cardner',
                        'Lower Llobregat' = 'LowerLLobregat', 'Upper Llobregat' = 'UpperLLobregat')
hc.X = to_simplex(as.matrix(hc[,6:19]))

hc.pca = comp_pca(hc.X)
hc.pca.df = data.frame(label = hc.label) %>%
  cbind(as.data.frame(hc.pca$scores))
ggplot(hc.pca.df, aes(x=PC1, y=PC2, col = hc.label)) +
  theme_bw() +
  geom_point()

hc.power_pca = comp_power_pca(hc.X)
hc.power_pca.df = data.frame(label = hc.label) %>%
  cbind(as.data.frame(hc.power_pca$scores))
ggplot(hc.power_pca.df, aes(x=PC1, y=PC2, col = hc.label)) +
  theme_bw() +
  geom_point()

hc.apca = comp_apca(hc.X)
hc.apca.df = data.frame(label = hc.label) %>%
  cbind(as.data.frame(hc.apca$scores))
ggplot(hc.apca.df, aes(x=PC1, y=PC2, col = hc.label)) +
  theme_bw() +
  geom_point()

hc.psas = psa('s', hc.X)
hc.psas.df = data.frame(label = hc.label) %>%
  cbind(as.data.frame(hc.psas$scores))
ggplot(hc.psas.df, aes(x=PC1, y=PC2, col = hc.label)) +
  theme_bw() +
  geom_point()

hc.psao = psa('o', hc.X)
hc.psao.df = data.frame(label = hc.label) %>%
  cbind(as.data.frame(hc.psao$scores))
ggplot(hc.psao.df, aes(x=PC1, y=PC2, col = hc.label)) +
  theme_bw() +
  geom_point()

L = 8
hc.pca.accuracy = rep(NA, L)
for(i in 1:L){
  qda.fit <- qda(label ~ ., data = hc.pca.df[,1:(i+1)])
  qda.prediction <- predict(qda.fit, newdata = hc.pca.df)
  hc.pca.accuracy[i] = mean(qda.prediction$class == hc.label)
}

hc.psas.accuracy = rep(NA, L)
for(i in 1:L){
  # df = as.data.frame(hc.psas$Xhat_reduced[[i+1]]) %>%
  #   mutate(label = hc.label)
  # qda.fit <- qda(label ~ ., data = df)
  # qda.prediction <- predict(qda.fit, newdata = df)

  rda.fit = rda(hc.psas$Xhat_reduced[[i+1]], grouping = hc.label,
                lambda = 0.2, gamma = 0.2)
  rda.prediction = predict(rda.fit, hc.psas$Xhat_reduced[[i+1]])

  # qda.fit <- qda(label ~ ., data = hc.psas.df[,1:(i+1)])
  # qda.prediction <- predict(qda.fit, newdata = hc.psas.df)
  hc.psas.accuracy[i] = mean(rda.prediction$class == hc.label)
}

hc.psao.accuracy = rep(NA, L)
for(i in 1:L){
  # df = as.data.frame(hc.psao$Xhat_reduced[[i+1]]) %>%
  #   mutate(label = hc.label)
  # qda.fit <- qda(label ~ ., data = df)
  # qda.prediction <- predict(qda.fit, newdata = df)

  rda.fit = rda(hc.psao$Xhat_reduced[[i+1]], grouping = hc.label,
                lambda = 0.2, gamma = 0.2)
  rda.prediction = predict(rda.fit, hc.psao$Xhat_reduced[[i+1]])

  # qda.fit <- qda(label ~ ., data = hc.psao.df[,1:(i+1)])
  # qda.prediction <- predict(qda.fit, newdata = hc.psao.df)
  hc.psao.accuracy[i] = mean(rda.prediction$class == hc.label)
}

hc.apca.accuracy = rep(NA, L)
for(i in 1:L){
  qda.fit <- qda(label ~ ., data = hc.apca.df[,1:(i+1)])
  qda.prediction <- predict(qda.fit, newdata = hc.apca.df)
  hc.apca.accuracy[i] = mean(qda.prediction$class == hc.label)
}

accuracy = as.data.frame(rbind(hc.pca.accuracy,
                               hc.psas.accuracy,
                               hc.psao.accuracy,
                               hc.apca.accuracy)) %>%
  mutate(Method = factor(c('PCA','PSAS','PSAO','Log-ratio'), levels = c('PCA','Log-ratio','PSAS','PSAO')))

accuracy %>%
  melt(id.vars = 'Method') %>%
  mutate(variable = as.numeric(as.factor(variable))) %>%
  ggplot(aes(x=variable, y=value, group = Method, color = Method)) +
  theme_bw() +
  geom_line() +
  geom_point() +
  ggtitle('Hydrochemical, QDA') +
  scale_x_continuous(name = 'Rank', breaks = 1:L) +
  scale_y_continuous(name = 'Accuracy', limits = c(0,1))
ggsave('revision round 2/cf_hydrochemical_qda.jpeg', width = 8, height = 4)

L = 10
hc.gamma = 2
ctrl = trainControl(method = 'cv', number = 5)
hc.pca.accuracy.svm = rep(NA, L)
for(i in 1:L){
  svm.fit <- svm(label ~ ., data = hc.pca.df[,1:(i+1)],
                 kernel = "linear", cost = 1, gamma = hc.gamma)
  svm.prediction <- predict(svm.fit, newdata = hc.pca.df)
  hc.pca.accuracy.svm[i] = mean(svm.prediction == hc.label)
}


hc.psas.accuracy.svm = rep(NA, L)
for(i in 1:L){
  df = as.data.frame(hc.psas$Xhat_reduced[[i+1]]) %>%
    mutate(label = hc.label)
  svm.fit <- svm(label ~ ., data = df,
                 kernel = "linear", cost = 1, gamma = hc.gamma)
  svm.prediction <- predict(svm.fit, newdata = df)

  # svm.fit <- svm(label ~ ., data = hc.psas.df[,1:(i+1)],
  #                kernel = "linear", cost = 1)
  # svm.prediction <- predict(svm.fit, newdata = hc.psas.df)

  hc.psas.accuracy.svm[i] = mean(svm.prediction == hc.label)
}

hc.psao.accuracy.svm = rep(NA, L)
for(i in 1:L){
  df = as.data.frame(hc.psao$Xhat_reduced[[i+1]]) %>%
    mutate(label = hc.label)
  svm.fit <- svm(label ~ ., data = df,
                 kernel = "linear", cost = 1, gamma = hc.gamma)
  svm.prediction <- predict(svm.fit, newdata = df)

  # svm.fit <- svm(label ~ ., data = hc.psao.df[,1:(i+1)],
  #                kernel = "linear", cost = 1)
  # svm.prediction <- predict(svm.fit, newdata = hc.psao.df)

  hc.psao.accuracy.svm[i] = mean(svm.prediction == hc.label)
}

hc.apca.accuracy.svm = rep(NA, L)
for(i in 1:L){
  svm.fit <- svm(label ~ ., data = hc.apca.df[,1:(i+1)],
                 kernel = "linear", cost = 1, gamma = c(0.1,0.2))
  svm.prediction <- predict(svm.fit, newdata = hc.apca.df)
  hc.apca.accuracy.svm[i] = mean(svm.prediction == hc.label)
}

################################################################################
set.seed(1)
hc.pca.accuracy.svm = rep(NA, L)
for(i in 1:L){
  svm.fit = tune.svm(label ~., data = hc.pca.df[,1:(i+1)],
                     kernel = 'radial',
                     gamma = 2**seq(-4,4,length.out=17),
                     cost = 1,
                     tunecontrol = tune.control(sampling = "cross", cross = 5))

  hc.pca.accuracy.svm[i] = 1-svm.fit$best.performance
}

hc.power_pca.accuracy.svm = rep(NA, L)
for(i in 1:L){
  svm.fit = tune.svm(label ~., data = hc.power_pca.df[,1:(i+1)],
                     kernel = 'radial',
                     gamma = 2**seq(-4,4,length.out=17),
                     cost = 1,
                     tunecontrol = tune.control(sampling = "cross", cross = 5))

  hc.power_pca.accuracy.svm[i] = 1-svm.fit$best.performance
}

hc.psas.accuracy.svm = rep(NA, L)
for(i in 1:L){
#  df = as.data.frame(hc.psas$Xhat_reduced[[i+1]]) %>%
  df = as.data.frame(hc.psas$scores[,1:i]) %>%
    mutate(label = hc.label)
  svm.fit = tune.svm(label ~., data = df,
                     kernel = 'radial',
                     gamma = 2**seq(-4,4,length.out=17),
                     cost = 1,
                     tunecontrol = tune.control(sampling = "cross", cross = 5))

  hc.psas.accuracy.svm[i] = 1-svm.fit$best.performance
}

hc.psao.accuracy.svm = rep(NA, L)
for(i in 1:L){
  # df = as.data.frame(hc.psao$Xhat_reduced[[i+1]]) %>%
  df = as.data.frame(hc.psao$scores[,1:i]) %>%
    mutate(label = hc.label)
  svm.fit = tune.svm(label ~., data = df,
                     kernel = 'radial',
                     gamma = 2**seq(-4,4,length.out=17),
                     cost = 1,
                     tunecontrol = tune.control(sampling = "cross", cross = 5))

  hc.psao.accuracy.svm[i] = 1-svm.fit$best.performance
}

hc.apca.accuracy.svm = rep(NA, L)
for(i in 1:L){
  svm.fit = tune.svm(label ~., data = hc.apca.df[,1:(i+1)],
                     kernel = 'radial',
                     gamma = 2**seq(-4,4,length.out=17),
                     cost = 1,
                     tunecontrol = tune.control(sampling = "cross", cross = 5))

  hc.apca.accuracy.svm[i] = 1-svm.fit$best.performance
}
################################################################################

accuracy.svm = as.data.frame(rbind(hc.pca.accuracy.svm,
                                   hc.power_pca.accuracy.svm,
                                   hc.psas.accuracy.svm,
                                   hc.psao.accuracy.svm,
                                   hc.apca.accuracy.svm)) %>%
  mutate(Method = factor(c('PCA','Power transform','PSAS','PSAO','Log-ratio'),
                         levels = c('PSAS','PSAO','PCA','Power transform','Log-ratio')))

accuracy.svm %>%
  melt(id.vars = 'Method') %>%
  mutate(variable = as.numeric(as.factor(variable))) %>%
  ggplot(aes(x=variable, y=value, group = Method, color = Method)) +
  theme_bw() +
  geom_line() +
  geom_point() +
#  ggtitle('Hydrochemical, SVM') +
  scale_x_continuous(name = 'Rank', breaks = 1:L) +
  scale_y_continuous(name = 'Accuracy', limits = c(0.4,1), breaks = c(0,0.2,0.4,0.6,0.8,1)) +
  scale_color_manual(values = c('red','orange','blue','dodgerblue','darkgray'))
ggsave('revision round 2/cf_hydrochemical_svm.jpeg', width = 8, height = 4)


as.data.frame(hc.psas$Xhat_reduced$`r=2`) %>%
  mutate(label = hc.label ) %>%
  ggtern(aes(x=V1, y=V2, z=V3, col = label)) +
  geom_point()

parallel_coord(hc.X, hc.label) +
  theme(legend.position='right') +
  scale_color_discrete(name = 'Tributary')
ggsave('revision round 2/hydrochemical_parallel_coord.jpeg', width = 8, height = 4)
