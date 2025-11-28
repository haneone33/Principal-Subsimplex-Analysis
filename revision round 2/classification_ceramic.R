ceramic = read.csv("revision round 2/Chemical Composion of Ceramic.csv")
ceramic = ceramic[-2,]
ceramic.ls = strsplit(ceramic$Ceramic.Name, '-')
ceramic.label = rep(NA, length(ceramic.ls))
for(i in 1:length(ceramic.ls)){
  x = ceramic.ls[[i]]
  if(x[1] == 'FLQ'){
    ceramic.label[i] = paste0('FLQ','_',x[3])
  }else{
    ceramic.label[i] = paste0(x[1],'_',x[4])
  }
}
ceramic.label = as.factor(ceramic.label)
table(ceramic.label)
ceramic.X = as.matrix(ceramic[,3:19])
ceramic.X[,1:8] = ceramic.X[,1:8]/100
ceramic.X[,9:17] = ceramic.X[,9:17]/1000000
ceramic.X = to_simplex(ceramic.X)

ceramic.pca = comp_pca(ceramic.X)
ceramic.pca.df = data.frame(label = ceramic.label) %>%
  cbind(as.data.frame(ceramic.pca$scores))
ggplot(ceramic.pca.df, aes(x=PC1, y=PC2, col = ceramic.label)) +
  theme_bw() +
  geom_point()

ceramic.apca = comp_apca(ceramic.X)
ceramic.apca.df = data.frame(label = ceramic.label) %>%
  cbind(as.data.frame(ceramic.apca$scores))
ggplot(ceramic.apca.df, aes(x=PC1, y=PC2, col = ceramic.label)) +
  theme_bw() +
  geom_point()

ceramic.psas = psa('s', ceramic.X)
ceramic.psas.df = data.frame(label = ceramic.label) %>%
  cbind(as.data.frame(ceramic.psas$scores))
ggplot(ceramic.psas.df, aes(x=PC1, y=PC2, col = ceramic.label)) +
  theme_bw() +
  geom_point()

ceramic.psao = psa('o', ceramic.X)
ceramic.psao.df = data.frame(label = ceramic.label) %>%
  cbind(as.data.frame(ceramic.psao$scores))
ggplot(ceramic.psao.df, aes(x=PC1, y=PC2, col = ceramic.label)) +
  theme_bw() +
  geom_point()

L = 10
ceramic.pca.accuracy = rep(NA, L)
for(i in 1:L){
  qda.fit <- qda(label ~ ., data = ceramic.pca.df[,1:(i+1)])
  qda.prediction <- predict(qda.fit, newdata = ceramic.pca.df)
  ceramic.pca.accuracy[i] = mean(qda.prediction$class == ceramic.label)
}

ceramic.psas.accuracy = rep(NA, L)
for(i in 1:L){
  qda.fit <- qda(label ~ ., data = ceramic.psas.df[,1:(i+1)])
  qda.prediction <- predict(qda.fit, newdata = ceramic.psas.df)
  ceramic.psas.accuracy[i] = mean(qda.prediction$class == ceramic.label)
}

ceramic.psao.accuracy = rep(NA, L)
for(i in 1:L){
  qda.fit <- qda(label ~ ., data = ceramic.psao.df[,1:(i+1)])
  qda.prediction <- predict(qda.fit, newdata = ceramic.psao.df)
  ceramic.psao.accuracy[i] = mean(qda.prediction$class == ceramic.label)
}

ceramic.apca.accuracy = rep(NA, L)
for(i in 1:L){
  qda.fit <- qda(label ~ ., data = ceramic.apca.df[,1:(i+1)])
  qda.prediction <- predict(qda.fit, newdata = ceramic.apca.df)
  ceramic.apca.accuracy[i] = mean(qda.prediction$class == ceramic.label)
}

accuracy = as.data.frame(rbind(ceramic.pca.accuracy,
                               ceramic.psas.accuracy,
                               ceramic.psao.accuracy,
                               ceramic.apca.accuracy)) %>%
  mutate(Method = factor(c('PCA','PSAS','PSAO','Log-ratio'), levels = c('PCA','Log-ratio','PSAS','PSAO')))

accuracy %>%
  melt(id.vars = 'Method') %>%
  mutate(variable = as.numeric(as.factor(variable))) %>%
  ggplot(aes(x=variable, y=value, group = Method, color = Method)) +
  theme_bw() +
  geom_line() +
  geom_point() +
  ggtitle('ceramic, QDA') +
  scale_x_continuous(name = 'Rank', breaks = 1:L) +
  scale_y_continuous(name = 'Accuracy', limits = c(0.5,1))
ggsave('revision round 2/cf_ceramic_qda.jpeg', width = 8, height = 4)


L = 10
ceramic.pca.accuracy.svm = rep(NA, L)
for(i in 1:L){
  svm.fit <- svm(label ~ ., data = ceramic.pca.df[,1:(i+1)],
                 kernel = "radial", cost = 1, gamma = 0.1)
  svm.prediction <- predict(svm.fit, newdata = ceramic.pca.df)
  ceramic.pca.accuracy.svm[i] = mean(svm.prediction == ceramic.label)
}

ceramic.psas.accuracy.svm = rep(NA, L)
for(i in 1:L){
  svm.fit <- svm(label ~ ., data = ceramic.psas.df[,1:(i+1)],
                 kernel = "radial", cost = 1, gamma = 0.1)
  svm.prediction <- predict(svm.fit, newdata = ceramic.psas.df)
  ceramic.psas.accuracy.svm[i] = mean(svm.prediction == ceramic.label)
}

ceramic.psao.accuracy.svm = rep(NA, L)
for(i in 1:L){
  svm.fit <- svm(label ~ ., data = ceramic.psao.df[,1:(i+1)],
                 kernel = "radial", cost = 1, gamma = 0.1)
  svm.prediction <- predict(svm.fit, newdata = ceramic.psao.df)
  ceramic.psao.accuracy.svm[i] = mean(svm.prediction == ceramic.label)
}

ceramic.apca.accuracy.svm = rep(NA, L)
for(i in 1:L){
  svm.fit <- svm(label ~ ., data = ceramic.apca.df[,1:(i+1)],
                 kernel = "radial", cost = 1, gamma = 0.1)
  svm.prediction <- predict(svm.fit, newdata = ceramic.apca.df)
  ceramic.apca.accuracy.svm[i] = mean(svm.prediction == ceramic.label)
}

accuracy.svm = as.data.frame(rbind(ceramic.pca.accuracy.svm,
                                   ceramic.psas.accuracy.svm,
                                   ceramic.psao.accuracy.svm,
                                   ceramic.apca.accuracy.svm)) %>%
  mutate(Method = factor(c('PCA','PSAS','PSAO','Log-ratio'), levels = c('PCA','Log-ratio','PSAS','PSAO')))

accuracy.svm %>%
  melt(id.vars = 'Method') %>%
  mutate(variable = as.numeric(as.factor(variable))) %>%
  ggplot(aes(x=variable, y=value, group = Method, color = Method)) +
  theme_bw() +
  geom_line() +
  geom_point() +
  ggtitle('ceramic, SVM') +
  scale_x_continuous(name = 'Rank', breaks = 1:L) +
  scale_y_continuous(name = 'Accuracy', limits = c(0.5,1))
ggsave('revision round 2/cf_ceramic_svm.jpeg', width = 8, height = 4)

