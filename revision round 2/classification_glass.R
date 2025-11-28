glass = read.table("revision round 2/glass.data", sep=',')
glass.label = as.factor(glass[,11])
glass.X = to_simplex(as.matrix(glass[glass.label != '6',3:10]))
colnames(glass.X) = c('Na','Mg','Al','Si','K','Ca','Ba','Fe')
glass.label = factor(glass.label[glass.label != '6'], levels = c('1','2','3','5','7'))


glass.pca = comp_pca(glass.X)
glass.pca.df = data.frame(label = glass.label) %>%
  cbind(as.data.frame(glass.pca$scores))
ggplot(glass.pca.df, aes(x=PC1, y=PC2, col = glass.label)) +
  theme_bw() +
  geom_point()

glass.apca = comp_apca(glass.X)
glass.apca.df = data.frame(label = glass.label) %>%
  cbind(as.data.frame(glass.apca$scores))
ggplot(glass.apca.df, aes(x=PC1, y=PC2, col = glass.label)) +
  theme_bw() +
  geom_point()

glass.psas = psa('s', glass.X)
glass.psas.df = data.frame(label = glass.label) %>%
  cbind(as.data.frame(glass.psas$scores))
ggplot(glass.psas.df, aes(x=PC1, y=PC2, col = glass.label)) +
  theme_bw() +
  geom_point()

glass.psao = psa('o', glass.X)
glass.psao.df = data.frame(label = glass.label) %>%
  cbind(as.data.frame(glass.psao$scores))
ggplot(glass.psao.df, aes(x=PC1, y=PC2, col = glass.label)) +
  theme_bw() +
  geom_point()

L = 7
glass.pca.accuracy = rep(NA, L)
for(i in 1:L){
  qda.fit <- qda(label ~ ., data = glass.pca.df[,1:(i+1)])
  qda.prediction <- predict(qda.fit, newdata = glass.pca.df)
  glass.pca.accuracy[i] = mean(qda.prediction$class == glass.label)
}

glass.psas.accuracy = rep(NA, L)
for(i in 1:L){
  qda.fit <- qda(label ~ ., data = glass.psas.df[,1:(i+1)])
  qda.prediction <- predict(qda.fit, newdata = glass.psas.df)
  glass.psas.accuracy[i] = mean(qda.prediction$class == glass.label)
}

glass.psao.accuracy = rep(NA, L)
for(i in 1:L){
  qda.fit <- qda(label ~ ., data = glass.psao.df[,1:(i+1)])
  qda.prediction <- predict(qda.fit, newdata = glass.psao.df)
  glass.psao.accuracy[i] = mean(qda.prediction$class == glass.label)
}

glass.apca.accuracy = rep(NA, L)
for(i in 1:L){
  qda.fit <- qda(label ~ ., data = glass.apca.df[,1:(i+1)])
  qda.prediction <- predict(qda.fit, newdata = glass.apca.df)
  glass.apca.accuracy[i] = mean(qda.prediction$class == glass.label)
}

accuracy = as.data.frame(rbind(glass.pca.accuracy,
                               glass.psas.accuracy,
                               glass.psao.accuracy,
                               glass.apca.accuracy)) %>%
  mutate(Method = factor(c('PCA','PSAS','PSAO','Log-ratio'), levels = c('PCA','Log-ratio','PSAS','PSAO')))

accuracy %>%
  melt(id.vars = 'Method') %>%
  mutate(variable = as.numeric(as.factor(variable))) %>%
  ggplot(aes(x=variable, y=value, group = Method, color = Method)) +
  theme_bw() +
  geom_line() +
  geom_point() +
  ggtitle('Glass, QDA') +
  scale_x_continuous(name = 'Rank', breaks = 1:L) +
  scale_y_continuous(name = 'Accuracy', limits = c(0,1))
ggsave('revision round 2/cf_glass_qda.jpeg', width = 8, height = 4)


L = 7
glass.pca.accuracy.svm = rep(NA, L)
for(i in 1:L){
  svm.fit <- svm(label ~ ., data = glass.pca.df[,1:(i+1)],
                 kernel = "radial", cost = 1, gamma = 0.1)
  svm.prediction <- predict(svm.fit, newdata = glass.pca.df)
  glass.pca.accuracy.svm[i] = mean(svm.prediction == glass.label)
}

glass.psas.accuracy.svm = rep(NA, L)
for(i in 1:L){
  svm.fit <- svm(label ~ ., data = glass.psas.df[,1:(i+1)],
                 kernel = "radial", cost = 1, gamma = 0.1)
  svm.prediction <- predict(svm.fit, newdata = glass.psas.df)
  glass.psas.accuracy.svm[i] = mean(svm.prediction == glass.label)
}

glass.psao.accuracy.svm = rep(NA, L)
for(i in 1:L){
  svm.fit <- svm(label ~ ., data = glass.psao.df[,1:(i+1)],
                 kernel = "radial", cost = 1, gamma = 0.1)
  svm.prediction <- predict(svm.fit, newdata = glass.psao.df)
  glass.psao.accuracy.svm[i] = mean(svm.prediction == glass.label)
}

glass.apca.accuracy.svm = rep(NA, L)
for(i in 1:L){
  svm.fit <- svm(label ~ ., data = glass.apca.df[,1:(i+1)],
                 kernel = "radial", cost = 1, gamma = 0.1)
  svm.prediction <- predict(svm.fit, newdata = glass.apca.df)
  glass.apca.accuracy.svm[i] = mean(svm.prediction == glass.label)
}

accuracy.svm = as.data.frame(rbind(glass.pca.accuracy.svm,
                                   glass.psas.accuracy.svm,
                                   glass.psao.accuracy.svm,
                                   glass.apca.accuracy.svm)) %>%
  mutate(Method = factor(c('PCA','PSAS','PSAO','Log-ratio'), levels = c('PCA','Log-ratio','PSAS','PSAO')))

accuracy.svm %>%
  melt(id.vars = 'Method') %>%
  mutate(variable = as.numeric(as.factor(variable))) %>%
  ggplot(aes(x=variable, y=value, group = Method, color = Method)) +
  theme_bw() +
  geom_line() +
  geom_point() +
  ggtitle('Glass, SVM') +
  scale_x_continuous(name = 'Rank', breaks = 1:L) +
  scale_y_continuous(name = 'Accuracy', limits = c(0,1))
ggsave('revision round 2/cf_glass_svm.jpeg', width = 8, height = 4)

