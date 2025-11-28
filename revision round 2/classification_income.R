income = read.csv("revision round 2/income.csv", sep=',')
income.label = as.factor(income$Income_category)
income.X = to_simplex(as.matrix(income[,2:6]))
income.X = to_simplex(as.matrix(income[income.label != '5',2:6]))
income.label = factor(income.label[income.label != '5'], levels = c('1','2','3','4'))


income.pca = comp_pca(income.X)
income.pca.df = data.frame(label = income.label) %>%
  cbind(as.data.frame(income.pca$scores))
ggplot(income.pca.df, aes(x=PC1, y=PC2, col = income.label)) +
  theme_bw() +
  geom_point()

income.apca = comp_apca(income.X)
income.apca.df = data.frame(label = income.label) %>%
  cbind(as.data.frame(income.apca$scores))
ggplot(income.apca.df, aes(x=PC1, y=PC2, col = income.label)) +
  theme_bw() +
  geom_point()

income.psas = psa('s', income.X)
income.psas.df = data.frame(label = income.label) %>%
  cbind(as.data.frame(income.psas$scores))
ggplot(income.psas.df, aes(x=PC1, y=PC2, col = income.label)) +
  theme_bw() +
  geom_point()

income.psao = psa('o', income.X)
income.psao.df = data.frame(label = income.label) %>%
  cbind(as.data.frame(income.psao$scores))
ggplot(income.psao.df, aes(x=PC1, y=PC2, col = income.label)) +
  theme_bw() +
  geom_point()

L = 4
income.pca.accuracy = rep(NA, L)
for(i in 1:L){
  qda.fit <- qda(label ~ ., data = income.pca.df[,1:(i+1)])
  qda.prediction <- predict(qda.fit, newdata = income.pca.df)
  income.pca.accuracy[i] = mean(qda.prediction$class == income.label)
}

income.psas.accuracy = rep(NA, L)
for(i in 1:L){
  qda.fit <- qda(label ~ ., data = income.psas.df[,1:(i+1)])
  qda.prediction <- predict(qda.fit, newdata = income.psas.df)
  income.psas.accuracy[i] = mean(qda.prediction$class == income.label)
}

income.psao.accuracy = rep(NA, L)
for(i in 1:L){
  qda.fit <- qda(label ~ ., data = income.psao.df[,1:(i+1)])
  qda.prediction <- predict(qda.fit, newdata = income.psao.df)
  income.psao.accuracy[i] = mean(qda.prediction$class == income.label)
}

income.apca.accuracy = rep(NA, L)
for(i in 1:L){
  qda.fit <- qda(label ~ ., data = income.apca.df[,1:(i+1)])
  qda.prediction <- predict(qda.fit, newdata = income.apca.df)
  income.apca.accuracy[i] = mean(qda.prediction$class == income.label)
}

accuracy = as.data.frame(rbind(income.pca.accuracy,
                               income.psas.accuracy,
                               income.psao.accuracy,
                               income.apca.accuracy)) %>%
  mutate(Method = factor(c('PCA','PSAS','PSAO','Log-ratio'), levels = c('PCA','Log-ratio','PSAS','PSAO')))

accuracy %>%
  melt(id.vars = 'Method') %>%
  mutate(variable = as.numeric(as.factor(variable))) %>%
  ggplot(aes(x=variable, y=value, group = Method, color = Method)) +
  theme_bw() +
  geom_line() +
  geom_point() +
  ggtitle('income, QDA') +
  scale_x_continuous(name = 'Rank', breaks = 1:L) +
  scale_y_continuous(name = 'Accuracy', limits = c(0,1))
ggsave('revision round 2/cf_income_qda.jpeg', width = 8, height = 4)

L = 4
income.pca.accuracy.svm = rep(NA, L)
for(i in 1:L){
  svm.fit <- svm(label ~ ., data = income.pca.df[,1:(i+1)],
                 kernel = "radial", cost = 1, gamma = 0.1)
  svm.prediction <- predict(svm.fit, newdata = income.pca.df)
  income.pca.accuracy.svm[i] = mean(svm.prediction == income.label)
}

income.psas.accuracy.svm = rep(NA, L)
for(i in 1:L){
  svm.fit <- svm(label ~ ., data = income.psas.df[,1:(i+1)],
                 kernel = "radial", cost = 1, gamma = 0.1)
  svm.prediction <- predict(svm.fit, newdata = income.psas.df)
  income.psas.accuracy.svm[i] = mean(svm.prediction == income.label)
}

income.psao.accuracy.svm = rep(NA, L)
for(i in 1:L){
  svm.fit <- svm(label ~ ., data = income.psao.df[,1:(i+1)],
                 kernel = "radial", cost = 1, gamma = 0.1)
  svm.prediction <- predict(svm.fit, newdata = income.psao.df)
  income.psao.accuracy.svm[i] = mean(svm.prediction == income.label)
}

income.apca.accuracy.svm = rep(NA, L)
for(i in 1:L){
  svm.fit <- svm(label ~ ., data = income.apca.df[,1:(i+1)],
                 kernel = "radial", cost = 1, gamma = 0.1)
  svm.prediction <- predict(svm.fit, newdata = income.apca.df)
  income.apca.accuracy.svm[i] = mean(svm.prediction == income.label)
}

accuracy.svm = as.data.frame(rbind(income.pca.accuracy.svm,
                                   income.psas.accuracy.svm,
                                   income.psao.accuracy.svm,
                                   income.apca.accuracy.svm)) %>%
  mutate(Method = factor(c('PCA','PSAS','PSAO','Log-ratio'), levels = c('PCA','Log-ratio','PSAS','PSAO')))

accuracy.svm %>%
  melt(id.vars = 'Method') %>%
  mutate(variable = as.numeric(as.factor(variable))) %>%
  ggplot(aes(x=variable, y=value, group = Method, color = Method)) +
  theme_bw() +
  geom_line() +
  geom_point() +
  ggtitle('income, SVM') +
  scale_x_continuous(name = 'Rank', breaks = 1:L) +
  scale_y_continuous(name = 'Accuracy', limits = c(0,1))
ggsave('revision round 2/cf_income_svm.jpeg', width = 8, height = 4)

