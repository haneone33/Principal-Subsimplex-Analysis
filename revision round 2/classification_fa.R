library(devtools)
load_all('../psacomp/')
library(MASS)
library(e1071)
library(compositions)
library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(GGally)

fa = read.csv('revision round 2/fatty acid.csv')
fa.label = as.factor(fa$Species)
fa.X = to_simplex(as.matrix(fa[,-1]))

fa.pca = comp_pca(fa.X)
fa.pca.df = data.frame(label = fa.label) %>%
  cbind(as.data.frame(fa.pca$scores))
ggplot(fa.pca.df, aes(x=PC1, y=PC2, col = fa.label)) +
  theme_bw() +
  geom_point()

fa.apca = comp_apca(fa.X)
fa.apca.df = data.frame(label = fa.label) %>%
  cbind(as.data.frame(fa.apca$scores))
ggplot(fa.apca.df, aes(x=PC1, y=PC2, col = fa.label)) +
  theme_bw() +
  geom_point()

fa.psas = psa('s', fa.X)
fa.psas = readRDS('revision round 2/fa_psas.RDS')
fa.psas.df = data.frame(label = fa.label) %>%
  cbind(as.data.frame(fa.psas$scores))
ggplot(fa.psas.df, aes(x=PC1, y=PC2, col = fa.label)) +
  theme_bw() +
  geom_point()

fa.psao = psa('o', fa.X)
fa.psao = readRDS('revision round 2/fa_psao.RDS')
fa.psao.df = data.frame(label = fa.label) %>%
  cbind(as.data.frame(fa.psao$scores))
ggplot(fa.psao.df, aes(x=PC1, y=PC2, col = fa.label)) +
  theme_bw() +
  geom_point()

L = 10
fa.pca.accuracy = rep(NA, L)
for(i in 1:L){
  qda.fit <- qda(label ~ ., data = fa.pca.df[,1:(i+1)])
  qda.prediction <- predict(qda.fit, newdata = fa.pca.df)
  fa.pca.accuracy[i] = mean(qda.prediction$class == fa.label)
}

fa.psas.accuracy = rep(NA, L)
for(i in 1:L){
  qda.fit <- qda(label ~ ., data = fa.psas.df[,1:(i+1)])
  qda.prediction <- predict(qda.fit, newdata = fa.psas.df)
  fa.psas.accuracy[i] = mean(qda.prediction$class == fa.label)
}

fa.psao.accuracy = rep(NA, L)
for(i in 1:L){
  qda.fit <- qda(label ~ ., data = fa.psao.df[,1:(i+1)])
  qda.prediction <- predict(qda.fit, newdata = fa.psao.df)
  fa.psao.accuracy[i] = mean(qda.prediction$class == fa.label)
}

fa.apca.accuracy = rep(NA, L)
for(i in 1:L){
  qda.fit <- qda(label ~ ., data = fa.apca.df[,1:(i+1)])
  qda.prediction <- predict(qda.fit, newdata = fa.apca.df)
  fa.apca.accuracy[i] = mean(qda.prediction$class == fa.label)
}

accuracy = as.data.frame(rbind(fa.pca.accuracy,
                               fa.psas.accuracy,
                               fa.psao.accuracy,
                               fa.apca.accuracy)) %>%
  mutate(Method = factor(c('PCA','PSAS','PSAO','Log-ratio'), levels = c('PCA','Log-ratio','PSAS','PSAO')))

accuracy %>%
  melt(id.vars = 'Method') %>%
  mutate(variable = as.numeric(as.factor(variable))) %>%
  ggplot(aes(x=variable, y=value, group = Method, color = Method)) +
  theme_bw() +
  geom_line() +
  geom_point() +
  ggtitle('Fatty Acid, QDA') +
  scale_x_continuous(name = 'Rank', breaks = 1:L) +
  scale_y_continuous(name = 'Accuracy', limits = c(0,1))
ggsave('revision round 2/cf_fatty_acid_qda.jpeg', width = 8, height = 4)




L = 10
fa.pca.accuracy.svm = rep(NA, L)
for(i in 1:L){
  svm.fit <- svm(label ~ ., data = fa.pca.df[,1:(i+1)],
                 kernel = "radial", cost = 1)
  svm.prediction <- predict(svm.fit, newdata = fa.pca.df)
  fa.pca.accuracy.svm[i] = mean(svm.prediction == fa.label)
}

hc.psas.accuracy.svm = rep(NA, L)
for(i in 1:L){
  df = as.data.frame(hc.psas$Xhat_reduced[[i+1]][,1:i]) %>%
    mutate(label = hc.label)
  svm.fit <- svm(label ~ ., data = df,
                 kernel = "radial", cost = 1)
  svm.prediction <- predict(svm.fit, newdata = df)

  # svm.fit <- svm(label ~ ., data = hc.psas.df[,1:(i+1)],
  #                kernel = "radial", cost = 1)
  # svm.prediction <- predict(svm.fit, newdata = hc.psas.df)

  hc.psas.accuracy.svm[i] = mean(svm.prediction == hc.label)
}

hc.psao.accuracy.svm = rep(NA, L)
for(i in 1:L){
  df = as.data.frame(hc.psao$Xhat_reduced[[i+1]][,1:i]) %>%
    mutate(label = hc.label)
  svm.fit <- svm(label ~ ., data = df,
                 kernel = "radial", cost = 1)
  svm.prediction <- predict(svm.fit, newdata = df)

  # svm.fit <- svm(label ~ ., data = hc.psao.df[,1:(i+1)],
  #                kernel = "radial", cost = 1)
  # svm.prediction <- predict(svm.fit, newdata = hc.psao.df)

  hc.psao.accuracy.svm[i] = mean(svm.prediction == hc.label)
}

fa.apca.accuracy.svm = rep(NA, L)
for(i in 1:L){
  svm.fit <- svm(label ~ ., data = fa.apca.df[,1:(i+1)],
                 kernel = "radial", cost = 1)
  svm.prediction <- predict(svm.fit, newdata = fa.apca.df)
  fa.apca.accuracy.svm[i] = mean(svm.prediction == fa.label)
}

accuracy.svm = as.data.frame(rbind(fa.pca.accuracy.svm,
                                   fa.psas.accuracy.svm,
                                   fa.psao.accuracy.svm,
                                   fa.apca.accuracy.svm)) %>%
  mutate(Method = factor(c('PCA','PSAS','PSAO','Log-ratio'), levels = c('PCA','Log-ratio','PSAS','PSAO')))

accuracy.svm %>%
  melt(id.vars = 'Method') %>%
  mutate(variable = as.numeric(as.factor(variable))) %>%
  ggplot(aes(x=variable, y=value, group = Method, color = Method)) +
  theme_bw() +
  geom_line() +
  geom_point() +
  ggtitle('Fatty Acid, SVM') +
  scale_x_continuous(name = 'Rank', breaks = 1:L) +
  scale_y_continuous(name = 'Accuracy', limits = c(0,1))
ggsave('revision round 2/cf_fatty_acid_svm.jpeg', width = 8, height = 4)
