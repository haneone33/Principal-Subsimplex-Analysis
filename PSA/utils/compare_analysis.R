compare_analysis <- function(X, X.label = NULL){
  X = as.matrix(X)
  X.label = as.factor(X.label)
  
  X.psas = tryCatch(psa('s', X), error = function(e) e)
  X.psao = tryCatch(psa('o', X), error = function(e) e)
  X.pca = tryCatch(comp_pca(X), error = function(e) e)
  X.apca = tryCatch(comp_apca(X), error = function(e) e)
  X.power_pca = tryCatch(power_pca(X), error = function(e) e)
  
  res = list(X = as.matrix(X), label = X.label,
             psas = X.psas, psao = X.psao,
             pca = X.pca, apca = X.apca, power_pca = X.power_pca)
  return(res)
}
