#' @title Comparison of PSA and Benchmark Methods
#' @description A wrapper function of PSA and benchmark methods for convenience
#'
#' @param X a data matrix
#' @return a list of data matrix `X` and outcomes of PSA-S, PSA-O, PCA, log-ratio PCA,
#' and power transform PCA with power 1/2.
#'
#' @seealso [psa()], [comp_pca()], [comp_apca()], [comp_power_pca()]
#' @export

compare_analysis <- function(X){
  X = as.matrix(X)
  if(is.null(colnames(X))){
    colnames(X) = paste0('V',1:ncol(X))
  }

  X.psas = tryCatch(psa('s', X), error = function(e) e)
  X.psao = tryCatch(psa('o', X), error = function(e) e)
  X.pca = tryCatch(comp_pca(X), error = function(e) e)
  X.apca = tryCatch(comp_apca(X), error = function(e) e)
  X.power_pca = tryCatch(comp_power_pca(X), error = function(e) e)

  res = list(X = as.matrix(X),
             psas = X.psas, psao = X.psao,
             pca = X.pca, apca = X.apca, power_pca = X.power_pca)
  return(res)
}

