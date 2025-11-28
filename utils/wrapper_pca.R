#' @title Principal Component Analysis for PSA Comparison
#' @description A wrapper function of `princomp()` for comparison to PSA. Manipulates
#' result of `princomp()` in a similar format to `psa()`.
#'
#' @param X a data matrix.
#' @return A list
#' + Xhat a list of lower dimensional representations with respect to the original basis
#' + scores a matrix of scores.
#' + rss a vector of residual sums of squares.
#' + modes a list of modes of variation. The \eqn{r}th element of the list is the difference of
#' rank \eqn{r} approximations to rank \eqn{r-1} approximations.
#' + loadings a matrix of loading vectors.
#' + center mean of the data
#' @seealso [psa()]
#'
#' @export
comp_pca <- function(X){
  n = nrow(X)
  d = ncol(X)-1
  var.names = colnames(X)

  res = list(Xhat = list(), scores = NULL, RSS = 0, loadings = NULL, center = NULL)
  X.pca = stats::princomp(X)

  n = nrow(X)
  d = ncol(X)-1

  res$Xhat[['r=0']] = matrix(rep(X.pca$center, n), byrow=T, nrow = n, ncol = d+1)
  for(r in 1:d){
    res$Xhat[[paste0('r=',r)]] = res$Xhat[[paste0('r=',r-1)]] +
      as.matrix(X.pca$scores[,r]) %*% t(X.pca$loadings[,r])
  }
  res$scores = X.pca$scores
  colnames(res$scores) = paste0('PC',1:(d+1))
  res$RSS = colSums(res$scores^2)

  # res$modes[[paste0('r=',d+1)]] = matrix(NA, n, d+1)
  # colnames(res$modes[[paste0('r=',d+1)]]) = var.names
  # for(i in d:1){
  #   res$modes[[paste0('r=',i)]] = res$Xhat[[paste0('r=',i+1)]] - res$Xhat[[paste0('r=',i)]]
  # }
  # res$modes = res$modes[(d+1):1]

  res$loadings = X.pca$loadings
  class(res$loadings) <- 'matrix'
  colnames(res$loadings) = paste0('PC',1:(d+1))
  res$center = X.pca$center
  names(res$center) = colnames(X)

  return(res)
}

#' @title Power Transform Principal Component Analysis for PSA Comparison
#' @description A wrapper function of `princomp()` applied to power transformed data.
#'
#' @param X a data matrix.
#' @param alpha a scalar \eqn{alpha>0} for power transform.
#' @inherit comp_pca return
#' @seealso [comp_pca()], [psa()]
#'
#' @export
comp_power_pca <- function(X, alpha = 1/2){
  X = X ** alpha
  res = comp_pca(X)
  return(res)
}

#' @title Replace zeros by a fraction of columnwise minimum of nonzero observations
#'
#' @param X a data matrix. Columns are features and rows are samples.
#' @param alpha a positive scalar between 0 and 1. Zeros are replaced by alpha times
#'  the minimum of nonzero observations
#'
#' @export
replace_zero <- function(X, alpha = 0.5){
  for(j in 1:ncol(X)){
    x = X[,j]
    x.min = min(x[x>0])
    X[x<=0,j] = alpha*x.min
  }
  X = to_simplex(X)
  return(X)
}

#' @title Log-Ratio Principal Component Analysis for PSA Comparison
#' @description A wrapper function of `princomp.acomp()` for comparison to PSA.
#' Zeros are substituted by half of the overall nonzero minimum. Manipulates
#' result of `princomp.acomp()` in a similar format to `psa()`.
#'
#' @param X a data matrix.
#' @inherit comp_pca return
#' @seealso [comp_pca()], [psa()]
#'
#' @export
comp_apca <- function(X){

  if(min(X) < 0){
    stop('X contains negative values')
  }
  if(min(X) <= 0){
    warning('X contains zero values')
    X = replace_zero(X)
  }

  n = nrow(X)
  d = ncol(X)-1
  var.names = colnames(X)

  res = list(Xhat = list(), scores = NULL, RSS = 0,
             loadings = NULL, center = NULL)
  X.apca = stats::princomp(compositions::acomp(X))

  ## representations in clr coordinate
  res$Xhat[['r=0']] = matrix(rep(X.apca$center, n), byrow=T, nrow = n, ncol = d+1)
  colnames(res$Xhat[['r=0']]) = var.names
  for(r in 1:d){
    res$Xhat[[paste0('r=',r)]] = res$Xhat[[paste0('r=',r-1)]] +
      as.matrix(X.apca$scores[,r]) %*% t(X.apca$loadings[,r])
  }
  res$scores = X.apca$scores
  colnames(res$scores) = paste0('PC',1:(d+1))
  res$RSS = colSums(res$scores^2)

  ## representations in simplex coordinate
  for(r in 0:d){
    m = compositions::clrInv(res$Xhat[[paste0('r=',r)]])
    class(m) = 'numeric'
    res$Xhat[[paste0('r=',r)]] = m
  }

  # res$modes[[paste0('r=',d+1)]] = matrix(NA, n, d+1)
  # colnames(res$modes[[paste0('r=',d+1)]]) = var.names
  # for(i in d:1){
  #   res$modes[[paste0('r=',i)]] = res$Xhat[[paste0('r=',i+1)]] - res$Xhat[[paste0('r=',i)]]
  # }
  # res$modes = res$modes[(d+1):1]

  res$loadings = X.apca$loadings
  class(res$loadings) <- 'matrix'
  colnames(res$loadings) = paste0('PC',1:d)
  res$center = X.apca$center
  names(res$center) = colnames(X)

  return(res)
}
