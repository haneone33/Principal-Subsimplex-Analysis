#' @title Map data to the simplex
#' @description Map data to the simplex while preserving proportions.
#' @param X a data matrix
#' @return a matrix of the same dimension as `X`
#' @export
to_simplex <- function(X){
  if(min(X)<0){
    # warning('X contains negative values')
    X = pmax(X, 0)
  }
  if(is.null(dim(X))){
    X = X/sum(X)
  }else{
    X = t(apply(X, 1, function(x){
      x = x/sum(x)
    }))
  }
  return(X)
}
