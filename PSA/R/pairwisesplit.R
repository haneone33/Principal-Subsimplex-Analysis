#' @title Iteratively split components by given weight.
#' @description Splits the mass in a component into two, with the `weight` dictating how much of the mass goes to the first new component.
#' The compositional data has the same number of rows as `weight` and begins as all `1`, with one component.
#' @param weight A matrix specifying how the mass in the component is split.
#' @param component A vector of integers specifying the component to be split each time.
#' @noRd
pairwisesubpartition <- function(weight, component = rep(1, ncol(weight))){
  X <- matrix(1, nrow = nrow(weight))
  for (i in 1:length(component)){
    X <- pairwisesubpartition_once(X, weight[, i], component[i])
  }
  return(X)
}

pairwisesubpartition_once <- function(X, weight, component = 1){
  if(length(weight)==1){weight <- rep(weight, nrow(X))}
  newcomp1 <- X[, component] * weight
  newcomp2 <- X[, component] * (1 - weight)
  out <- cbind(newcomp1, newcomp2, X[, -component, drop = FALSE])
  colnames(out)[1:2] <- paste0(colnames(X)[component], c("a", "b"))
  return(out)
}
