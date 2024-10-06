#' @title Iteratively replace components with scaled components, and rescale to a composition.
#' @description Replaces the mass in a component with two new components that are multiples of the original component, with the `weight1` the multiple for the first new component and `weight2` the multiple for the second new component
#' The compositional data has the same number of rows as `weight` and begins as all `1`, with one component.
#' @param weight1,weight2 A matrix specifying how the mass in the component is split.
#' @param component A vector of integers specifying the component to be split each time.
pairwisereplace <- function(weight1, weight2, component = rep(1, ncol(weight1))){
  X <- matrix(1, nrow = nrow(weight1))
  for (i in 1:length(component)){
    X <- pairwisereplace_once(X, weight1[, i], weight2[, i], component[i])
  }
  return(X)
}

pairwisereplace_once <- function(X, weight1, weight2, component = 1){
  if(length(weight1)==1){weight1 <- rep(weight1, nrow(X))}
  if(length(weight2)==1){weight2 <- rep(weight2, nrow(X))}
  newcomp1 <- X[, component] * weight1
  newcomp2 <- X[, component] * weight2
  out <- cbind(newcomp1, newcomp2, X[, -component, drop = FALSE])
  out <- out/rowSums(out)
  colnames(out)[1:2] <- paste0(colnames(X)[component], c("a", "b"))
  return(out)
}
