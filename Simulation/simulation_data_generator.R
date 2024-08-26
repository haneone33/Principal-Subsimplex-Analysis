library(compositions)
library(ggplot2)

#' @title Multivariate normal random sample
#' 
#' @param n number of points to generate.
#' @param center mean.
#' @param Sigma variance matrix. Alternatively, a scalar can be given.
rmnorm <- function(n, center, Sigma){
  if(is.null(dim(Sigma))){
    Sigma = diag(Sigma, length(center))
  }
  p = length(center)
  X = matrix(rnorm(n*p, mean = 0), p, n)
  X = t(chol(Sigma)) %*% X
  X = sweep(X, 1, center, '+')
  return(t(X))
}

#' @title Create a single cluster
#' @description Create a normally-distributed sample on a simplex.
#' @param n number of points to generate.
#' @param center center of the cluster.
#' @param Sigma variance matrix. Alternatively, a scalar can be given
make_cluster <- function(center, n, Sigma, margin = 0){
  X = rmnorm(n, center, Sigma)
  X = to_simplex(X)
  return(X)
}

#' @title Make a toy example
#' 
#' @param centers list of cluster centers.
#' @param n size of each cluster.
#' @param sigma2 vector of variances of clusters.
#' 
#' @return a matrix of sample points.
#' @export
make_toy_example <- function(centers, n = 40, sigma2 = 0.04^2, margin = 1e-8){
  if(is.null(dim(sigma2))){
    sigma2 = rep(sigma2, length(centers))
  }
  k = length(centers)
  p = length(centers[[1]])
  X = matrix(NA,0,p)
  for(i in 1:k){
    X = rbind(X, make_cluster(centers[[i]], n, sigma2[i], margin))
  }
  return(X)
}


