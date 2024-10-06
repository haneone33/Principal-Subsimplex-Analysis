#' @title Principal Nested Simplices Analysis
#' @description Estimate PSA-S or PSA-O of given data matrix.
#' 
#' @param type `s` for PSA-S or `o` for PSA-O.
#' @param X A data matrix.
#' @param testweights A vector of weights to try.
#' @return A list
#' + vertices a list of matrix representing vertices of the lower dimensional subsimplex. 
#'The $r$th element of the list corresponds to the rank $r-1$ subsimplex.
#' + pts a list of lower dimensional representations with respect to the reduced basis `vertices`
#' + pts.approx a list of lower dimensional representations with respect to the original basis
#' + scores a matrix of scores.
#' + rss a vector of residual sums of squares.
#' + modes a list of modes of variation. The $r$th element of the list is the difference of 
#' rank $r$ approximations to rank $r-1$ approximations.
#' + loadings a matrix of loading vectors.
#' + const.info a data frame of merged vertices and merging weight at each merge.
#' + dendrogram.input additional information to apply `plotdendrogram2()`.
#' @seealso [plotdendrogram2()]
#' 
#' @export

psa <- function(type, X, testweights = seq(0, 1, length.out = 100)){
  n = nrow(X)
  d = ncol(X)-1
  var.names = colnames(X)
  
  res = list(vertices = list(),       # vertices of low rank subsimplices
             pts = list(),            # low rank approximations w.r.t. vertices
             pts.approx = list(),     # low rank approximations w.r.t. original vertices
             scores = list(),         # scalars
             rss = list(),            # residual sum of squares
             modes = list(),          # residuals (vectors; in Delta_d scale)
             loadings = NULL,         # loading vectors (matrix)
             const.info = list(),     # merge indices and weight
             dendrogram.input = NULL) # raw output from psa-s and psa-o
  
  if(type == 's'){
    
    out = pssa(X, testweights = testweights)
    res$vertices = out$vertices
    res$pts = out$pts
    for(i in names(out$pts)){
      res$pts.approx[[i]] = out$pts[[i]] %*% out$vertices[[i]]
    }
    res$scores = out$scores
    res$rss = sapply(res$scores, function(vec) sum(vec^2))
    res$modes = out$modes
    res$const.info = out$merges
    res$dendrogram.input = out
    
  }else if(type == 'o'){
    
    out = psoa(X, testweights = testweights)
    res$vertices = out$vertices
    res$pts = out$pts
    for(i in names(out$pts)){
      res$pts.approx[[i]] = divL1(divL2(out$pts[[i]]) %*% out$vertices[[i]])
    }
    res$scores = out$scores
    res$rss = sapply(res$scores, function(vec) sum(vec^2))
    res$modes = out$modes
    res$const.info = out$merges
    res$dendrogram.input = out
    
  }else{
    stop('type should be either `s` or `o`')
  }
  
  res$vertices = res$vertices[(d+1):1]
  res$pts = res$pts[(d+1):1]
  res$pts.approx = res$pts.approx[(d+1):1]
  res$scores = res$scores[(d+1):1]
  res$rss = res$rss[(d+1):1]
  res$modes = res$modes[(d+1):1]
  
  loadings.ls = get_loading(res$dendrogram.input)
  res$loadings = do.call('cbind', lapply(loadings.ls, function(l) l$'v2'-l$'v1'))
  res$loadings = res$loadings[,ncol(res$loadings):1]
  
  res$scores = do.call('cbind', res$scores)
  colnames(res$scores) = paste0('Comp.',1:ncol(res$scores))
  colnames(res$loadings) = paste0('Comp.',1:ncol(res$loadings))
  res$rss[d+1]=0
  return(res)
}
