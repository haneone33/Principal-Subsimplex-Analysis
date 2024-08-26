psa <- function(type, X, ...){
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
             const.info = list(),     # (merge indices and weight) or (chosen index and normal vector)
             dendrogram.input = NULL) # raw output from psa-s and psa-o
  
  if(type == 's'){
    
    out = pssa(X, ...)
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
    
    out = psoa(X, ...)
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
    
  }else if(type == 'g'){
    
    out = gssa(X, ...)
    res$vertices = out$vertices
    res$pts = out$Xhat_rescaled
    res$pts.approx = out$Xhat
    res$scores = out$scores
    res$rss = sapply(res$scores, function(vec) sum(vec^2))
    res$modes[[paste0('r=',d+1)]] = matrix(NA, n, d+1)
    colnames(res$modes[[paste0('r=',d+1)]]) = var.names
    for(i in d:1){
      res$modes[[paste0('r=',i)]] = res$pts.approx[[paste0('r=',i+1)]] - res$pts.approx[[paste0('r=',i)]]
    }
    res$const.info = list(index = out$idx, a = out$a)
    
    
  }else{
    stop('type should be one of `s`, `o`, and `g`.')
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
