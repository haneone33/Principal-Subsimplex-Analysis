to_2d <- function(X){
  if(!inherits(X, 'matrix')){
    X = as.matrix(as.data.frame(X))
  }
  X = X[,1:3, drop=F]/rowSums(X[,1:3, drop=F])
  A = matrix(c(-1,1,0,
               0,0,sqrt(3)),
             byrow=T, 2,3)
  return(as.data.frame(t(A%*%t(X))))
}

empty_tern <- function(axis){
  if(is.null(axis)){
    axis = c('V1','V2','V3')
  }
  ggplot(data = NULL, aes(V1, V2)) +
    theme_void() +
    xlim(-1.2, 1.2) + ylim(-0.25, 1.9) +
    geom_path(data = data.frame(V1 = c(0,-1,1,0), V2 = c(sqrt(3),0,0,sqrt(3))),
              linewidth = 0.7, color = 'black') +
    annotate('text', x = c(-1.1, 1.1, 0), y = c(-0.05, -0.05, sqrt(3)+0.1), label = axis)
}