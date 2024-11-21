## orthogonal projection onto the hyperplane spanned by the simplex
## for displaying power transform PCA
proj_simplex <- function(X){
  if(!inherits(X, 'matrix')){
    X = as.matrix(as.data.frame(X))
  }
  d = ncol(X) - 1
  nvec = t(t(rep(1, d+1)/sqrt(d+1)))
  return(X - (X %*% nvec - 1/sqrt(d+1)) %*% t(nvec))
}


draw_ternary_pc <- function(X, X.approx, X.pc1, X.center, X.colors){
  ## X, X.approx, X.pc1, X.center are data frames with three columns
  g.pc <-  empty_tern(colnames(X)[1:3]) +
    theme(plot.title = element_text(hjust = 0.5, size = 20)) +
    geom_path(data = to_2d(X.pc1), color = 'red', linewidth = 1) +            # PC1
    geom_point(data = to_2d(X.center), color = 'black', size = 3) +           # center
    geom_point(data = to_2d(X), shape = 16, color = X.colors, size = 2) +     # original data points
    geom_point(data = to_2d(X.approx), shape = 1, color = X.colors, size = 2) # 1 dimensional approximations
  
  return(g.pc)
}



ternary_pc <- function(X.res, X.label, type = 'data', t.grid = NULL){
  ## type is one of 'data', 'psas', 'psao', 'pca', 'power', 'logratio'

  X.palette = gg_color_hue(length(unique(X.label)))
  X.colors = X.palette[X.label]
  m = 30 # resolution of paths connecting points to their approximations
    
  if(type == 'data'){

    g.data <- empty_tern(colnames(X.res$X)[1:3]) +
      ggtitle('Data') +
      theme(plot.title = element_text(hjust = 0.5, size = 20)) +
      geom_point(data = to_2d(X.res$X), shape = 16, color = X.colors, size = 2)
    
    return(g.data)
    
  }else if(type == 'pca'){

    if(is.null(t.grid)){
      t.grid = seq(-5, 5, length.out = 101)
    }
    
    X.pc1 = sweep(t.grid %*% t(X.res$pca$loadings[,1]), 2, X.res$pca$center, '+')
    X.pc1 = as.data.frame(X.pc1)
    X.center = as.data.frame(t(X.res$pca$center))
    
    g <- draw_ternary_pc(X.res$X, X.res$pca$pts.approx$`r=2`, X.pc1, X.center, X.colors) +
      ggtitle('PCA')
    
    for(i in 1:nrow(X.res$X)){
      proj.path = rbind(X.res$X[i,], X.res$pca$pts.approx$`r=2`[i,])
      g <- g + geom_path(data = to_2d(proj.path), color = X.colors[i])
    }
    
    return(g)
    
  }else if(type == 'logratio'){

    if(is.null(t.grid)){
      t.grid = c(-exp(seq(5,-5,length.out = 51)), exp(seq(-5,5,length.out = 51)))
    }
    
    X.pc1 = sweep(t.grid %*% t(X.res$apca$loadings[,1]), 2, X.res$apca$center, '+')
    X.pc1 = as.data.frame(clrInv(X.pc1))
    X.center = as.data.frame(t(clrInv(X.res$apca$center)))
    
    g <- draw_ternary_pc(X.res$X, X.res$apca$pts.approx$`r=2`, X.pc1, X.center, X.colors) +
      ggtitle('Log-ratio PCA')
    
    for(i in 1:nrow(X.res$X)){
      x = as.numeric(clr(X.res$apca$pts.approx[[length(X.res$apca$pts.approx)]][i,]))
      y = as.numeric(clr(X.res$apca$pts.approx$`r=2`[i,]))
      proj.path = cbind(seq(x[1], y[1], length.out = m),
                        seq(x[2], y[2], length.out = m),
                        seq(x[3], y[3], length.out = m))
      proj.path = as.data.frame(clrInv(proj.path))
      g <- g + geom_path(data = to_2d(proj.path), color = X.colors[i])
    } 
    
    return(g)
    
  }else if(type == 'power'){

    alpha = 1/2
    if(is.null(t.grid)){
      t.grid = seq(-5, 5, length.out = 101)
    }
    
    X.pc1 = sweep(t.grid %*% t(X.res$power_pca$loadings[,1]), 2, X.res$power_pca$center, '+')
    X.pc1 = as.data.frame(X.pc1)
    X.center = as.data.frame(t(X.res$power_pca$center))
    
    g <- draw_ternary_pc(X.res$X, proj_simplex(X.res$power_pca$pts.approx$`r=2`),
                              proj_simplex(X.pc1), proj_simplex(X.center), X.colors) +
      ggtitle('Power Transform PCA')
    
    for(i in 1:nrow(X.res$X)){
      x = as.numeric(X.res$X[i,]^alpha) # power transform
      y = as.numeric(X.res$power_pca$pts.approx$`r=2`[i,]) # approximation
      z = as.numeric(X.res$X[i,]) # original point
      proj.path.xy = as.data.frame(rbind(x,y))
      proj.path.xz = as.data.frame(rbind(x,z))
      g <- g + 
        geom_point(data = to_2d(proj_simplex(t(x))), color = X.colors[i], pch = "*", size = 5) +
        geom_path(data = to_2d(proj_simplex(proj.path.xz)), color = X.colors[i]) +
        geom_path(data = to_2d(proj_simplex(proj.path.xy)), color = X.colors[i])
    }  
    
    ## draw boundary of orthant
    tt.grid = seq(0,pi/2,length.out = 51)
    boundary.12 = cbind(cos(tt.grid),sin(tt.grid),0)
    boundary.13 = cbind(cos(tt.grid),0,sin(tt.grid))
    boundary.23 = cbind(0,cos(tt.grid),sin(tt.grid))
    g <- g +
      geom_path(data = to_2d(proj_simplex(boundary.12)), color = 'black') +
      geom_path(data = to_2d(proj_simplex(boundary.23)), color = 'black') +
      geom_path(data = to_2d(proj_simplex(boundary.13)), color = 'black')
    
    return(g)
    
  }else if(type == 'psas'){
    
    if(is.null(t.grid)){
      t.grid = seq(0, 1, length.out = 101)
    }
    
    X.pc1 = cbind(t.grid, 1-t.grid) %*% X.res$psas$vertices$`r=2`
    X.pc1 = as.data.frame(X.pc1)
    X.center = as.data.frame(X.res$psas$vertices$`r=1`)
    
    g <- draw_ternary_pc(X.res$X, X.res$psas$pts.approx$`r=2`, X.pc1, X.center, X.colors) +
      ggtitle('PSA-S')
    
    for(i in 1:nrow(X.res$X)){
      x = as.numeric(X.res$X[i,])
      y = as.numeric(X.res$psas$pts.approx$`r=2`[i,])
      proj.path = cbind(seq(x[1], y[1], length.out = m),
                        seq(x[2], y[2], length.out = m),
                        seq(x[3], y[3], length.out = m))
      g <- g + geom_path(data = to_2d(proj.path), color = X.colors[i])
    }  
    
    return(g)
    
  }else if(type == 'psao'){
    
    if(is.null(t.grid)){
      t.grid = seq(0, 1, length.out = 101)
    }
    
    X.pc1 = cbind(t.grid, 1-t.grid) %*% X.res$psao$vertices$`r=2`
    X.pc1 = as.data.frame(X.pc1)
    X.center = as.data.frame(X.res$psao$vertices$`r=1`)
    
    g <- draw_ternary_pc(X.res$X, X.res$psao$pts.approx$`r=2`, X.pc1, X.center, X.colors) +
      ggtitle('PSA-O')
    
    for(i in 1:nrow(X.res$X)){
      x = as.numeric(X.res$X[i,])
      y = as.numeric(X.res$psao$pts.approx$`r=2`[i,])
      proj.path = cbind(seq(x[1], y[1], length.out = m),
                        seq(x[2], y[2], length.out = m),
                        seq(x[3], y[3], length.out = m))
      g <- g + geom_path(data = to_2d(proj.path), color = X.colors[i])
    }  
    
    return(g)
    
  }else{
    stop('Incorrect type')
  }
  
}
