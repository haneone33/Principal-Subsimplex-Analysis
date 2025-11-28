## orthogonal projection onto the hyperplane spanned by the simplex for displaying power transform PCA
proj_simplex <- function(X){
  if(!inherits(X, 'matrix')){
    X = as.matrix(as.data.frame(X))
  }
  d = ncol(X) - 1
  nvec = t(t(rep(1, d+1)/sqrt(d+1)))
  return(X - (X %*% nvec - 1/sqrt(d+1)) %*% t(nvec))
}

## rotate 3d vector on the hyperplane spanned by the simplex into the 2d view
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

## create an empty ternary ggplot object
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



ternary_pc <- function(X, X.label, type = 'data', X.res = NULL, t.grid = NULL){
  ## type is one of 'data', 'psas', 'psao', 'pca', 'power', 'logratio'

  X.palette = gg_color_hue(length(unique(X.label)))
  X.colors = X.palette[X.label]
  m = 30 # resolution of paths connecting points to their approximations

  if(type == 'data'){

    g.data <- empty_tern(colnames(X)[1:3]) +
      ggtitle('Data') +
      theme(plot.title = element_text(hjust = 0.5, size = 20)) +
      geom_point(data = to_2d(X), shape = 16, color = X.colors, size = 2)

    return(g.data)

  }else if(type == 'pca'){

    if(is.null(t.grid)){
      t.grid = seq(-5, 5, length.out = 101)
    }

    X.pc1 = sweep(t.grid %*% t(X.res$loadings[,1]), 2, X.res$center, '+')
    X.pc1 = as.data.frame(X.pc1)
    X.center = as.data.frame(t(X.res$center))

    g <- draw_ternary_pc(X, X.res$Xhat$`r=1`, X.pc1, X.center, X.colors) +
      ggtitle('PCA')

    for(i in 1:nrow(X)){
      proj.path = rbind(X[i,], X.res$Xhat$`r=1`[i,])
      g <- g + geom_path(data = to_2d(proj.path), color = X.colors[i])
    }

    return(g)

  }else if(type == 'logratio'){

    if(is.null(t.grid)){
      t.grid = c(-exp(seq(5,-5,length.out = 51)), exp(seq(-5,5,length.out = 51)))
    }

    X.pc1 = sweep(t.grid %*% t(X.res$loadings[,1]), 2, X.res$center, '+')
    X.pc1 = as.data.frame(clrInv(X.pc1))
    X.center = as.data.frame(t(clrInv(X.res$center)))

    g <- draw_ternary_pc(X, X.res$Xhat$`r=1`, X.pc1, X.center, X.colors) +
      ggtitle('Log-ratio PCA')

    for(i in 1:nrow(X)){
      x = as.numeric(clr(X.res$Xhat[[length(X.res$Xhat)]][i,]))
      y = as.numeric(clr(X.res$Xhat$`r=1`[i,]))
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

    X.pc1 = sweep(t.grid %*% t(X.res$loadings[,1]), 2, X.res$center, '+')
    X.pc1 = as.data.frame(X.pc1)
    X.center = as.data.frame(t(X.res$center))

    g <- draw_ternary_pc(X, proj_simplex(X.res$Xhat$`r=1`),
                         proj_simplex(X.pc1), proj_simplex(X.center), X.colors) +
      ggtitle('Power Transform PCA')

    for(i in 1:nrow(X)){
      x = as.numeric(X[i,]**alpha) # power transform
      y = as.numeric(X.res$Xhat$`r=1`[i,]) # approximation
      z = as.numeric(X[i,]) # original point
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

    X.pc1 = cbind(t.grid, 1-t.grid) %*% t(X.res$Vhat$`r=1`)
    X.pc1 = as.data.frame(X.pc1)
    X.center = as.data.frame(t(X.res$Vhat$`r=0`))

    g <- draw_ternary_pc(X, X.res$Xhat$`r=1`, X.pc1, X.center, X.colors) +
      ggtitle('PSA-S')

    for(i in 1:nrow(X)){
      x = as.numeric(X[i,])
      y = as.numeric(X.res$Xhat$`r=1`[i,])
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

    X.pc1 = cbind(t.grid, 1-t.grid) %*% t(X.res$Vhat$`r=1`)
    X.pc1 = as.data.frame(X.pc1)
    X.center = as.data.frame(t(X.res$Vhat$`r=0`))

    g <- draw_ternary_pc(X, X.res$Xhat$`r=1`, X.pc1, X.center, X.colors) +
      ggtitle('PSA-O')

    for(i in 1:nrow(X)){
      x = as.numeric(X[i,])
      y = as.numeric(X.res$Xhat$`r=1`[i,])
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
