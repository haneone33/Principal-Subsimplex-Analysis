## orthogonal projection onto the hyperplane spanned by the simplex
proj_simplex <- function(X){
  if(!inherits(X, 'matrix')){
    X = as.matrix(as.data.frame(X))
  }
  d = ncol(X) - 1
  nvec = t(t(rep(1, d+1)/sqrt(d+1)))
  return(X - (X %*% nvec - 1/sqrt(d+1)) %*% t(nvec))
}

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


ternary_pc <- function(X, X.approx, X.colors, X.pc1, X.center){
  
  g.pc <-  empty_tern(colnames(X)[1:3]) +
    theme(plot.title = element_text(hjust = 0.5, size = 15)) +
    ## PC1 and center
    geom_path(data = to_2d(X.pc1), color = 'red', linewidth = 1) +
    geom_point(data = to_2d(X.center), color = 'black', cex = 3) +
    ## original data points
    geom_point(data = to_2d(X), shape = 16, color = X.colors, cex = 2) +
    ## 1 dimensional approximations
    geom_point(data = to_2d(X.approx), shape = 1, color = X.colors, cex = 2)
  
  return(g.pc)
  
}


compare_ternary <- function(X.res, X.label, main = NULL){
  
  m = 101
  t.grid.pca = seq(-5, 5, length.out = 101)
  t.grid.apca = c(-exp(seq(5,-5,length.out = 51)), exp(seq(-5,5,length.out = 51)))
  t.grid.psa = seq(0, 1, length.out = 101)
  X.palette = gg_color_hue(length(unique(X.label)))
  X.colors = X.palette[X.label]
  
  ###############################
  ## data
  ###############################
  
  g.data <- empty_tern(colnames(X.res$X)[1:3]) +
    ggtitle('Data') +
    theme(plot.title = element_text(hjust = 0.5, size = 15)) +
    geom_point(data = to_2d(X.res$X), shape = 16, color = X.colors, cex = 2)
  
  ################################
  ## PCA
  ################################
  
  pca.pc1 = sweep(t.grid.pca %*% t(X.res$pca$loadings[,1]), 2, X.res$pca$center, '+')
  pca.pc1 = as.data.frame(pca.pc1)
  pca.center = as.data.frame(t(X.res$pca$center))
  
  g.pca <- ternary_pc(X.res$X, X.res$pca$pts.approx$`r=2`, X.colors, pca.pc1, pca.center) +
    ggtitle('PCA')
  
  for(i in 1:nrow(X.res$X)){
    proj.path = rbind(X.res$X[i,], X.res$pca$pts.approx$`r=2`[i,])
    g.pca <- g.pca + geom_path(data = to_2d(proj.path), color = X.colors[i])
  }  
  
  ####################################
  ## APCA
  ####################################
  
  apca.pc1 = sweep(t.grid.apca %*% t(X.res$apca$loadings[,1]), 2, X.res$apca$center, '+')
  apca.pc1 = as.data.frame(clrInv(apca.pc1))
  apca.center = as.data.frame(t(clrInv(X.res$apca$center)))
  
  g.apca <- ternary_pc(X.res$X, X.res$apca$pts.approx$`r=2`, X.colors, apca.pc1, apca.center) +
    ggtitle('Log-ratio PCA')
  
  for(i in 1:nrow(X.res$X)){
    x = as.numeric(clr(X.res$apca$pts.approx[[length(X.res$apca$pts.approx)]][i,]))
    y = as.numeric(clr(X.res$apca$pts.approx$`r=2`[i,]))
    proj.path = cbind(seq(x[1], y[1], length.out = m),
                      seq(x[2], y[2], length.out = m),
                      seq(x[3], y[3], length.out = m))
    proj.path = as.data.frame(clrInv(proj.path))
    g.apca <- g.apca + geom_path(data = to_2d(proj.path), color = X.colors[i])
  }  
  
  ################################
  ## Power PCA
  ################################
  alpha = 1/2
  
  power_pca.pc1 = sweep(t.grid.pca %*% t(X.res$power_pca$loadings[,1]), 2, X.res$power_pca$center, '+')
  power_pca.pc1 = as.data.frame(power_pca.pc1)
  power_pca.center = as.data.frame(t(X.res$power_pca$center))
  
  g.power_pca <- ternary_pc(X.res$X, proj_simplex(X.res$power_pca$pts.approx$`r=2`), X.colors,
                            proj_simplex(power_pca.pc1), proj_simplex(power_pca.center)) +
    ggtitle('Power Transform PCA')
  
  for(i in 1:nrow(X.res$X)){
    x = as.numeric(X.res$X[i,]^alpha) # power transform
    y = as.numeric(X.res$power_pca$pts.approx$`r=2`[i,]) # approximation
    z = as.numeric(X.res$X[i,]) # original point
    proj.path.xy = as.data.frame(rbind(x,y))
    proj.path.xz = as.data.frame(rbind(x,z))
    g.power_pca <- g.power_pca + 
      geom_point(data = to_2d(proj_simplex(t(x))), color = X.colors[i], pch = "*", size = 5) +
      geom_path(data = to_2d(proj_simplex(proj.path.xz)), color = X.colors[i]) +
      geom_path(data = to_2d(proj_simplex(proj.path.xy)), color = X.colors[i])
  }  
  
  t.grid = seq(0,pi/2,length.out = 51)
  boundary.12 = cbind(cos(t.grid),sin(t.grid),0)
  boundary.13 = cbind(cos(t.grid),0,sin(t.grid))
  boundary.23 = cbind(0,cos(t.grid),sin(t.grid))
  g.power_pca <- g.power_pca +
    geom_path(data = to_2d(proj_simplex(boundary.12)), color = 'black') +
    geom_path(data = to_2d(proj_simplex(boundary.23)), color = 'black') +
    geom_path(data = to_2d(proj_simplex(boundary.13)), color = 'black')
  
  ###############################################
  ## PSA-S
  ###############################################
  
  psas.pc1 = cbind(t.grid.psa, 1-t.grid.psa) %*% X.res$psas$vertices$`r=2`
  psas.pc1 = as.data.frame(psas.pc1)
  psas.center = as.data.frame(X.res$psas$vertices$`r=1`)
  
  g.psas <-  ternary_pc(X.res$X, X.res$psas$pts.approx$`r=2`, X.colors, psas.pc1, psas.center) +
    ggtitle('PSA-S')
  
  for(i in 1:nrow(X.res$X)){
    x = as.numeric(X.res$X[i,])
    y = as.numeric(X.res$psas$pts.approx$`r=2`[i,])
    proj.path = cbind(seq(x[1], y[1], length.out = m),
                      seq(x[2], y[2], length.out = m),
                      seq(x[3], y[3], length.out = m))
    g.psas <- g.psas + geom_path(data = to_2d(proj.path), color = X.colors[i])
  }  
  
  ###############################################
  ## PSA-O
  ###############################################
  
  psao.pc1 = cbind(t.grid.psa, 1-t.grid.psa) %*% X.res$psao$vertices$`r=2`
  psao.pc1 = as.data.frame(psao.pc1)
  psao.center = as.data.frame(X.res$psao$vertices$`r=1`)
  
  g.psao <-  ternary_pc(X.res$X, X.res$psao$pts.approx$`r=2`, X.colors, psao.pc1, psao.center) +
    ggtitle('PSA-O')
  
  for(i in 1:nrow(X.res$X)){
    x = as.numeric(X.res$X[i,])
    y = as.numeric(X.res$psao$pts.approx$`r=2`[i,])
    proj.path = cbind(seq(x[1], y[1], length.out = m),
                      seq(x[2], y[2], length.out = m),
                      seq(x[3], y[3], length.out = m))
    g.psao <- g.psao + geom_path(data = to_2d(proj.path), color = X.colors[i])
  }  
  
  #################
  ## binding
  #################
  
  g = plot_grid(g.data, g.psas, g.psao, g.pca, g.power_pca, g.apca, nrow = 2)
  if(!is.null(main)){
    g.title = ggdraw() + draw_label(main, size = 20)
    g <- plot_grid(g.title, g, ncol=1, rel_heights=c(0.1, 1))
  }
  
  g
}
