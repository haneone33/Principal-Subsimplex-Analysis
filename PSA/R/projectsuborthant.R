#' @title Project data to a lower suborthant
#'
#' @param X a data matrix whose coordinates are computed with respect to `V`.
#' @param V matrix for vertices of the reduced simplex.
#' @param v1,v2 vertices to be merged.
#' @param w weight for merge.
#' @details
#' The order of v1 and v2 specify the sign of the scores: points on the v2 side of the orthant have positive scores
#' @return
#'  + `X` projected data points in the new coordinates
#'  + `V` new vertices
#'  + `scores` Signed scores of the projection
#'  + `mergedirection` The direction of the merge, with positive pointing towards positive score points.
#' @noRd

projectsuborthant <- function(X, V = diag(ncol(X)), v1, v2, w){
  #define normal vector to new suborthant, pointing away from v1:
  # normal vector n is in the plane of v1 and v2 (i.e. n = a*v1 + b*v2), and perpendicular to w * v1 + (1-w) * v2
  # t(a*v1 + b*v2) * (w * v1 + (1-w) * v2) = 0
  # ==> a*w + b*(1-w) = 0
  # ==> a = b * (w - 1)/w
  # also a^2 + b^2 = 1 (unit vector)
  # ==> b^2 * {((w-1)/w)^2 + 1} = 1
  # ==> b^2 = 1/{((w-1)/w)^2 + 1} = w^2/{(w-1)^2 + w^2}
  # ==> b = w / sqrt((w-1)^2 + w^2)
  # ==> a = {(w-1)/w} * w / sqrt((w-1)^2 + w^2) = (w-1) / sqrt((w-1)^2 + w^2)
  # above formulas are correct when w = 0 too (a should be -1, and b zero)
  e1 <- e2 <- rep(0, ncol(X)) #specifying basis vectors corresponding to v1 and v2
  e1[v1] <- e2[v2] <- 1
  sizew1mw <- sqrt((w-1)^2 + w^2)
  nvec <- ((w-1)/sizew1mw) * e1 + (w/sizew1mw) * e2

  #project X to sphere
  Xsph <- divL2(X)

  #project Xsph to suborthant
  Xproj <- project.subsphere(Xsph, nvec)

  #Signed scores
  absscores <- acos_clamped(rowSums(Xproj * Xsph))
  sscores <- drop(sign(Xsph %*% nvec) * absscores)

  #return new X to simplex, first converting to new vertices, using current basis as e1, e2 etc
  newVwrtolddiag <- rbind((w/sizew1mw) * e1 + ((1-w)/sizew1mw) * e2,
                diag(ncol(X))[-c(v1, v2), , drop = FALSE])
  newX <- divL1(Xproj %*% t(newVwrtolddiag))
  colnames(newX) <- c(sprintf("(%s,%s)", colnames(X)[v1],  colnames(X)[v2]),
                      colnames(X)[-c(v1, v2)])

  #Record merge direction in terms of input V # the next two bits are compeletely separate from the X above
  mergedirection <- V[v2, ] - V[v1, ]

  #record new vertex in V
  newV <- rbind( ((w/sizew1mw) * V[v1, ] + ((1-w)/sizew1mw) * V[v2, ] ),
                 V[-c(v1, v2), , drop = FALSE]
               )
  # sum(nvec * newV[1, ]) == 0


  return(list(
    X = newX,
    V = newV,
    scores = sscores,
    mergedirection = mergedirection
  ))
}

# For suborthant merging: A function that evaluates each proposed merge and returns a single non-negative value for the quality of the merge. Arguments of `X`, `v1`, `v2`, and `w` at least.
evaluatesuborthant <- function(X, v1, v2, w){
  proj <- projectsuborthant(X = X, V = diag(ncol(X)), v1 = v1, v2 = v2, w = w)
  sum(proj$scores^2)
}

# X is matrix, each row a location
divL2 <- function(X){
  l2sizes <- sqrt(rowSums(X^2))
  return(X / l2sizes)
}

# warning ONLY for non-negative entries
divL1 <- function(X){
  l1sizes <- rowSums(abs(X))
  return(X / l1sizes)
}

#based on project.subsphere in shapes R package, commit 45d4dde. The computation matches equation 4 of Jung et al 2012. Here we know r = pi/2, so sin(r) = 1
# center is the normal vector of the suborthant
# I've added a fix when a point in X is coincident with +/- center
# to arbitrarily put the projected point as the projection of the first coincident point
project.subsphere = function(X, center)
{
   rho <- acos_clamped(X %*% center)
   sinrho1 <- sin(rho - pi/2)
   sinrho2 <- sin(rho)
   b <- sinrho1 %x% t(center)#kronecker product is %x%. Each row should be an element of sin(rho-pi/2) * center
   # c <- sinrho1 %*% t(center) Could also do this - 50% slower, but 25% less memory usage. Makes sense to do this if going to parallel.
   Xproj <- X + b
   Xproj <- Xproj/drop(sinrho2)
   zerothreshold <- 1E-10
   zeropts <- (sin(rho) < zerothreshold)
   if (any(zeropts)){
      ref <- match(FALSE, zeropts)
      Xproj[zeropts, ] <- matrix(Xproj[ref, ], byrow = TRUE, ncol = ncol(Xproj),
                                 nrow = sum(zeropts))
   }
   return(Xproj)
}

#version where dot products slightly above 1 are clamped to 1, so don't get NA
acos_clamped <- function(x){
  big1 <- abs(x) > 1
  stopifnot(all(abs(x[big1]) - 1 < .Machine$double.eps * 1E6))
  x[big1] <- sign(x)[big1]
  acos(x)
}

