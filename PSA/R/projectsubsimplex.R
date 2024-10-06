#' @title Project to a Subsimplex
#' @details Following Version 4 of the descriptions.
#' The result of the projection are points with mass (i.e. value) in the newly created vertex equal to the sum of the mass (i.e. value) in the components specified by `v1` and `v2`.
#' The new vertex is computed as `w * V[v1, ] + (1-w) * V[v2, ]`, which is the convex combination of two vertices and is *not* a vector with L2-norm equal to 1.
#' Scores a signed as positive if the difference between the input point and the projected point has a positive inner-product with the merge direction `V[v2, ] - V[v1, ]`.
#' @inheritParams projectsuborthant
#' @noRd
projectsubsimplex <- function(X, V = diag(ncol(X)), v1, v2, w){
  #get new points without lowering dimension
  Xv1v2 <- X[, v1] + X[, v2]
  newXtmp <- X
  newXtmp[, v1] <- w * Xv1v2
  newXtmp[ , v2] <- (1-w) * Xv1v2

  #scores
  sscores <- sqrt(2) * (w * X[, v2] - (1-w) * X[, v1]) #skipping sqrt(^2) means sign included here :)

  #projected X expressed in lower dimension basis vector
  newX <- cbind(Xv1v2, X[, -c(v1, v2), drop = FALSE])
  colnames(newX)[[1]] <- paste0("", sprintf("(%s,%s)", colnames(X)[v1],  colnames(X)[v2])) #the paste0("", ) is fix so a string of "" comes out if colnames are NULL, rather than forming an error

  #merge direction
  mergedirection <- V[v2, ] - V[v1, ]
  #record new vertex in V
  newV <- rbind(w * V[v1, ] + (1-w) * V[v2, ],
                 V[-c(v1, v2), , drop = FALSE]
               )

  return(list(
    X = newX,
    V = newV,
    scores = sscores,
    mergedirection = mergedirection
  ))
}

# A function that evaluates each proposed merge and returns a single non-negative value for the quality of the merge. Arguments of `X`, `v1`, `v2`, and `w` at least.
evaluatesubsimplex <- function(X, v1, v2, w){
  proj <- projectsubsimplex(X = X, V = diag(ncol(X)), v1 = v1, v2 = v2, w = w)
  sum(proj$scores^2)
}
