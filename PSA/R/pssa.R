#' @title Principle Sub Orthant Analysis for the Simplex
#' @inheritParams recurse_tolowerdimension
pssa <- function(X, V = diag(ncol(X)), testweights = seq(0, 1, length.out = 100),
                 merges = NULL, maxsteps = min(ncol(X), nrow(merges)) - 1){
  results <- recurse_tolowerdimension(X, V = V,
                                       evaluator = evaluatesubsimplex,
                                       projector = projectsubsimplex,
                                       testweights = testweights,
                                       merges = merges,
                                       maxsteps = maxsteps
  )
  return(results)
}
