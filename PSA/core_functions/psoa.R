#' @title Principle Sub Orthant Analysis for the Simplex
#' @inheritParams recurse_tolowerdimension
#'
#' @export
psoa <- function(X, V = diag(ncol(X)), testweights = seq(0, 1, length.out = 100),
                 merges = NULL, maxsteps = min(ncol(X), nrow(merges)) - 1){
  results <- recurse_tolowerdimension(X, V = V,
                                       evaluator = evaluatesuborthant,
                                       projector = projectsuborthant,
                                       testweights = testweights,
                                       merges = merges,
                                       maxsteps = maxsteps
  )
  return(results)
}
