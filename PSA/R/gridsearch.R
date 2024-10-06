#' @title Find the best merge by grid search
#' @param X matrix of data points for a simplex with vertices at the basis coordinates.
#' @inheritParams gridsearch_pair
#' @return The best pair and weight as a single row data frame (data frame for easier row binding and because the v1 and v2 are best stored as integers).
#' @noRd
gridsearch <- function(X, evaluator, testweights = seq(0, 1, length.out = 100), ...){
  if (ncol(X) == 2){
    pairs <- matrix(as.integer(c(1,2)), nrow = 2)
    res <- gridsearch_pair(X, pairs[1], pairs[2], testweights = testweights, evaluator = evaluator, ...)
    res <- t(as.matrix(res))
  } else {
    pairs <- utils::combn(ncol(X), 2)
    res <- apply(pairs, MARGIN = 2, function(v1v2){gridsearch_pair(X, v1v2[1], v1v2[2], testweights = testweights, evaluator = evaluator, ...)})
    res <- do.call(rbind, res)
  }
  best <- which.min(res[,"value"])
  rownames(pairs) <- c("v1", "v2")
  out <- data.frame(v1 = pairs[1, best], v2 = pairs[2, best], w = unlist(res[best, "w"]))
  return(out)
}


#' @title Find the best weight for merging two vertices by grid search
#' @param X matrix of data points for a simplex with vertices at the basis coordinates.
#' @param v1 an index of a vertex (a basis vector)
#' @param v2 an index of a vertex (a basis vector)
#' @param testweights A vector of weights to try.
#' @param evaluator A function that evaluates each proposed merge and return a single non-negative value for the quality of the merge. Arguments of `X`, `v1`, `v2`, and `w` at least
#' @param ... Passed to `evaluator`
#' @noRd
gridsearch_pair <- function(X, v1, v2,
                             testweights = seq(0, 1, length.out = 100),
                             evaluator,
                             ...){
  objvals <- vapply(testweights, function(w){
         evaluator(X=X, v1 = v1, v2 = v2, w = w, ...)},
         FUN.VALUE = 1.2)
  bestobjval <- min(objvals)
  bestw <- testweights[which.min(objvals)]
  return(list(
    w = bestw,
    value = bestobjval
  ))
}
