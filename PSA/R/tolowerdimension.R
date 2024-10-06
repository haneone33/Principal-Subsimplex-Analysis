#' @title Finds and applies best merge
#' @param X matrix of data points for a simplex with vertices at the basis coordinates. Each row an observation.
#' @param V matrix of vertices in the basis of the *original* data. Each row a vertex corresponding to the *columns* of `X`.
#' @param merge If supplied, used instead of [`gridsearch()`]. `bestmerge` must have same structure as result from [`gridsearch()`]
#' @inheritParams gridsearch
#' @param projector A function that takes result of `gridsearch()` and applies a projection. Returning `X`, `V`, `scores`, `mergedirection`.
#' @return A list
#' + `X` Projected data points in with the new vertices as the coordinates
#' + `V` New set of vertices with row names.
#' + `scores` Signed scores for the merge
#' + `mergedirection` The unit vector parallel to `v2 - v1` and in the direction of positive scores, and in the original coordinates.
#' + `v1` The index of `v1` used in `projector` (corresponding to a row of the old vertices matrix).
#' + `v2` The index of `v2` used in `projector` (corresponding to a row of the old vertices matrix).
#' + `w` The weight of the merge between `v1` and `v2`. As in `w * v1 + (1-w) * v2`, where `v1` and `v2` are the first and second vertices given by `mergedvertices`.
#' @details If `merge` is not supplied and the *largest* score of the found merge is *negative*, then swaps the signs of the scores and records this by swapping the returned `v1` with `v2`, changing w to 1-w, and flipping the direction of `mergedirection`. This should mean that the returned merge information can be reapplied on new data and the direction will be the same.
#'
#' If `merge` is non-`NULL` then `evaluator` is ignored - no search occurs.
#' @noRd
tolowerdimension <- function(X, V = diag(1, ncol(X)), evaluator, projector, testweights = seq(0, 1, length.out = 100), merge = NULL){
  findmerge <- is.null(merge)
  if (findmerge){
    merge <- gridsearch(X = X, evaluator = evaluator, testweights = testweights)
  }
  lowerD <- projector(X = X, V = V, v1 = merge$v1, v2 = merge$v2, w = merge$w)
  if (findmerge){#reorder v1 and v2 so that largest score is positive
     flipdirection <- (sign(lowerD$scores[which.max(abs(lowerD$scores))]) < 0)
     if (flipdirection){
       lowerD$scores <- -1 * lowerD$scores
       lowerD$mergedirection <- -1 * lowerD$mergedirection
       merge[c("v1", "v2")] <- rev(merge[c("v1", "v2")])
       merge["w"] <- 1 - merge["w"]
     }
  }

  return(list(
    X = lowerD$X,
    V = lowerD$V,
    scores = drop(lowerD$scores),
    mergedirection = lowerD$mergedirection,
    v1 = merge$v1,
    v2 = merge$v2,
    w = merge$w
  ))
}

