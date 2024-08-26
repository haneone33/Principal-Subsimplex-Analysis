#' @title Recurse through dimensions
#' @description Given a merge `evaluator` and a merge `projector`, recurse through dimensions keeping track of the appropriate data.
#' @inheritParams tolowerdimension
#' @param merges If non-NULL, then a data frame of prespecified merges.
#' @details
#' + `Vbirth` keeps track of the merge that created each vertex. The merge is recorded acording to the destination dimension of the merge, so a value of `i` means the merge that created that vertex was creating the data for dimension `d - (i-1)`, where `d` is the starting dimension.
#' @export
recurse_tolowerdimension <- function(X, V = diag(ncol(X)),
          evaluator, projector,
          testweights = seq(0, 1, length.out = 100),
          merges = NULL,
          maxsteps = NULL
          ){
  X <- as.matrix(X)
  stopifnot(all(X >= 0))
  stopifnot(all(abs(rowSums(X) - 1) < 1E-10))
  # prepare vertex names
  if (all(V == diag(ncol(X)))){
    if (!is.null(colnames(X))){colnames(V) <- rownames(V) <- colnames(X)}
    else if (!is.null(colnames(V))){colnames(X) <- rownames(V) <- colnames(V)}
    else if (!is.null(rownames(V))){colnames(X) <- colnames(V) <- rownames(V)}
    else {colnames(X) <- colnames(V) <- rownames(V) <- paste0("V", 1:nrow(V))}
  } else {
    if (!is.null(colnames(X))){rownames(V) <- colnames(X)}
    else if (!is.null(rownames(V))){colnames(X) <- rownames(V)}
    else {colnames(X) <- rownames(V) <- paste0("V", 1:nrow(V))}
  }

  #computing maxsteps here rather than in argument to avoid lazy computation after I've modified merges already
  if (is.null(maxsteps)){maxsteps <- min(ncol(X), nrow(merges)) - 1}

  #initialise recursion, in the below lists content for slot i is either the value for d - (i-1) dimension or a value of the transition *to* d - (i-1) dimension, with the first slot being zero for these transition values.
  pts <- scores <- vertices <- modes <- list()
  pts[[1]] <- X
  vertices[[1]] <- V
  scores[[1]] <- rep(NA_real_, nrow(X))
  modes[[1]] <- matrix(NA_real_, nrow = nrow(X), ncol = ncol(X))
  inmerges <- merges
  merges <- data.frame(v1 = NA_integer_, v2 = NA_integer_, w = NA_real_)
  Vbirth <- list()
  Vbirth[[1]] <- -1 * (1:nrow(V)) #Vertices that started have a 'birth' of -1 * index according to `hclust()` help
  names(Vbirth[[1]]) <- rownames(V)
  mergedirections <- matrix(NA, nrow = 1+maxsteps, ncol = ncol(V))
  colnames(mergedirections) <- colnames(V)

  for (i in 1 + (1:maxsteps)){#i is the dimension less
    res <- tolowerdimension(X = pts[[i-1]],
                            V = vertices[[i-1]],
                            evaluator = evaluator,
                            projector = projector,
                            testweights = testweights,
                            merge = inmerges[i, ])
    pts[[i]] <- res$X
    vertices[[i]] <- res$V
    scores[[i]] <- res$scores
    mergedirections[i, ] <- t(res$mergedirection)
    modes[[i]] <- res$scores %x% t(res$mergedirection)
    colnames(modes[[i]]) <- names(res$mergedirection)
    merges[i, ] <- data.frame(v1 = res$v1, v2 = res$v2, w = res$w)
    Vbirth[[i]] <- c(i, Vbirth[[i-1]][-c(res$v1, res$v2)]) #assumes the first vertex is the new one, the i means values refer to the d - (i-1) destination dimension of the merge that created the vertex.
    print(sprintf("Finished creating dimension %i.", ncol(X) - i + 1))
  }

  names(pts) <- paste0("r=", seq(ncol(X), by = -1, length.out = maxsteps+1))
  names(vertices) <- names(scores) <- names(modes) <- names(Vbirth) <-  names(pts)
  rownames(merges) <- rownames(mergedirections) <- names(pts)

  return(list(
    pts = pts,
    vertices = vertices,
    scores = scores,
    modes = modes,
    merges = merges,
    mergedirections = mergedirections,
    Vbirth = Vbirth
  ))
}
