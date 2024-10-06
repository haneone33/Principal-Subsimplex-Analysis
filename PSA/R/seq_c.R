# special function that produces a list of locations equidistant from each other and the half-distant to the edge - for annotation locations
seq_c <- function(from, to, length.out){
  celledges <- seq(from, to, length.out = length.out + 1)
  mids <- (celledges[-1] + celledges[-length(celledges)])/2
  return(mids)
}
