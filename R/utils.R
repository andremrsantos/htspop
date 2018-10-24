# Wrapper to convert a distance vector into `dist`
.as_dist <- function(vec, pops) {
  D <- matrix(0, length(pops), length(pops))
  D[lower.tri(D)] <- vec
  return(as.dist(D))
}