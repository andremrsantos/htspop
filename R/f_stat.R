#' Compute Reich's f-statistics
#' 
#' Computes *f2*, *f3*, and *f4* statistics unbiased estimateve used to measure
#' populations strutcture and treeness as described by
#' [Reich et al. 2009](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2842210/).
#' The *f-statistic* are described the following statistics:
#' \deqn{F4(W, X; Y,Z) = E[(W - X) * (Y - Z)]}
#' \deqn{F3(W; X, Y) = E[(W - X) * (W - Y)]}
#' \deqn{F2(X, Y) = E[(X - Y)^2]}
#' 
#' @param ac A allele count matrix including the populations to investigate.
#' @param pops A vector of populations to compute the statistic.
#' 
#' @name f_stat
NULL
#> NULL

#' @rdname  f_stat
#' @export
f4_stat <- function(ac, pops) {
  if (class(ac) != "allele_count")
    stop("You must provide a `allele_count` to use this function")
  if (!is.atomic(pops) || length(pops) != 4)
    stop("`pops` must be a vector of 4 entries W, X, Y, Z for f4(W, X; Y, Z)")
  ## Compute statistic
  freq <- allele_frequency(ac, cols = pops)
  return((freq[,1] - freq[,2]) * (freq[,3] - freq[,4]))
}

#' @rdname  f_stat
#' @export
f3_stat <- function(ac, pops) {
  if (class(ac) != "allele_count")
    stop("You must provide a `allele_count` to use this function")
  if (!is.atomic(pops) || length(pops) != 3)
    stop("`pops` must be a vector of 3 entries W, X, Y for f3(W; X, Y)")
  ## Compute statistic
  freq <- allele_frequency(ac, cols = pops)
  c <- count(ac[, pops[1]])
  n <- total(ac[, pops[1]])
  hc <- c * (n - c) / (n**3 - n)
  return((freq[,1]-freq[,2]) * (freq[,1]-freq[,3]) - hc)
}

#' @rdname  f_stat
#' @export
f2_stat <- function(ac, pops) {
  if (class(ac) != "allele_count")
    stop("You must provide a `allele_count` to use this function")
  if (!is.atomic(pops) || length(pops) != 2)
    stop("`pops` must be a vector of 2 entries X, Y for f2(X, Y)")
  ## Compute statistic
  freq <- allele_frequency(ac, cols = pops)
  c <- count(ac[, pops])
  n <- total(ac[, pops])
  hc <- matrixStats::rowSums2(c * (n - c) / (n**3 - n))
  return((freq[,1] - freq[,2])**2 - hc)
}