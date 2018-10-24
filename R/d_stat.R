#' Compute D-statistics
#' 
#' Computes the unbiased estimative of *D-statististic*, also known as
#' 4-population test, used to measure populations strutcture as described by
#' [Patterson et al. 2012](http://www.genetics.org/content/192/3/1065).
#' Giving 4 populations (*W*, *X*, *Y*, and *Z*), *D-statistic* compares the
#' probability of a *BABA* event, meaning *W* and *Y* agree and *X* and *Z*
#' agree, to the probability of a *ABBA* event, meaning *W* and *Z* agree.
#' When *D-statistic* present positive values it favors *BABA*, when negative
#' present *ABBA*, and zero favors *AABB* (*W* and *X* agree). The
#' *D-statistic* describe the following equation:
#' \deqn{D = (P(BABA) - P(ABBA)) / (P(BABA) + P(ABBA))}
#' 
#' @param ac A allele count matrix including the populations to investigate.
#' @param pops A vector of populations to compute the statistic.
#' 
#' @export
d_stat <- function(ac, pops) {
  fq <- allele_frequency(ac, cols = pops)
  ## Compute values
  num <- (fq[,1] - fq[,2]) * (fq[,3] - fq[,4])
  dem <-
    (fq[,1] + fq[,2] - 2 * fq[,1] * fq[,2]) *
    (fq[,3] + fq[,4] - 2 * fq[,3] * fq[,4])
  return(cbind(num, dem))
}