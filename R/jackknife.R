#' A blocked jacknife estimator
#' 
#' It performs a blocked jacknife estimator of mean and standard error of the
#' input vector. Given a matrix as input, it considers the first column as
#' numerator and second column as denominator. It can also perform a weighted
#' estimation by counting the number of non-NA values per block.
#' *This implementation is largely inspired on
#' (Treemix)[https://bitbucket.org/nygcresearch/treemix/wiki/Home]
#' implementation of the algorithm*.
#' 
#' @param x A numeric vector or matrix to estimate mean and standard error.
#' @param block A integer scalar indicating the size of blocks (greater
#'   than zero).
#' @param weighted A boolean scalar indicating if should perform a weighted
#'   estimation.
#'
#' @export
jackknife <- function(x, block = 10, weighted = TRUE) {
  if (!is.numeric(x))
    stop("`x` must be a numeric vector or a matrix with at two columns.")
  if (block < 1) 
    stop("block must be a integer greater than 0.")
  ## Add denominator when needed
  if ((is.matrix(x) && ncol(x) < 2) || !is.matrix(x))
    x <- cbind(x, 1)
  ## Compute blocks
  nb <- ceiling(nrow(x)/block)
  num <- .blocked_sums(x[,1], block, nb)
  den <- .blocked_sums(x[,2], block, nb)
  ## Run statistic
  numsum <- matrixStats::sum2(num)
  densum <- matrixStats::sum2(den)
  val <- numsum / densum
  jackvals <- (numsum - num) / (densum - den)
  jack_mean <- NA
  jack_stde <- NA
  if (weighted) {
    nsnv <- blocked_sum(!is.na(x[,1]), block)
    weight <- 1 - nsnv / matrixStats::sum2(nsnv)
    jack_mean <- nb * val - matrixStats::sum2(weight * jackvals)
    hj <- matrixStats::sum2(nsnv) / nsnv
    jack_stde <- (hj * val - (hj - 1) * jackvals - jack_mean)**2 / (hj - 1)
    jack_stde <- sqrt(matrixStats::mean2(jack_stde, na.rm = TRUE))
  } else {
    psedovals <- nb * val - (nb - 1) * jackvals
    jack_mean <- matrixStats::mean2(psedovals)
    jack_stde <- matrixStats::sum2((pseudovals - jack_mean)**2)
    jack_stde <- sqrt(jack_stde / (nb * (nb-1)))
  }
  return(list(
    n = nrow(x),
    n_blocks = nb,
    mean = jack_mean,
    se = jack_stde,
    z = jack_mean / jack_stde))
}