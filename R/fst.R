#' Calculate the pairwise Fst statistic
#' 
#' Given a `allele_count` matrix, it computes the *Fst* statistic and its
#' bootstrap replications between all pairs of population using several
#' methods.
#' 
#' The estimator methods include *Reich et al. (2019)* (`reich`),
#' *Weir & Cockerham (1984)* (`weir_cockerham` or `wc`), *Hudson et al. (1992)*
#' (`hudson`) and *Wright (1951)* (`wright`). These methods use different
#' approaches to estimate *Fst*. The *Wright* method is a biased estimator and
#' can be affected by different populations sizes, *Hudson* and *Reich* methods
#' are unbiased fast and tend to present reliable results on large variant
#' sizes, but limited precision snpwise. In other hand, *Weir & Cockerham*
#' method is the most common estimator (implemented by tools such as *plink*,
#' *vcftools*, and *scikit-allele*) and uses a population size correction to
#' reduce bias, but is a little slower than the other methods.
#' 
#' Please note that *Reich* method requires population sizes larger than 2 (
#' more than one diploid individual).
#' 
#' @name fst
#' @aliases NULL
NULL
#> NULL
#' @rdname fst
#' @export
fst_methods = c("reich", "weir_cockerham", "wc", "hudson", "wright")
#> c("reich", "weir_cockerham", "wc", "hudson", "wright")
# Validate and identify the fst function
.fst_func <- function(method=fst_methods) {
  func <- switch(
    match.arg(method),
    reich = .pairwise_fst_reich,
    weir_cockerham = .pairwise_fst_weir_cockerham,
    wc = .pairwise_fst_weir_cockerham,
    hudson = .pairwise_fst_hudson,
    wright = .pairwise_fst_wright)
  return(func)
}
#' @rdname fst
#' 
#' @param ac allele count matrix for which to compute `Fst`.
#' @param method fst estimator method. Using `reich` as default.
#' @param snpwise a boolean scalar idicating if should return a the statistic
#' for each variant (when `TRUE`) or the summary (when `FALSE`).
#' 
#' @examples
#' c <- matrix(sample(1:10, 250, replace = TRUE), ncol = 10, nrow = 25)
#' n <- matrix(10, ncol = 10, nrow = 25)
#' ac <- allele_count(c, n)
#' ## Summarized Fst
#' fst(ac)
#' fst(ac, "wc")
#' ## Snpwise Fst
#' fst(ac, snpwise = TRUE)
#' fst(ac, "wc", snpwise = TRUE)
#' 
#' @export
fst <- function(ac, method = fst_methods, snpwise = FALSE) {
  if (class(ac) != "allele_count")
    stop("You must provide a `allele_count` to compute fst.")
  if (!purrr::is_scalar_logical(snpwise))
    stop("`snpwise` must be a logical scalar")
  ## Extract population names
  pnames <- population_names(ac)
  ## Process output
  M <- .fst_func(method)(count(ac), total(ac))
  if (snpwise) {
    if (is.list(M))
      M <- M$den / M$num
    colnames(M) <- combn(pnames, 2, function(x) paste(x, collapse = "/"))
    rownames(M) <- rownames(ac)
    return(M)
  }
  if (is.list(M)) {
    M$den[M$num <= 0] <- 0
    M$num[M$num <= 0] <- 0
    M <-
      matrixStats::colSums2(M$den, na.rm = TRUE) /
      matrixStats::colSums2(M$num, na.rm = TRUE)
  } else {
    M <- matrixStats::colMeans2(M, na.rm = TRUE)
  }
  return(.as_dist(M, pnames))
}
#' @rdname fst
#' 
#' @param boots number of bootstrap replicates to perform
#' @param block statistic blocking size
#' @param as_matrix if `TRUE`, then return a matrix whith bootstrap interations
#' as rows and pairwise comparisons as columns.
#' 
#' @examples
#' ## bootstrap fst
#' bootstrap_fst(ac)
#' 
#' @export
bootstrap_fst <- function(
    ac, method = fst_methods,
    boots = 100L,
    block = 10L,
    as_matrix = FALSE) {
  ## Check input and validate input
  if (class(ac) != "allele_count")
    stop("You must provide a `allele_count` to compute fst.")
  if (!purrr::is_scalar_logical(as_matrix))
    stop("`as_matrix` must be a logical scalar")
  if (!purrr::is_scalar_numeric(boots) || !purrr::is_scalar_numeric(block))
    stop("`boots` and `block` must be a logical scalar")
  
  pnames <- population_names(ac)
  nblock <- ceiling(nrow(ac)/block)
  M <- .fst_func(method)(count(ac), total(ac))
  num <- NULL
  den <- NULL
  if (is.list(M)) { ## Deal with reich return
    num <- .blocked_col_sums(M$num, block, nblock)
    den <- .blocked_col_sums(M$den, block, nblock)
  } else {
    num <- .blocked_col_sums(M, block, nblock)
    den <- .blocked_col_sums(!is.na(M), block, nblock)
  }
  simu = .bootstrap_col_rates(num, den, boots)
  if (as_matrix) {
    colnames(simu) <- combn(pnames, 2, function(x) paste(x, collapse = "/"))
    return(simu)
  }
  return(purrr::map(seq_len(boots), function(i) .as_dist(simu[i,], pnames)))
}