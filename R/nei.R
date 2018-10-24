#' Calculate the pairwise Nei distance
#' 
#' Given a `allele_count` matrix, it computes *Nei* distance and its bootstrap
#' replications between all pairs of population using different methods.
#' 
#' The estimator methods include two implementations of the *standard Nei's*
#' genetic distance (`fast` and `safe`), as described in *Nei (1972)*, and
#' *Nei's Da* (`da`), described in *Nei et al. (1983)*. Both Nei's distance
#' (*standard* and *Da*) measures genetic differences caused by mutation and
#' genetic drift with *Da* result in more reliable results for microsatellite
#' data, but more sensible to small population sizes.
#' 
#' Two versions of *standard Nei's distance* were implemented, `safe` and
#' `fast`. The `fast` implementation is 2 times faster than `safe`, but can
#' introduce errors when there is data missing, if no missing data is present
#' it results the same result as the `safe` implementation.
#' 
#' @name nei
#' @aliases NULL
NULL
#> NULL

#' @rdname nei
#' @export
nei_methods = c("fast", "safe", "da");
#> c("fast", "safe", "da")

#' @rdname nei
#' 
#' @param ac allele count matrix for which to compute nei's distance.
#' @param method nei's distance estimator method, with `fast` as default.
#' 
#' @examples
#' c <- matrix(sample(1:10, 250, replace = TRUE), ncol = 10, nrow = 25)
#' n <- matrix(10, ncol = 10, nrow = 25)
#' ac <- allele_count(c, n)
#' ## Summarized Nei
#' nei(ac)
#' nei(ac, "safe")
#' 
#' @export
nei <- function(ac, method = nei_methods) {
  if (class(ac) != "allele_count")
    stop("You must provide a `allele_count`")
  method <- match.arg(method)
  func <- switch(
    method,
    fast = function (f) {
      f[is.na(f)] <- 0
      return(.fast_nei(f))
    },
    safe = .nei,
    da = .nei_da)
  return(as.dist(func(allele_frequency(ac))))
}

#' @rdname nei
#' 
#' @param boots number of bootstrap replicates to perform
#' @param block statistic blocking size
#' @param as_matrix if `TRUE`, then return a matrix whith bootstrap interations
#' as rows and pairwise comparisons as columns.
#' 
#' @examples
#' bootstrap_nei(ac)
#' 
#' @export
bootstrap_nei <- function(
  ac, method = nei_methods,
  boots = 100L,
  block = 10L,
  as_matrix = FALSE) {
  if (class(ac) != "allele_count")
    stop("You must provide a `allele_count`")
  method <- match.arg(method)
  func <- switch(
    method,
    fast = function(f, nboots, block) {
      f[is.na(f)] <- 0
      .bootstrap_fast_nei
    },
    safe = function(f, nboots, block) {
      nei <- nei_matrix(f)
      bootstrap_nei(nei$jxy, nei$jx, block, nboots)
    },
    da = function(f, nboots, block) {
      nei <- nei_da_matrix(f)
      nbl <- ceiling(nrow(nei) / block)
      num <- .blocked_col_sums(nei, block, nbl)
      den <- .blocked_col_sums(!is.na(nei), block, nbl)
      .bootstrap_col_rates(num, den, nboots)
    })
  pops <- .pop_names(ac)
  M <- func(allele_frequency(ac), nboots, block)
  if (!as_matrix) {
    if (method == "fast"){
      colnames(M) <- rownames(M) <- pops
      M <- purrr::map(seq_len(nboots), ~ as.dist(M[,,.x]))
    } else {
      M <- purrr::map(seq_len(nboots), ~ .as_dist(M[.x,], pops))
    }
    return(M)
  }
  if (as_matrix && method == "fast") {
    ltri <- lower.tri(M[,,1])
    M <- matrix(M[ltri], ncol = sum(ltri), byrow = TRUE)
  }
  colnames(M) <- combn(pops, 2, function(x) paste(x, collapse = "/"))
  return(M)
}