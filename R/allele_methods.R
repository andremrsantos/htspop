#' Allele Count Manipulation
#' 
#' Functions to explore and manipule the data from allele count matrices.
#' 
#' @param ac Allele count matrix to analyze
#' 
#' @name allele_methods
#' @aliases NULL
NULL
#> NULL

#' @rdname allele_methods
#' @export
count_allele <- function(ac) {
  if (class(ac) != "allele_count")
    stop("`ac` must be an `allele_count`.")
  return(ac$count)
}

#' @rdname allele_methods
#' @export
total_allele <- function(ac) {
  if (class(ac) != "allele_count")
    stop("`ac` must be an `allele_count`.")
  return(ac$n)
}


#' @rdname allele_methods
#' @export
invert_allele <- function(ac) {
  if (class(ac) != "allele_count")
    stop("`ac` must be an `allele_count`.")
  new_mtx <- list(count = ac$n - ac$count, n = ac$n)
  return(structure(new_mtx, class = "allele_count"))
}

#' @rdname allele_methods
#' @export
population_names <- function(ac) {
  if (class(ac) != "allele_count")
    stop("`ac` must be an `allele_count`.")
  names <- colnames(ac)
  if (is.null(names))
    return(seq_len(ncol(ac)))
  return(names)
}
