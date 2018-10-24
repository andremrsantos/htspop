#' Filter Variants from Allele Count
#' 
#' Implements common strategies to filter variants of an allele frequency
#' matrix, including `complete_variants`, which exclude variants with missing
#' genotype, and `informative_variants`, which excludes variants were allgroups
#' are reference or alternative homozygous and keeping only informative sites.
#' These functions return a `boolean` vector of the subset and can also be used
#' with `filter_variants` to return the subsetted allele count or frequency
#' matrix.
#' 
#' @param x A allele count or frequency matrix to be filtered.
#' @param cols A vector indicating the subset of columns considered.
#' 
#' @export
filter_variants <- function(x, fn, cols = NULL){
  if (!is.function(fn))
    stop("`fn` must be a function returning a boolean vector")
  return(x[fn(x, cols = cols),])
}

#' @rdname  filter_variants
#' @export
complete_variants <- function(x, cols = NULL) {
  if (class(x) == "allele_count")
    return(matrixStats::rowAlls(total(x) > 0, cols = cols))
  if (is.matrix(x) && is.double(x))
    return(!matrixStats::rowAnyNAs(x, cols = TRUE))
  stop("You must provide a allele count or frequency matrix to filter.")
}

#' @rdname filter_variants
#' @export
informative_variants <- function(x, cols = NULL) {
  if (class(x) == "allele_count") {
    rmeans <- matrixStats::rowMeans2(count(x), cols = cols, na.rm = TRUE)
    return(rmeans != 0 & rmeans != 2)
  } 
  if (is.matrix(x) && is.double(x)){
    rmeans <- matrixStats::rowMeans2(x, cols = cols, na.rm = TRUE)
    return(rmeans != 0 & rmeans != 1)
  }
  stop("You must provide a allele count or frequency matrix to filter.")
}