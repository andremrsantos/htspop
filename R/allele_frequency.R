#' The allele frequency from the allele count
#' 
#' Given a allele count, it computes the allele frequency as numeric.
#' 
#' @param ac A allele count matrix to compute allele frequency
#' @param rows,cols A vector indicating the subset of rows and columns to be 
#'   evaluated. If TRUE, no subseting is done.
#' 
#' @export
allele_frequency <- function(ac, rows = TRUE, cols = TRUE) {
  if (class(ac) != "allele_count")
    stop("You must provide a allele count to compute allele frequency")
  subac <- ac[rows, cols]
  return(count_allele(subac) / total_allele(subac))
}