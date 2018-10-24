## Private functions
rsum_ <- function(mtx, idx)
  matrixStats::rowSums2(mtx, cols = idx, na.rm = TRUE)

map_matrix_ <- function(...) 
  as.matrix(dplyr::bind_cols(purrr::map(...)))

#' Allele Count
#' 
#' Creates an allele count matrix used on other analysis. The structure
#' encapsulates the allele count and number of individual with variant as rows
#' and groups as columns. It allows subsetting (`[]`), to check dimensions 
#' (with `dim`, `nrow`, and `ncol`), and renaming them.
#'  
#' @param count Matrix giving the number of alternative (or mutant) alleles 
#' with variants as rows and samples as columns.
#' @param n Matrix giving the total number of alleles observed with variants as 
#' rows and samples as columns.
#' 
#' @return Allele count matrix
#' 
#' @examples
#' ac <- allele_count(matrix(0:5, ncol = 2), matrix(rep(6, 6), ncol = 2))
#' 
#' print(ac[1:2, ]) # allows subsetting
#' 
#' nrow(ac)         # number of variants
#' rownames(ac)     # variants identifier
#' 
#' ncol(ac)         # number of groups
#' colnames(ac)     # groups names
#' @export
allele_count <- function(count, n) {
  ## Sanity check data
  if (!is.numeric(count) || !is.matrix(count))
    stop("`count` must be a numeric matrix.")
  if (!is.numeric(n) || !is.matrix(n))
    stop("`n` must be a numeric matrix.")
  if (any(dim(count) != dim(n)))
    stop("different dimensions found.")
  if (any(count > n))
    stop("`n` must be greater than `count`.")
  ## Instanciate class
  structure(list(count=count, n=n), class="allele_count")
}

#' @rdname allele_count
#' 
#' @param genotype Sample genotype matrix with the variants as rows and samples
#' as column in the format `012`, where `0` is the reference homozygous, `1`
#' the heterozygous and `2` the alternative homozygous.
#' @param groups List of sample vector or a single vector of groups used to
#' compute the allele count. When is *NULL* (default value) it considers all
#' samples as a single group, when is either a vector it considers only the
#' corresponding samples, and when is a list it computes the allele count for
#' each sub-vector.
#' 
#' @examples
#' ## From genotype matrix
#' geno_sample <- matrix(sample(0:2, 100, replace = TRUE), ncol = 10)
#' ac_matrix_from_genotype(geno_sample)      ## all samples
#' ac_matrix_from_genotype(geno_sample, 1:5) ## only the first 5
#' ## for three groups
#' ac_matrix_from_genotype(
#'   geno_sample,
#'   list(a = 1:5, b = 3:7, c = 6:10))
#' 
#' @export
allele_count_from_genotype <- function(genotype, groups = NULL) {
  if (!is.numeric(genotype) || !is.matrix(genotype))
    stop("`genotype` must be a numeric matrix.")
  ## Given a list of indexes
  if (is.list(groups)) {
    ac <- list(
      count = map_matrix_(groups, rsum_, mtx = genotype),
      n = 2 * map_matrix_(groups, rsum_, mtx = !is.na(genotype)))
    return(structure(ac, class = "allele_count"))
  }
  ## Default
  ac <- list(
    count = t(rsum_(genotype, groups)),
    n = 2 * t(rsum_(!is.na(genotype), groups)))
  return(structure(ac, class = "allele_count"))
}

#' @export
`[.allele_count` <- function(ac, i, j) {
  subac <- list(
    count = ac$count[i, j, drop=FALSE], 
    n = ac$n[i, j, drop=FALSE])
  structure(subac, class = "allele_count")
}

#' @export
dim.allele_count <- function(ac) {
  return(dim(ac$count))
}

#' @export
dimnames.allele_count <- function(ac) {
  return(dimnames(ac$count))
}

#' @export
`dimnames<-.allele_count` <- function(ac, ...) {
  `dimnames<-`(ac$count, ...)
  `dimnames<-`(ac$n, ...)
}

#' @export
print.allele_count <- function(ac, n = 10) {
  dim <- dim(ac)
  rowidx <- seq_len(pmin(n, dim[1]))
  colidx <- seq_len(dim[2])
  submtx <- ac[rowidx, colidx]
  
  to_print <- matrix(paste0(submtx$count, "/", submtx$n), ncol = ncol(submtx))
  dimnames(to_print) <- dimnames(submtx)
  
  cat(paste0("Allele Count Matrix (", dim[1], ", ", dim[2], ")\n"))
  print(to_print)
  if (length(rowidx) < dim[1]) cat("[...]\n")
}