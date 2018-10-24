#include "blocked_stats.hpp"
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::uvec blocked_sample(const uint bsize, const uint nblocks, const uint length) {
  arma::umat mtx = arma::reshape(
    arma::regspace<arma::uvec>(0, length), bsize, nblocks);
  arma::uvec sam = arma::randi<arma::uvec>(
    nblocks, arma::distr_param(0, nblocks-1));
  return(arma::vectorise(mtx.cols(sam)));
}

// [[Rcpp::export(name=".blocked_sums")]]
arma::vec blocked_sums(
    const arma::vec& val,
    const uint bsize,
    const uint nblock,
    const bool rm_na) {
  const uint limit = std::min(val.size(), bsize * nblock);

  arma::vec aux = arma::vec(nblock, arma::fill::zeros);
  for (uint i = 0; i < limit; i++) {
    if (rm_na && !arma::is_finite(val(i)))
      continue;
    aux(i/bsize) += val(i);
  }
  return(aux);
}

// [[Rcpp::export(name=".blocked_means")]]
arma::vec blocked_means(
    const arma::vec& val,
    const uint bsize,
    const uint nblock,
    const bool rm_na) {
  const uint limit = std::min(val.size(), bsize * nblock);

  arma::vec num = arma::vec(nblock, arma::fill::zeros);
  arma::vec den = arma::vec(nblock, arma::fill::zeros);
  for (uint i = 0; i < limit; i++) {
    if (rm_na && !arma::is_finite(val(i)))
      continue;
    num(i/bsize) += val(i);
    den(i/bsize) += 1;
  }
  return(num/den);
}

// [[Rcpp::export(name=".blocked_col_sums")]]
arma::mat blocked_col_sums(
    const arma::mat& matrix,
    const uint bsize,
    const uint nblock,
    const bool rm_na) {
  const uint limit = std::min(matrix.n_rows, bsize * nblock);

  arma::rowvec aux;
  arma::mat blocked(nblock, matrix.n_cols, arma::fill::zeros);
  for (uint i = 0; i < limit; i++) {
    aux = matrix.row(i);
    if (rm_na)
      aux(arma::find_nonfinite(aux)).zeros();
    blocked.row(i/bsize) += aux;
  }
  return(blocked);
}

// [[Rcpp::export(name=".blocked_col_means")]]
arma::mat blocked_col_means(
    const arma::mat& matrix,
    const uint bsize,
    const uint nblock,
    const bool rm_na) {
  const uint limit = std::min(matrix.n_rows, bsize * nblock);

  arma::rowvec aux_num, aux_den;
  arma::uvec is_na;
  arma::mat num(nblock, matrix.n_cols, arma::fill::zeros);
  arma::mat den(nblock, matrix.n_cols, arma::fill::zeros);
  for (int i = 0; i < limit; i++) {
    aux_num = matrix.row(i);
    aux_den = arma::rowvec(matrix.n_cols, arma::fill::ones);
    if (rm_na) {
      is_na = arma::find_nonfinite(aux_num);
      aux_num(is_na).zeros();
      aux_den(is_na).zeros();
    }
    num.row(i/bsize) += aux_num;
    den.row(i/bsize) += aux_den;
  }
  return(num/den);
}