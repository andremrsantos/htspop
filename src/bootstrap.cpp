#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "_tracker.h"

// [[Rcpp::export(name=".bootstrap_col_rates")]]
arma::mat bootstrap_col_rates(
    const arma::mat& num,
    const arma::mat& den,
    const uint nboots,
    const bool rm_na = false) {
  const uint nrow = num.n_rows;
  const uint ncol = num.n_cols;
  // Clean NA
  arma::mat aux_num = num;
  arma::mat aux_den = den;
  if (rm_na) {
    arma::uvec is_na;
    // remove na
    is_na = arma::find_nonfinite(num + den);
    aux_num.elem(is_na).zeros();
    aux_den.elem(is_na).zeros();
    // remove den == 0
    is_na = arma::find(den == 0);
    aux_num.elem(is_na).zeros();
    aux_den.elem(is_na).zeros();
  }

  Tracker t(nboots);
  arma::uvec index;
  arma::mat simu = arma::mat(nboots, ncol);
  for (uint i = 0; i < nboots; i++) {
    index = arma::randi<arma::uvec>(nrow, arma::distr_param(0, nrow-1));
    simu.row(i) =
      arma::sum(aux_num.rows(index), 0) / arma::sum(aux_den.rows(index), 0);
    t.it();
  }
  return(simu);
}

// [[Rcpp::export(name=".bootstrap_col_means")]]
arma::mat bootstrap_col_means(
    const arma::mat& val,
    const uint nboots,
    const bool rm_na = true) {
  const int nrow = val.n_rows;
  const int ncol = val.n_cols;

  arma::uvec is_na;
  arma::mat num = val;
  arma::mat den = arma::mat(nrow, ncol, arma::fill::ones);
  if (rm_na) {
    is_na = arma::find_nonfinite(val);
    num.elem(is_na).zeros();
    den.elem(is_na).zeros();
  }

  Tracker t(nboots);
  arma::uvec index;
  arma::mat simu = arma::mat(nboots, ncol);
  for (uint i = 0; i < nboots; i++) {
    index = arma::randi<arma::uvec>(nrow, arma::distr_param(0, nrow-1));
    simu.row(i) =
      arma::sum(num.rows(index), 0) / arma::sum(den.rows(index), 0);
    t.it();
  }
  return(simu);
}