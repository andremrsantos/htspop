// [[Rcpp::depends(RcppArmadillo)]]
#include "RcppArmadillo.h"

#include "blocked_stats.hpp"
#include "_tracker.h"

arma::mat crossprod(const arma::mat& x) {
  return(arma::trans(x) * x);
}

arma::mat nei_jx(const arma::mat& f) {
  return(1 + 2 * (arma::square(f) - f));
}

arma::vec nei_jxy(const arma::vec& fx, const arma::vec& fy) {
  return(1 - fx - fy + 2 * fx % fy);
}

arma::mat nei_jxy(const arma::mat& f) {
  uint vi = 0;
  arma::mat jxy(f.n_rows, f.n_cols * (f.n_cols-1)/2);
  for (uint i = 0; i < f.n_cols; i++) {
    for (uint j = i + 1; j < f.n_cols; j++) {
      jxy.col(vi) = nei_jxy(f.col(i), f.col(j));
      vi++;
    }
  }
  return(jxy);
}
// [[Rcpp::export(name=".nei_matrix")]]
Rcpp::List nei_matrix(const arma::mat& f) {
  arma::mat jx;
  return(Rcpp::List::create(
      Rcpp::Named("jx") = nei_jx(f),
      Rcpp::Named("jxy") = nei_jxy(f)));
}

arma::mat nei(const arma::mat& jxy, const arma::mat& jx) {
  arma::uvec non_na;
  arma::vec vjx, vjy, vjxy;
  arma::mat dist(jx.n_cols, jx.n_cols);
  for(uint i = 0; i < jx.n_cols; i++) {
    for(uint j = i + 1; j < jx.n_cols; j++) {
      vjxy = jxy.col((j-i) + (2 * jx.n_cols * i - i*(i+1))/2 - 1);
      vjx = jx.col(i);
      vjy = jx.col(j);

      non_na = arma::find_finite(jxy);
      dist(i, j) =
        log(arma::sum(vjx(non_na)))/2 +
        log(arma::sum(vjy(non_na)))/2 -
        log(arma::sum(jxy(non_na)));
      dist(j, i) = dist(i, j);
    }
    dist(i, i) = 0;
  }
  return(dist);
}

// [[Rcpp::export(name=".nei")]]
arma::mat nei(const arma::mat& f) {
  return(nei(f, nei_jx(f)));
}

// [[Rcpp::export(name=".nei_da_matrix")]]
arma::mat nei_da_matrix(const arma::mat& f) {
  uint vi;
  arma::mat mat = arma::mat(f.n_rows, f.n_cols * (f.n_cols-1)/2);
  for (uint i = 0; i < f.n_cols; i++) {
    for (uint j = i + 1; j < f.n_cols; j++) {
      vi = ((j-i) + (2 * f.n_cols * i - i*(i+1))/2 - 1);
      mat.col(vi) = arma::sqrt(f.col(i) % f.col(j)) +
        arma::sqrt((1 - f.col(i)) % (1 - f.col(j)));
    }
  }
  return(mat);
}

// [[Rcpp::export(name=".nei_da")]]
arma::mat nei_da(const arma::mat& f) {
  arma::mat da = nei_da_matrix(f);
  arma::vec aux;
  arma::mat dist(f.n_cols, f.n_cols);
  for (int i = 0; i < f.n_cols; i++) {
    for (int j = i + 1; j < f.n_cols; j++) {
      aux = da.col((j-i) + (2 * f.n_cols * i - i*(i+1))/2 - 1);
      dist(i, j) = 1 - arma::mean(aux(arma::find_finite(aux)));
      dist(j, i) = dist(i, j);
    }
  }
  return(dist);
}

arma::mat fast_nei(const arma::mat& f, const arma::mat& jmt) {
  arma::mat num, den;

  num = arma::log(crossprod(f) + crossprod(1-f));
  den = arma::log(crossprod(arma::sum(jmt, 0)))/2;
  return(den - num);
}

// [[Rcpp::export(name=".fast_nei")]]
arma::mat fast_nei(const arma::mat& f) {
  return(fast_nei(f, nei_jx(f)));
}

// [[Rcpp::export(name=".bootstrap_fast_nei")]]
arma::cube bootstrap_fast_nei(
    const arma::mat& f,
    const uint boots,
    const uint block) {
  const uint nrow = f.n_rows;
  const uint ncol = f.n_cols;
  const uint nblk = std::ceil((double) nrow / block);
  const arma::mat jmt = nei_jx(f);

  Tracker t(boots);
  arma::uvec sample;
  arma::cube bootstrap(ncol, ncol, boots);
  for (uint i = 0; i < boots; i++) {
    sample = blocked_sample(block, nblk, nrow);
    bootstrap.slice(i) = fast_nei(f.rows(sample), jmt.rows(sample));
    t.it();
  }
  return(bootstrap);
}

arma::umat pairs(const uint n) {
  uint vi = 0;
  arma::umat pairs = arma::umat(n * (n-1)/2, 2);
  for (uint i = 0; i < n; i++) {
    for (uint j = i + 1; j < n; j++) {
      pairs(vi, 0) = i;
      pairs(vi, 1) = j;
      vi++;
    }
  }
  return(pairs);
}
// [[Rcpp::export(name=".bootstrap_nei")]]
arma::mat bootstrap_nei(const arma::mat& jxy, const arma::mat& jx, const uint block, const uint boots) {
  const uint nblks = std::ceil((double) jxy.n_rows / block);

  Tracker t(boots);
  arma::uvec sample;
  arma::umat comb = pairs(jx.n_cols);
  arma::mat bootstrap(boots, jxy.n_cols);
  for (int i = 0; i < boots; i++) {
    sample = blocked_sample(block, nblks, jxy.n_rows);
    bootstrap.row(i) =
      arma::log(arma::sum(jx(sample, comb.col(0)), 0))/2 +
      arma::log(arma::sum(jx(sample, comb.col(1)), 0))/2 -
      arma::log(arma::sum(jxy.rows(sample), 0));
    t.it();
  }
  return(bootstrap);
}
// Alternative implementation, but yield worst results
// // [[Rcpp::export(name = ".fnei")]]
// arma::mat fnei(const arma::mat& f, const arma::mat& non_na) {
//   arma::mat num, den;
//   num = arma::log(crossprod(f) + crossprod(1-f));
//   den = arma::log(arma::trans(1 + 2 * (arma::square(f)-f)) * non_na)/2;
//   return(den + den.t() - num);
// }