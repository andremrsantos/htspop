#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export(name=".fst_reich")]]
arma::mat fst_reich(
    const arma::mat& fq,
    const arma::mat& ht,
    const arma::mat& n) {
  arma::mat fst = arma::mat(fq.n_rows, 2);
  fst.col(0) = arma::square(fq.col(0)-fq.col(1)) - arma::sum(ht/n, 1);
  fst.col(1) = fst.col(0) + arma::sum(ht, 1);
  return(fst);
}

// [[Rcpp::export(name=".pairwise_fst_reich")]]
List pairwise_fst_reich(
    const arma::mat& count,
    const arma::mat& n) {
  const int npop = count.n_cols;
  const int nsnp = count.n_rows;
  const int npair = npop * (npop - 1) /2;
  const arma::mat fq = count / n;
  const arma::mat ht = count % (n - count) / (n % (n - 1));

  uint vi = 0;
  arma::uvec idx;
  arma::mat aux;
  arma::mat num = arma::mat(nsnp, npair);
  arma::mat den = arma::mat(nsnp, npair);
  for (uint i = 0; i < npop; i++) {
    for (uint j = i + 1; j < npop; j++) {
      idx = arma::uvec({i, j});
      aux = fst_reich(fq.cols(idx), ht.cols(idx), n.cols(idx));
      num.col(vi) = aux.col(0);
      den.col(vi) = aux.col(1);
      vi++;
    }
  }
  
  return List::create(
    Rcpp::Named("num") = num, 
    Rcpp::Named("den") = den);
}

// [[Rcpp::export(name = ".fst_weir_cockerham")]]
arma::vec fst_weir_cockerham(
  const arma::mat& fq,
  const arma::mat& n) {
  arma::vec nt, nc, fb, var, t1, t2;
  nt = arma::sum(n, 1);
  nc = (nt - arma::sum(arma::square(n), 1)/nt);
  fb = arma::sum(fq % n, 1) / nt;
  var = 2 * arma::sum(n % arma::square(fq.each_col() - fb), 1) / nt;
  t1 = var - (fb % (1-fb) - var/2)/(nt - 1);
  t2 =
    (2 * nc-1) / (nt - 1) % (fb % (1 - fb)) +
    (1 + (nt - 2 * nc)/(nt - 1)) % (var/2);
  return (t1/t2);
}

// [[Rcpp::export(name = ".pairwise_fst_weir_cockerham")]]
arma::mat pairwise_fst_weir_cockerham(
    const arma::mat& count,
    const arma::mat& n) {
  const int npop = count.n_cols;
  const int nsnp = count.n_rows;
  const int npair = npop * (npop - 1) /2;
  const arma::mat fq = count / n;
  
  uint vi = 0;
  arma::uvec idx;
  arma::mat fst = arma::mat(nsnp, npair);
  for (uint i = 0; i < npop; i++) {
    for (uint j = i + 1; j < npop; j++) {
      idx = arma::uvec({i, j});
      fst.col(vi) = fst_weir_cockerham(fq.cols(idx), n.cols(idx));
      vi++;
    }
  }
  return(fst);
}

// [[Rcpp::export(name=".fst_hudson")]]
arma::vec fst_hudson(
    const arma::mat& fq,
    const arma::mat& n) {
  arma::vec num = arma::sum(n / (n-1) % fq % (1-fq), 1);
  arma::vec den = fq.col(0) % (1-fq.col(1)) + fq.col(1) % (1-fq.col(0));
  return(1-num/den);
}

// [[Rcpp::export(name=".pairwise_fst_hudson")]]
arma::mat pairwise_fst_hudson(
    const arma::mat& count,
    const arma::mat& n) {
  const int npop = count.n_cols;
  const int nsnp = count.n_rows;
  const int npair = npop * (npop - 1) / 2;
  const arma::mat fq = count / n;

  uint vi = 0;
  arma::uvec idx;
  arma::mat fst = arma::mat(nsnp, npair);
  for (uint i = 0; i < npop; i++) {
    for (uint j = i + 1; j < npop; j++) {
      idx = arma::uvec({i, j});
      fst.col(vi) = fst_hudson(fq.cols(idx), n.cols(idx));
      vi++;
    }
  }
  return fst;
}

// [[Rcpp::export(name=".fst_wright")]]
arma::vec fst_wright(const arma::mat& fq) {
  return(arma::square(fq.col(0) - fq.col(1))/
         (arma::sum(fq, 1)) % (2 - arma::sum(fq, 1)));
}

// [[Rcpp::export(name=".pairwise_fst_wright")]]
arma::mat pairwise_fst_wright(
    const arma::mat& count,
    const arma::mat& total) {
  const int npop = count.n_cols;
  const int nsnp = count.n_rows;
  const int npair = npop * (npop - 1) / 2;
  const arma::mat fq = count/total;

  uint vi = 0;
  arma::uvec idx;
  arma::mat fst = arma::mat(nsnp, npair);
  for (uint i = 0; i < npop; i++) {
    for (uint j = i + 1; j < npop; j++) {
      idx = arma::uvec({i, j}); 
      fst.col(vi) = fst_wright(fq.cols(idx));
      vi++;
    }
  }
  return(fst);
}