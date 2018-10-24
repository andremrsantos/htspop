#include "RcppArmadillo.h"

// Generate random block samples
arma::uvec blocked_sample(
    const uint bsize,
    const uint nblocks,
    const uint length);
// Compute the sum value within blocks
arma::vec blocked_sums(
    const arma::vec& val,
    const uint bsize,
    const uint nblock,
    const bool rm_na = true);
// Compute the mean value within blocks
arma::vec blocked_means(
    const arma::vec& val,
    const uint bsize,
    const uint nblock,
    const bool rm_na = true);
// Compute columns sum value within blocks
arma::mat blocked_col_sums(
    const arma::mat& matrix,
    const uint bsize,
    const uint nblock,
    const bool rm_na = true);
// Compute column mean value within blocks
arma::mat blocked_col_means(
    const arma::mat& matrix,
    const uint bsize,
    const uint nblock,
    const bool rm_na = true);