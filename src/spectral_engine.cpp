#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// Evaluate the SLD composite log-density for all observations.
//
// Z:           N x C indicator matrix (0/1).
// Z_mis:       N x C logical missingness mask.
// pi_c:        length-C vector of marginal probabilities.
// A_star:      C x C hollow projection matrix.
// item_starts: length-J vector of 0-based start indices for each item.
// item_sizes:  length-J vector of category counts per item.
//
// Returns an N-vector of per-observation composite log-densities.

// [[Rcpp::export]]
arma::vec eval_spectral_density_cpp(const arma::mat& Z,
                                    const arma::umat& Z_mis,
                                    const arma::vec& pi_c_raw,
                                    const arma::mat& A_star,
                                    const arma::ivec& item_starts,
                                    const arma::ivec& item_sizes) {
    int N = Z.n_rows;
    int C = Z.n_cols;
    int J = item_starts.n_elem;

    // Floor probabilities
    arma::vec pi_c = arma::clamp(pi_c_raw, 1e-15, 1.0);

    // Expectation matrix E_c (N x C, each row = pi_c)
    arma::mat E_c(N, C);
    for (int c = 0; c < C; c++) {
        E_c.col(c).fill(pi_c(c));
    }

    // Impute missing entries and compute residuals
    arma::mat Z_imp = Z;
    for (int i = 0; i < N; i++) {
        for (int c = 0; c < C; c++) {
            if (Z_mis(i, c)) {
                Z_imp(i, c) = E_c(i, c);
            }
        }
    }

    arma::mat R_c = Z_imp - E_c;

    // Shift: S = R * A_star
    arma::mat S_c = R_c * A_star;

    // eta = log(E_c) + S_c, clamped to [-30, 30]
    arma::mat eta = arma::log(E_c) + S_c;
    eta = arma::clamp(eta, -30.0, 30.0);

    // Per-item softmax and log-density accumulation
    arma::vec ld(N, arma::fill::zeros);

    for (int j = 0; j < J; j++) {
        int start = item_starts(j);
        int size  = item_sizes(j);

        for (int i = 0; i < N; i++) {
            // Find row max for this item's columns
            double row_max = eta(i, start);
            for (int c = 1; c < size; c++) {
                double val = eta(i, start + c);
                if (val > row_max) row_max = val;
            }

            // Softmax
            double denom = 0.0;
            for (int c = 0; c < size; c++) {
                denom += std::exp(eta(i, start + c) - row_max);
            }

            // Accumulate log P(observed category)
            for (int c = 0; c < size; c++) {
                if (Z(i, start + c) > 0.5) {
                    double p_hat = std::exp(eta(i, start + c) - row_max) / denom;
                    if (p_hat < 1e-15) p_hat = 1e-15;
                    ld(i) += std::log(p_hat);
                }
            }
        }
    }

    return ld;
}
