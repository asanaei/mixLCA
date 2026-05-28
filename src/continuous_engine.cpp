#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::vec eval_continuous_density_cpp(const arma::mat& Y, const arma::vec& mu, const arma::mat& Sigma) {
    int N = Y.n_rows;
    arma::vec ld(N, arma::fill::zeros);
    double log2pi = std::log(2.0 * arma::datum::pi);

    for(int i = 0; i < N; i++) {
        arma::vec yi = Y.row(i).t();
        arma::uvec obs = arma::find_finite(yi);
        int p = obs.n_elem;

        if (p == 0) continue;

        arma::vec y_obs = yi.elem(obs);
        arma::vec mu_obs = mu.elem(obs);
        arma::mat S_obs = Sigma.submat(obs, obs);
        
        S_obs = (S_obs + S_obs.t()) / 2.0;

        arma::vec eigval;
        arma::mat eigvec;
        arma::eig_sym(eigval, eigvec, S_obs);

        if(eigval.min() <= 0.0) {
            S_obs.diag() += 1e-5;
            arma::eig_sym(eigval, eigvec, S_obs);
        }

        double log_det = arma::sum(arma::log(eigval));
        arma::mat S_inv = eigvec * arma::diagmat(1.0 / eigval) * eigvec.t();

        arma::vec d_vec = y_obs - mu_obs;
        double quad = arma::as_scalar(d_vec.t() * S_inv * d_vec);

        ld(i) = -0.5 * (p * log2pi + log_det + quad);
    }
    return ld;
}
