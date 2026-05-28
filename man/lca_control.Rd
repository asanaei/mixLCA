% Generated manually pending roxygenise(); regenerate with devtools::document().
\name{lca_control}
\alias{lca_control}
\title{Control Parameters for mixLCA}
\usage{
lca_control(max_iter = 500L, tol = 1e-6, n_starts = 1L,
            seed = 110L, kmeans_nstart = 1L)
}
\arguments{
\item{max_iter}{Maximum EM iterations per start.}
\item{tol}{Convergence tolerance on absolute log-likelihood change.}
\item{n_starts}{Integer: number of random starting points.}
\item{seed}{Base random seed; start \emph{s} uses \code{seed + s}.
The global \code{.Random.seed} is never modified.}
\item{kmeans_nstart}{Integer: random initializations for the internal
\code{kmeans} used to seed continuous-indicator starting values.}
}
\value{
A list of control values to pass to \code{\link{fit_lca}}.
}
\description{
Bundles optimiser settings into a single list for use with
\code{\link{fit_lca}}.
}
