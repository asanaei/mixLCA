% Generated manually pending roxygenise(); regenerate with devtools::document().
\name{fit_lca}
\alias{fit_lca}
\title{Fit a Latent Class Model with mixLCA}
\usage{
fit_lca(
  data,
  continuous = NULL,
  categorical = NULL,
  concomitant = NULL,
  n_classes = 2L,
  dependence = "full",
  penalty = "auto",
  control = lca_control(),
  verbose = TRUE,
  cat_direct_effects = NULL,
  spectral_rank = 0L,
  spectral_pool = FALSE,
  init_model = NULL
)
}
\arguments{
\item{data}{Data frame of observations.}
\item{continuous}{Character vector of continuous indicator variable
names, or NULL.}
\item{categorical}{Character vector of categorical indicator variable
names, or NULL.}
\item{concomitant}{Character vector of concomitant predictor names,
or a one-sided formula (e.g.\ \code{~ age * income + poly(bmi, 2)}),
or NULL. Missing values in concomitant predictors are not permitted.}
\item{n_classes}{Integer >= 2.}
\item{dependence}{One of \code{"none"}, \code{"full"}, or \code{"penalized"}.}
\item{penalty}{Non-negative scalar or the string \code{"auto"}.}
\item{control}{List from \code{\link{lca_control}}.}
\item{verbose}{Logical.}
\item{cat_direct_effects}{List of two-element character vectors specifying
parent->child direct effects between categorical indicators.}
\item{spectral_rank}{Scalar or length-K integer vector of SLD ranks.}
\item{spectral_pool}{Logical: pool Burt matrices across classes.}
\item{init_model}{Optional \code{mixLCA} object for warm-starting.}
}
\value{
An opaque object of class \code{mixLCA}.
}
\description{
Estimates a finite mixture model with concomitant predictors,
mixed-type measurement indicators, and optional Spectral Local
Dependence (SLD). Parallel multi-start execution honours the
user's active \code{future::plan()}.
}
\seealso{
\code{\link{lca_control}}, \code{\link{distal}}, \code{\link{fit_indices}},
\code{\link{bvr_tests}}, \code{\link{enumerate_lca}}.
}
