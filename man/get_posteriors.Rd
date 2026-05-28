% Generated manually pending roxygenise(); regenerate with devtools::document().
\name{get_posteriors}
\alias{get_posteriors}
\alias{get_posteriors.mixLCA}
\title{Extract Posterior Class Probabilities}
\usage{
get_posteriors(x, ...)

\method{get_posteriors}{mixLCA}(x, ...)
}
\arguments{
\item{x}{A fitted \code{mixLCA} object.}
\item{...}{Unused.}
}
\value{
Numeric matrix of posterior class probabilities (N x K).
}
\description{
Accessor returning the posterior probability matrix. Use this rather
than \code{model$posteriors} so downstream code is insulated against
future restructuring of internal fields.
}
