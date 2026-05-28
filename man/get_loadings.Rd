% Generated manually pending roxygenise(); regenerate with devtools::document().
\name{get_loadings}
\alias{get_loadings}
\alias{get_loadings.mixLCA}
\title{Extract Spectral Loadings}
\usage{
get_loadings(x, ...)

\method{get_loadings}{mixLCA}(x, ...)
}
\arguments{
\item{x}{A fitted \code{mixLCA} object.}
\item{...}{Unused.}
}
\value{
A data frame from \code{\link{spectral_loadings}}, or \code{NULL}
if the model has no Spectral Local Dependence component.
}
\description{
Accessor returning the SLD loadings table. Use this rather than
\code{model$cat_spectral_params} so downstream code is insulated
against future restructuring of internal fields.
}
