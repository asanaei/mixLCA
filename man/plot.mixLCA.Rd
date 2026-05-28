% Generated manually pending roxygenise(); regenerate with devtools::document().
\name{plot.mixLCA}
\alias{plot.mixLCA}
\title{Plot Method for mixLCA}
\usage{
\method{plot}{mixLCA}(x,
                     type = c("profiles", "bvr", "distal",
                              "uncertainty", "convergence",
                              "categorical", "spectral_scree",
                              "spectral_loadings"),
                     ...)
}
\arguments{
\item{x}{A fitted \code{mixLCA} object.}
\item{type}{One of \code{"profiles"}, \code{"bvr"}, \code{"distal"},
\code{"uncertainty"}, \code{"convergence"}, \code{"categorical"},
\code{"spectral_scree"}, or \code{"spectral_loadings"}.}
\item{...}{Forwarded to the underlying plotting routine
(e.g.\ \code{data} for \code{"bvr"} and \code{"distal"},
\code{variable} for \code{"distal"}, \code{ci} for \code{"profiles"},
\code{dimension} for \code{"spectral_loadings"}).}
}
\value{
A \code{ggplot} object.
}
\description{
Single dispatch surface for all built-in visualisations.
}
