% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mi.R
\name{MetaMi}
\alias{MetaMi}
\title{Generate imputations for a single covariate}
\usage{
MetaMi(t, x, v, m = 5, variable = NULL, type, fixed = TRUE, to_impute = NULL)
}
\arguments{
\item{t}{a tibble of effect sizes (one column called 't')}

\item{x}{a tibble of covariates}

\item{v}{a tibble of estimation variances (one column called 'v')}

\item{m}{number of imputations to generate}

\item{variable}{a string indicating which variable in \code(x) to impute}

\item{type}{string indicating covariate type ('binary', 'continuous')}

\item{fixed}{boolean indicates fixed/random effects model}

\item{to_impute}{optional argument, list of row indices to impute}
}
\value{
mids object of imputed datasets
}
\description{
Generate imputations for a single covariate
}
