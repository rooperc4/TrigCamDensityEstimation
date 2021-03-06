% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DensityFunctions.R, R/DensityFunctions2.R
\name{logistic_model}
\alias{logistic_model}
\title{A logistic model function}
\usage{
logistic_model(species = "Fish", range, ob_dens, weights = 1,
  x0_limit = 2, near_limit = 1)

logistic_model(species = "Fish", range, ob_dens, weights = 1,
  x0_limit = 2, near_limit = 1)
}
\arguments{
\item{x0_limit}{lower limit of the x0 variable or the inflection point of the curve}

\item{near_limit}{closest range to fit the model or the closest range to a fish observable in both cameras}

\item{start_L}{Maximum length}

\item{x0}{initial size}

\item{start_K}{carrying capacity}

\item{start_L}{Maximum length}

\item{x0}{initial size}

\item{start_K}{carrying capacity}

\item{x0_limit}{lower limit of the x0 variable or the inflection point of the curve}

\item{near_limit}{closest range to fit the model or the closest range to a fish observable in both cameras}
}
\description{
This function allows you to fit a logistic function.

This function allows you to fit a logistic function.
}
\examples{
logistic_model()
logistic_model()
}
\keyword{logisitic}
\keyword{model}
