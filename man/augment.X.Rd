% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper_functions.R
\name{augment.X}
\alias{augment.X}
\title{Augment X for missing data approach for MNAR}
\usage{
augment.X(X)
}
\arguments{
\item{X}{matrix of covariates, including some missing values (NAs)}
}
\value{
Matrix same size as X, with augmented columns and zeros in the missing spots
}
\description{
Augment X for missing data approach for MNAR
}
\examples{
set.seed(1)
n <- 100
X <- matrix(runif(n*2, 0, 1), ncol=2)
X[sample(1:length(X), round(.1*length(X)))] <- NA
X.new <- augment.X(X)
sum(is.na(X.new))
}
