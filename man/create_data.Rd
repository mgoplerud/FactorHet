% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_data.R
\name{create_data}
\alias{create_data}
\title{Create the (Sparse) Design Matrix for Analysis}
\usage{
create_data(
  design,
  penalty_for_regression,
  warn_missing = TRUE,
  verif_row = TRUE,
  remove_cols = Inf
)
}
\arguments{
\item{design}{Data frame with levels of each factor assigned to observation i}

\item{penalty_for_regression}{An object from "create_penalty"}
}
\description{
Create a sparse design matrix with weighted sum to zero constraints.
}
\keyword{internal}
