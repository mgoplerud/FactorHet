% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict_functions.R
\name{predict.FactorHet}
\alias{predict.FactorHet}
\title{Predict after using FactorHet}
\usage{
\method{predict}{FactorHet}(
  object,
  newdata = NULL,
  type = "posterior",
  by_cluster = FALSE,
  return = "prediction",
  ...
)
}
\arguments{
\item{object}{A model estimated using \code{FactorHet} or
\code{FactorHet_mbo}.}

\item{newdata}{Data to estimate predictions on; default \code{NULL} uses
estimation data.}

\item{type}{How to generate the prediction? If \code{"posterior"} (default),
use the posterior probabilities for each observation in the training data
for weighting the clusters. If \code{"posterior_predictive"}, use the
cluster memberships implied by the moderator parameters, e.g. pi_{ik}. If
using this for evaluating out-of-sample accuracy, using
\code{"posterior_predictive"} may be sensible.}

\item{by_cluster}{Default of \code{FALSE}. If \code{TRUE}, return predicted
probabilities by cluster.}

\item{return}{Default of \code{"prediction"} returns the prediction. Use
\code{"detailed"} to return a variety of additional terms. This is mostly
used internally for estimating marginal effects and other quantities of
interest.}

\item{...}{Miscellanous options used internally and not documented.}
}
\description{
Get predictions of chosing particular profiles for either the data used
during estimation or some out-of-sample dataset.
}
\examples{
data(immigration)
# Fit a model once for simplicity
fit_MBO <- FactorHet(Chosen_Immigrant ~ Plans + Ed + Country,
 design = immigration, lambda = 1e-4,
 initialize = FactorHet_init(nrep = 1),
 K = 2, group = ~ CaseID, task = ~ contest_no, 
 choice_order = ~ choice_id)
immigration$pred_FH <- predict(fit_MBO)  
}
