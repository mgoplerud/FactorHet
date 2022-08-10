#' Small Dataset on Immigration Preferences
#'
#' An example dataset of 100 randomly chosen respondents from the replication
#' data in Hainmuller and Hopkins (2015). Only a small selection of the factors
#' and moderators in the original experiment are included in this example
#' dataset. The full data can be downloaded from the replication archive in
#' "Source" below. More details on all variables can be found in the original
#' paper.
#'
#' The replication data for Goplerud et al. (2022) provides code to process and
#' analyze the original data using \code{FactorHet}.
#'
#' @format A data frame with 1000 rows and 10 variables:
#' \describe{
#'   \item{CaseID}{Unique identifier for respondent.}
#'   \item{contest_no}{Task number (1-5) for each respondent}
#'   \item{choice_id}{Identifier for the profile shown, i.e. was it the "left"
#'   or "right" profile.}
#'   \item{Chosen_Immigrant}{Immigrant profile chosen by respondent.}
#'   \item{Country}{Immigrant's country of origin}
#'   \item{Ed}{Immigrant's education level}
#'   \item{Plans}{Immigrant's employment plans after arrival}
#'   \item{Gender}{Immigrant's gender}
#'   \item{party_ID}{\bold{Respondent's} party identification.}
#'   \item{census_div}{Level of immigration in \bold{respondent's} ZIP code}
#' }
#' @references
#' 
#' Goplerud, Max, Kosuke Imai, and Nicole E. Pashley. 2022. "Estimating
#' Heterogeneous Causal Effects of High-Dimensional Treatments: Application to
#' Conjoint Analysis." arxiv preprint: \url{https://arxiv.org/abs/2201.01357}
#' 
#' Hainmueller, Jens and Daniel J. Hopkins. 2015. "The Hidden American
#' Immigration Consensus: A Conjoint Analysis of Attitudes Toward Immigrants."
#' \emph{American Journal of Political Science} 59(3):529-548.
#' @source \url{https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/25505}
"immigration"
