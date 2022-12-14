% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/postregression_functions.R
\name{marginal_AME}
\alias{marginal_AME}
\alias{marginal_ACE}
\alias{marginal_AMIE}
\title{Calculate Marginal Effects}
\usage{
marginal_AME(
  object,
  baseline = NULL,
  vcov = TRUE,
  design = NULL,
  ignore_restrictions = FALSE,
  vcov.type = NULL,
  average_position = TRUE,
  verbose = TRUE,
  plot = TRUE,
  extra_restriction = NULL
)

marginal_ACE(
  object,
  baseline,
  design = NULL,
  average_position = TRUE,
  ignore_restrictions = FALSE,
  extra_restriction = NULL,
  verbose = TRUE,
  plot = TRUE
)

marginal_AMIE(
  object,
  design = NULL,
  baseline = NULL,
  average_position = TRUE,
  ignore_restrictions = FALSE,
  verbose = FALSE,
  plot = TRUE
)
}
\arguments{
\item{object}{An object from FactorHet or FactorHet_mbo.}

\item{baseline}{For \code{marginal_AME}, the options are \code{NULL} (choose
baseline as first factor level), a named list for each factor and a
corresponding baseline level, or \code{NA} for no baseline (approximating
the "marginal means" from Leeper et al. [2020]). For \code{marginal_ACE}, a
named list of two factors with a baseline level must be provided. For
\code{marginal_AMIE}, you can provide either a named list of two factors
with a baseline level or \code{NULL} to estimate all pairwise AMIE.}

\item{vcov}{Calculate standard errors for the marginal effect. Not yet implemented for the AMIE.}

\item{design}{The dataset used to estimate the marginal effects. By default
(NULL), the estimation data is used.}

\item{ignore_restrictions}{Ignore any randomization restrictions when
calculating the marginal effects. See "Details" for more information.}

\item{vcov.type}{NULL uses the standard errors associated with estimation.
Pass arguments from \code{vcov.FactorHet}'s "se.method" argument.}

\item{average_position}{For factorial designs, should the marginal effect for the "left"
profile and the "right" profile be averaged? See Goplerud et al. (2022) for
discussion.}

\item{verbose}{"TRUE" returns information on the progress of estimating the marginal effects.}

\item{plot}{"FALSE" only returns the data of the estimated marginal effects.}

\item{extra_restriction}{Manually add randomization restrictions for
calculating the marginal effects. See "Details" for more information.}
}
\value{
\item{plot}{The ggplot2 visualization of the marginal effects.}
\item{data}{The data used to create the ggplot2 visualization.}
}
\description{
Calculate the average marginal [component] effect (AME or AMCE), the average
combination effect (ACE), or the average marginal interaction effect (AMIE)
with a FactorHet model.
}
\details{
\bold{Choice of Baseline}: For ACE and AMIE, a choice of baseline is
required. This should be provided as a named list with two elements for two
factors where each element of the list contains the baseline level. See the
examples for usage and Egami and Imai (2019) for details on the ACE and AMIE.

For AME, a choice of baseline corresponds to a "standard" AME (see Egami and
Imai 2019). The option \code{NULL} choses the first level of each factor. It
can be manually specified using a named list. If a named list is provided,
only AMEs for those named factors are calculated. This can be helpful if
there are many factors.

If \code{NA} is provided as the baseline level, the AME is calculated without
a baseline; while this does not correspond to a "proper" AME, it is designed
to approximate the "marginal means" discussed in Leeper et al. (2020). Note
that in the presence of randomization restrictions, the quantity estimated
with a \code{NA} baseline may not be centered around 0.5. Ignoring the
randomization restrictions may be useful in this scenario. Goplerud et al.
(2022)'s supporting information provides more discussion of this point.

\bold{Randomization Restrictions}: Randomization restrictions can be set in
one of two ways. By default, FactorHet checks whether for each
pairwise combinations of factors, some combination of levels do not occur at all
(e.g. "doctor" and "high school") or whether some included interactions are
extremely rare (see \code{rare_threshold} in \code{FactorHet_control}). Those are assumed
to be the randomization restrictions implied by the design as, for designs of
moderate size and uniform randomization, the probability of any combination
not occurring is very low. 

However, if this is inappropriate for a specific dataset, randomization
restrictions can be set automatically as follows. First, set
\code{ignore_restrictions = TRUE}. This will ignore all "data-driven"
estimates of randomization restrictions. Second, the argument
\code{extra_restriction} should be a named list where the name of each
element corresponds to a factor (e.g. "Job") and each element is a vector of
the levels that \emph{cannot} be used. When using this approach,
\code{marginal_AME} should be used only for one factor at a time. An example
is shown below.

\bold{Plots}: Note that for the ggplot2 visualizations of the ACE and AMIE,
gray squares indicate combinations that are excluded due to randomization
restrictions. White indicates baseline levels.
}
\examples{
data(immigration)
# Induce "fake" randomization restriction
immigration$joint_id <- paste(immigration$CaseID, immigration$contest_no)
remove_profiles <- subset(immigration, Plans == 'No plans' & Ed == 'GradDeg')
immigration <- subset(immigration, !(joint_id \%in\% remove_profiles$joint_id))
# Fit with one cluster and limited regularization for example only
fit_MBO <- FactorHet(Chosen_Immigrant ~ Plans + Ed + Country,
  design = immigration, lambda = 1e-4,
  K = 1, group = ~ CaseID, task = ~ contest_no, choice_order = ~ choice_id)
# Estimate AME of "promise" with randomization restriction
est_AME <- marginal_AME(fit_MBO, baseline = list('Ed' = 'GradDeg'))
# Estimate AME ignoring randomization restriction
est_AME_norr <- marginal_AME(fit_MBO, 
  baseline = list('Ed' = 'GradDeg'), ignore_restrictions = TRUE)
# Estimate AME by manually specifying randomization restrictions
est_AME_rr_manual <- marginal_AME(fit_MBO,
  baseline = list('Ed' = 'GradDeg'), ignore_restrictions = TRUE,
  extra_restriction = list('Plans' = 'No plans'))
stopifnot(isTRUE(all.equal(est_AME$data, est_AME_rr_manual$data)))
# Estimate without baseline
est_MM <- marginal_AME(fit_MBO, baseline = list('Ed' = NA))
# Estimate ACE and AMIE
est_ACE <- marginal_ACE(fit_MBO, baseline = list('Ed' = 'GradDeg', 'Plans' = 'Has contract'))
est_AMIE <- marginal_AMIE(fit_MBO, baseline = list('Ed' = 'GradDeg', 'Plans' = 'Has contract'))
}
\references{
Egami, Naoki and Kosuke Imai. 2019. "Causal Interaction in
  Factorial Experiments: Application to Conjoint Analysis." \emph{Journal of the
  American Statistical Association}. 114(526):529-540.

  Goplerud, Max, Kosuke Imai, and Nicole E. Pashley. 2022. "Estimating
  Heterogeneous Causal Effects of High-Dimensional Treatments: Application to
  Conjoint Analysis." arxiv preprint: \url{https://arxiv.org/abs/2201.01357}
 
  Leeper, Thomas J., Sara B. Hobolt, and James Tilley. 2020. "Measuring Subgroup
  Preferences in Conjoint Experiments." \emph{Political Analysis}. 28(2):207-221.
}
