#' Estimate heterogeneous effects in factorial and conjoint experiments
#' 
#' Fit a model to estimate heterogeneous effects in factorial or conjoint
#' experiments using a "mixture of experts" (i.e. a finite mixture of
#' regularized regressions with covariates affecting group assignment). Effects
#' are regularized using an overlapping group LASSO. \code{FactorHet_mbo} finds
#' an optimal lambda via Bayesian optimization whereas \code{FactorHet} requires
#' a lambda to be provided. \code{FactorHet_mbo} typically used in practice.
#' 
#' @param formula Formula specifying model. The syntax is \code{y ~ X1 + X2}
#'   where \code{y} is the outcome and \code{X1} and \code{X2} are factors.
#'   Interactions can be specified using \code{*} syntax. All main factors must
#'   be explicitly included.
#' @param design A data.frame containing the data to be analyzed.
#' @param K An integer specifying the number of groups; \code{K=1} specifies a
#'   model with a single group.
#' @param lambda A positive numeric value denoting regularization strength; this
#'   is scaled internally by the number of observations, see
#'   \code{\link{FactorHet_control}}. \code{FactorHet_mbo} calibrates through
#'   model-based optimization. "Details" provides more discussion of this
#'   approach.
#' @param moderator A formula of variables (moderators) that affect the prior
#'   probability of group membership. This is ignored when \code{K=1} or
#'   \code{moderator=NULL}.
#' @param group A formula of a single variable, e.g. \code{~ person_id}, that is
#'   used when there are repeated observations per individual.
#' @param task A formula of a single variable that indicates the task number
#'   performed by each individual. This is not used when \code{group} is
#'   unspecified.
#' @param choice_order A formula of a single variable that indicates which
#'   profile is on the "left" or "right" in a conjoint experiment.
#' @param weights A formula of a single variable that indicates the weights for
#'   each observation (e.g., survey weights). If \code{group} is specified, the
#'   weights must be constant inside of each value of group.
#' @param control An object from \code{\link{FactorHet_control}} that sets
#'   various model estimation options.
#' @param initialize An object from \code{\link{FactorHet_init}} that determines
#'   how the model is initialized.
#' @param verbose A logical value that prints intermediate information about
#'   model fitting. The default is \code{TRUE}.
#'   
#' @details 
#'   
#'   \bold{Caution:} Many settings in \link{FactorHet_control} can be modified
#'   to allow for slight variations in how the model is estimated. Some of these
#'   are faster but may introduce numerical differences across versions of
#'   \code{R} and machines. The default settings aim to mitigate this. One of
#'   the default settings (\code{FactorHet_control(step_SQUAREM=NULL)})
#'   considerably increases the speed of convergence and the quality of the
#'   optimum located at the expense of sometimes introducing numerical
#'   differences across machines. To address this, one could not use SQUAREM
#'   (\code{do_SQUAREM=FALSE}) or set it to use some fixed step-size (e.g.,
#'   \code{step_SQUAREM=-10}). If SQUAREM produces a large step, a message to
#'   this effect will be issued.
#'   
#'   \bold{Factorial vs. Conjoint Experiment:} A factorial experiment, i.e.
#'   without a forced-choice between profiles, can be modeled by ignoring the
#'   \code{choice_order} argument and ensuring that each \code{group} and
#'   \code{task} combination corresponds to exactly one observation in the
#'   design.
#'   
#'   \bold{Estimation:} All models are estimated using an AECM algorithm
#'   described in Goplerud et al. (2025). Calibration of the amount of
#'   regularization (i.e. choosing \eqn{\lambda}), should be done using
#'   \code{FactorHet_mbo}. This uses a small number (default 15) of attempts to
#'   calibrate the amount of regularization by minimizing a user-specific
#'   criterion (defaulting to the BIC), and then fits a final model using the
#'   \eqn{\lambda} that is predicted to minimize the criterion.
#'   
#'   Options for the model based optimization (\code{mbo}) can be set using
#'   \code{\link{FactorHet_mbo_control}}. Options for model estimation can be
#'   set using \code{\link{FactorHet_control}}.
#'   
#'   \bold{Ridge Regression:} While more experimental, ridge regression can be
#'   estimated by setting \code{lambda = 0} (in \code{FactorHet}) and then
#'   setting \code{prior_var_beta} in \code{\link{FactorHet_control}} or by using
#'   \code{FactorHet_mbo} and setting \code{mbo_type = "ridge"}.
#'   
#'   \bold{Moderators:} Moderators can be provided via the \code{moderator}
#'   argument. These are important when \code{K > 1} for ensuring the stability
#'   of the model. Repeated observations per individual can be specified by
#'   \code{group} and/or \code{task} if relevant for a force-choice conjoint.
#'   
#' @examples
#' # Use a small subset of the immigration data from Hainmueller and Hopkins
#' data(immigration)
#' \donttest{
#' set.seed(1)
#' # Fit with two groups and tune regularization via MBO
#' fit_MBO <- FactorHet_mbo(
#'   formula = Chosen_Immigrant ~ Country + Ed + Gender + Plans,
#'   design = immigration, group = ~ CaseID,
#'   task =  ~ contest_no, choice_order = ~ choice_id,
#'   # Only do one guess after initialization for speed
#'   mbo_control = FactorHet_mbo_control(iters = 1),
#'   K = 2)
#' # Plot the raw coefficients
#' cjoint_plot(fit_MBO)
#' # Check how MBO fared at calibrating regularization
#' visualize_MBO(fit_MBO)
#' # Visualize posterior distribution of group membership
#' posterior_FactorHet(fit_MBO)
#' # Get AMEs
#' AME(fit_MBO)
#' }
#' @return Returns an object of class \code{FactorHet}. Typical use will involve
#'   examining the patterns of estimated treatment effects.
#'   \code{\link{cjoint_plot}} shows the raw (logistic) coefficients.
#'   
#'   Marginal effects of treatments (e.g. average marginal effects) can be
#'   computed using \code{\link{AME}}, \code{\link{ACE}}, or \code{\link{AMIE}}.
#'   
#'   The impact of moderators on group membership can be examined using
#'   \code{\link{margeff_moderators}} or \code{\link{posterior_by_moderators}}.
#'   
#'   The returned object is a list containing the following elements:
#'   \describe{
#'   \item{parameters: }{Estimated model parameters. These are usually obtained
#'   via \code{\link{coef.FactorHet}}.}
#'   \item{K: }{The number of groups}
#'   \item{posterior: }{Posterior group probability for each observation. This
#'   is list of two data.frames one with posterior probabilities
#'   (\code{"posterior"}) and one (\code{"posterior_predictive"}) implied solely
#'   by the moderators, i.e. \eqn{\pi_{k}(X_i)} from Goplerud et al. (2025).}
#'   \item{information_criterion: }{Information on the BIC, degrees of freedom,
#'   log-likelihood, and number of iterations.}
#'   \item{internal_parameters: }{A list of many internal parameters. This is
#'   used for debugging or by other post-estimation functions.}
#'   \item{vcov: }{Named list containing the estimated variance-covariance
#'   matrix. This is usually extracted with \code{vcov}.}
#'   \item{lp_shortEM: }{If \code{"short EM"} is applied (only applicable if 
#'   \code{FactorHet}, not \code{FactorHet_mbo}, is used), it lists the
#'   log-posterior at the end of each short run.}
#'   \item{MBO: }{If \code{FactorHet_mbo} is used, information about the
#'   model-based optimization (MBO) is stored here. \code{\link{visualize_MBO}} 
#'   provides a quick graphical summary of the BIC at different \eqn{\lambda}.}
#'   }
#' @importFrom graphics plot lines
#' @export
FactorHet <- function(formula, design, K, lambda, 
  moderator = NULL, group = NULL, task = NULL, choice_order = NULL, 
  weights = NULL, control = FactorHet_control(), 
  initialize = FactorHet_init(), verbose = TRUE){
  
  if (!inherits(initialize, 'FactorHet_init')){
    stop('initialize must be made using FactorHet_init()')
  }
  if (initialize$verbose){
    verbose_function <- function(x){x}
  }else{
    verbose_function <- suppressMessages
  }
  
  if (lambda == 0 & control$log_method != 'standard'){
    control$log_method <- 'standard'
    verbose_function(message('Ridge regression (lambda = 0) requires standard projection.'))
  }
  controlFH_init <- control
  
  # Allow "short_EM" to be provided directly to
  # FactorHet_control as short-hand.
  if (!inherits(control$init_method, 'list')){
    if (control$init_method == 'short_EM'){
      initialize$short_EM <- TRUE
    }else if (control$init_method %in% c('mclust', 'mclust_aug', 'spectral', 'kmeans')){
      # #If using deterministic initalization, ONLY run the model once
      initialize$short_EM <- FALSE
      initialize$nrep <- 1
    }
  }else{
    initialize$short_EM <- FALSE
    initialize$nrep <- 1
  }
  
  short_EM <- initialize$short_EM
  nrep <- initialize$nrep
  do_SE <- control$calc_se
  
  if (K == 1 & initialize$force_rep == FALSE){
    if (short_EM){
      controlFH_init$init_method <- initialize$short_EM_init
      control$init_method <- initialize$short_EM_init
      short_EM <- FALSE
    }
    if (nrep > 1 & initialize$verbose){
      message('For K = 1, no need to repeatedly initialize. If really necessary set "force_rep = TRUE". Doing only one run')
    }
    nrep <- 1
  }
  mlist <- as.list(rep(NA, nrep))
  
  if (short_EM){
    controlFH_init$iterations <- initialize$short_EM_it
    controlFH_init$calc_df <- FALSE
    controlFH_init$init_method <- initialize$short_EM_init
    controlFH_init$beta_cg_it <- initialize$short_EM_cg_it
    controlFH_init$beta_method <- initialize$short_EM_beta_method
    controlFH_init$skip_check_rank <- TRUE
    controlFH_init$maxit_pi <- initialize$short_EM_pi
  }else{
    controlFH_init <- control
  }
  
  all_init_options <- list(
    formula = formula, design = design, lambda = lambda, K = K,
    moderator = moderator, weights = weights,
    group = group, task = task, choice_order = choice_order, control = controlFH_init
  )
  #If multiple repetitions, turn off standard errors for most runs 
  if (nrep > 1 & short_EM){
    all_init_options$control$calc_se <- FALSE
  }
  
  if (short_EM){message('Doing Short EM')}
  
  plot_func <- function(x){log(diff(x))/log(10)}
  
  for (i in 1:nrep){
    if (nrep > 1){message('|', appendLF = FALSE)}
    time_short <- proc.time()
    mlist[[i]] <- verbose_function(do.call('EM_analysis', args = all_init_options))
    time_short <- time_short - proc.time()
    
    if (mlist[[i]]$internal_parameters$diagnostic$basic){
      message(
        'FactorHet_control() has non-default options that may induce numerical differences across machines.'
      )
    }
    if (mlist[[i]]$internal_parameters$diagnostic$SQUAREM){
      msg <- c(
        'SQUAREM has proposed a large step size (before backtracking); this often results in',
        'good performance but may induce numerical differences across machines.'
      )
      message(paste(msg, collapse='\n'))
    }
    
    if (initialize$plot_repeat){
      if (i == 1){
        plot(log(diff(mlist[[i]]$internal_parameters$trajectory$ll[,1])), type = 'l', col = i)
      }else{
        lines(log(diff(mlist[[i]]$internal_parameters$trajectory$ll[,1])), col = i)
      }
    }
  }
  rm(all_init_options)
  
  ll.final <- sapply(mlist, logLik)

  if (nrep > 1 | short_EM){
    message('\n', appendLF = FALSE)
  }  
  
  if (short_EM){
    if (initialize$debug_repeat){
      print(cbind(ll.final, 
                  sapply(mlist, FUN=function(i){i$parameters$pi[1]}), 
                  sapply(mlist, FUN=function(i){max(abs(coef(i)))})))
    }
    
    best_short <- mlist[[which.max(ll.final)]]
    
    all_final_options <- list(K = K, formula = formula, design = design, 
                              lambda = lambda, moderator = moderator, weights = weights,
                              group = group, task = task, choice_order = choice_order, control = control
    )
    
    all_final_options$control$init_method <- list(beta = best_short$parameters$nullspace_beta, 
                                                  pi = best_short$parameters$pi, 
                                                  phi = best_short$parameters$phi, 
                                                  group_E.prob = as.matrix(best_short$posterior$posterior[,-1]))
    all_final_options$control$calc_se <- do_SE
    
    message('Starting Long EM')
    mlist <- suppressMessages(do.call('EM_analysis', args = all_final_options))
    mlist$lp_shortEM <- ll.final
    return(mlist)
  }else{
    if (initialize$debug_repeat){
      print(cbind(ll.final, sapply(mlist, FUN=function(i){i$parameters$pi[1]}), 
                  sapply(mlist, FUN=function(i){max(abs(coef(i)))})))
    }
    if (initialize$return_all){
      return(mlist)
    }else{
      mlist <- mlist[[which.max(ll.final)]]
      mlist$lp_shortEM <- ll.final
      return(mlist)
    }
  }
}
