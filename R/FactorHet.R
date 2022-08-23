#' Estimating Heterogeneous Effects in Factorial and Conjoint Experiments
#' 
#' Fit a model to estimate heterogeneous effects in factorial or conjoint
#' experiments using a "mixture of experts" (i.e. a finite mixture of
#' regularized regressions with covariates affecting cluster assignment).
#' Effects are regularized using an overlapping group LASSO.
#'
#' All models are estimated using an AECM algorithm described in Goplerud et al.
#' (2022).
#' 
#' @param formula Formula specifying model. The syntax is y ~ X1 + X2 where y is
#'   the outcome and X1 and X2 are factors. Interactions can be specified using
#'   * syntax. All main factors must be explicitly included.
#' @param design A data.frame containing the experimental data.
#' @param K The number of clusters to estimate.
#' @param lambda The regularization strength; it is not common to set this
#'   manually versus using FactorHet_mbo to optimize this term.
#' @param moderator A formula of variables (moderators) that affect the prior
#'   probability of cluster membership. Ignored when K = 1.
#' @param group A formula of a single variable, e.g. \code{~ person_id}, that is
#'   used when there are repeated observations per individual.
#' @param task A formula of a single variable that indicates the task number
#'   performed by each individual. This is not used when \code{group} is unspecified.
#' @param choice_order A formula of a single variable that indicates which
#'   profile is on the "left" or "right" in a conjoint experiment.
#' @param weights A vector of weights for each observation (e.g. survey
#'   weights). If group is specified, they must be constant inside of each value
#'   of group.
#' @param control An object from \link{FactorHet_control} that sets various model
#'   estimation options.
#' @param initialize An object from \link{FactorHet_init} that determines how the
#'   model is initialized.
#' @param verbose A logical term that prints intermediate information about
#'   model fitting.
#'   
#' @details \bold{Estimation:} Calibration of the amount of regularization
#'   should be done using \code{FactorHet_mbo}. This uses some number of
#'   attempts to calibrate the amount of regularization by minimizing some
#'   criterion (e.g. the BIC), and then fits a final model using the optimal
#'   estimated regularization. Note that by default \code{FactorHet_mbo}
#'   requires the suggested package \code{tgp} to be installed. Different
#'   methods can be used; see the documentation for \link{FactorHet_mbo_control}
#'   for more information.
#'
#'   Options for the model based optimization (mbo) can be set using
#'   \code{FactorHet_mbo_control()}. Options for model estimation can be set
#'   using \code{FactorHet_control()}.
#'   
#'   \bold{Ridge Regression:} While more experimental, ridge regression can be
#'   estimated by setting \code{lambda = 0} (in \code{FactorHet}) and then
#'   setting \code{prior_var_beta} in \code{FactorHet_control} or by using
#'   \code{FactorHet_mbo} and setting \code{mbo_type = "ridge"}.
#'   
#'   \bold{Moderators:} Moderators can be provided via the \code{moderator}
#'   argument. These are important when \code{K > 1} for ensuring the stability
#'   of the model. Repeated observations per individual can be specified by
#'   \code{group} and/or \code{task} if relevant for a force-choiced conjoint.
#' 
#' @examples
#' # Use a small subset of the immigration data from Hainmueller and Hopkins
#' data(immigration)
#' # Fit with two clusters and tune regularization via MBO
#' # Only do one iteration (iters = 1) for speed
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
#' # Visualize posterior distribution of cluster membership
#' posterior_FactorHet(fit_MBO)
#' # Get AMEs
#' marginal_AME(fit_MBO)
#' @return Returns an object of class \code{FactorHet}. Typical use will involve
#'   examining the patterns of estimated treatment effects. \code{cjoint_plot} shows
#'   the raw (logistic) coefficients. 
#'   
#'   Marginal effects (e.g. average marginal
#'   effect) can be computed using \code{marginal_AME}, \code{marginal_ACE}, or
#'   \code{marginal_AMIE}. 
#'   
#'   The effects of moderators can be examined using
#'   \code{moderator_AME} or \code{posterior_by_moderators}. Please see the corresponding
#'   documentation for details.
#'   
#'   The returned object is a list containing the following elements:
#'   \itemize{
#'   \item{parameters: }{Estimated model parameters. Usually extract these via
#'   \code{coef}.}
#'   \item{K: }{The number of clusters}
#'   \item{posterior: }{Posterior cluster probability for each observation; list
#'   of two data.frames one with posterior (\code{"posterior"}) and one
#'   (\code{"posterior_predictive"}) based on the probabilities \eqn{\pi_{ik}} implied
#'   by the moderators.}
#'   \item{information_criterion: }{Information on the information criterion estiamted.}
#'   \item{internal_parameters: }{Internal parameters; for debugging or use by other functions}
#'   \item{vcov: }{Named list containing the estimated variance-covariance
#'   matrix, usually extracted with \code{vcov}. Other objects for internal
#'   use.}
#'   \item{lp_shortEM: }{If "short EM" is applied, the log-posterior at the end
#'   of each short run.}
#'   \item{MBO:}{If \code{FactorHet_mbo} is used, information about the
#'   model-based optimization is stored here. Use \link{visualize_MBO} to see a
#'   graphical summary.}
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
  if (class(control$init_method) != 'list'){
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
  
  message('\n', appendLF = FALSE)
  
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
