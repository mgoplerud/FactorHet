#' Control for FactorHet estimation
#' 
#' Provides a set of control arguments to \code{\link{FactorHet}}. Arguments
#' around the initialization of the model (important when \code{K > 1}) can be
#' set via \code{\link{FactorHet_init}} and arguments for the model-based
#' optimization tuning of regularization strength \eqn{\lambda} can be found in
#' \code{\link{FactorHet_mbo_control}}. The parameters can be divided into ones
#' governing the model priors, model estimation, and miscellaneous settings. All
#' arguments have default values.
#' 
#' @param single_intercept A logical value or \code{NULL} that indicates whether
#'   a single intercept should be used across groups. The default is \code{NULL}
#'   which uses a single intercept if the study is a forced-choice conjoint
#'   (i.e., \code{choice_order} is used) and a varying intercept by group
#'   otherwise.
#' @param gamma A non-negative numerical value that determines whether
#'   sparsity-inducing prior be "spread" across groups in proportion to the
#'   average prior probability of membership. Default of 1; see Städler et al.
#'   (2010) and Goplerud et al. (2025) for more discussion.
#' @param adaptive_weight An argument that determines the weights given to
#'   different terms in the penalty function. The default (\code{"B&R"}) uses
#'   Bondell and Reich (2009), generalized appropriately if needed, see Goplerud
#'   et al. (2025) for discussion. If a matrix is provided (e.g. from a prior
#'   run of \code{\link{FactorHet}}), this can be used to set up an "adaptive
#'   overlapping group LASSO". \code{"none"} imposes no weights. To use a matrix
#'   and \emph{not} use Bondell and Reich weights, additional set \code{override_BR =
#'   TRUE}.
#' @param log_method An argument for specifying whether latent overlapping
#'   groups should be used when interactions are included. The default is
#'   \code{"log_ginv"}. Options beginning with \code{"log_"} employ latent
#'   overlapping groups (see Yan and Bien 2017 and the supporting information of
#'   Goplerud et al. 2025). The projection matrix can be either the generalized
#'   inverse extending Post and Bondell (2013) (\code{"log_ginv"}), a random
#'   matrix (\code{"log_random"}), or zero (\code{"log_0"}). \code{"standard"}
#'   does not implement overlapping groups.
#' @param prior_var_phi A numerical value that encodes the variance of
#'   multivariate normal prior on moderator coefficients. \bold{Note:} The
#'   moderators are not standardized internally and thus should be on broadly
#'   comparable scales to avoid differential amounts of regularization on
#'   different moderators. The default value is 4.
#' @param prior_var_beta A numerical value of normal prior on each treatment
#'   effect coefficient. The default is \code{Inf} when using sparse estimation.
#'   A different value can be set when using "ridge" regression, i.e.
#'   \code{lambda=0}.
#' @param iterations A numerical value setting the maximum number of iterations used in
#'   the algorithm. The default is 1000.
#' @param tolerance.parameters A numerical value setting the one convergence
#'   criterion: When no parameter changes by more than this
#'   amount, terminate the algorithm. Default is 1e-5.
#' @param tolerance.logposterior A numerical value setting the one convergence
#'   criterion: When the log-posterior changes by less than this amount,
#'   terminate the algorithm. Default is 1e-5.
#' @param maxit_pi An argument setting the maximum number of iterations used in
#'   each M-Step that updates the moderators. The default is \code{NULL} and
#'   uses default settings in optimizer. For \code{"lib_lbfgs"}, this optimizes
#'   until convergence is obtained.
#' @param optim_phi_controls A list of options for optimizer used in updating
#'   the moderator
#'   parameters. A method must be provided at minimum, e.g., \code{list(method =
#'   "lib_lbfgs")}. \code{"lib_lbfgs"} uses \code{\link[lbfgs]{lbfgs}} from the
#'   accompanying package. All other options use the base \code{\link{optim}}
#'   function in \code{R}. The maximum number of iterations should be specified
#'   via \code{maxit_pi}. All other options are specified through this argument.
#' @param repeat_beta An integer setting the number of times to repeat the E-M
#'   cycle for updating \eqn{\beta} before moving to update the moderator
#'   parameters \eqn{\phi}. The default is 1.
#' @param init_method An argument for initializing the algorithm. One set of
#'   options are different character values: \code{"kmeans"} (k-means clustering
#'   on the moderators), \code{"mclust"} (\code{"mclust"} on the moderators),
#'   \code{"random_pi"} (random probabilities of group membership for each
#'   person), \code{"random_member"} (random hard assignment),
#'   \code{"random_beta"} (random coefficients). This can be set with a named
#'   list with group membership probabilities. This should consist of a named
#'   list with a single element \code{"group_E.prob"} that is a data.frame which
#'   contains probabilities for each group/unit with the column names
#'   \code{"group"} and then \code{"group_[0-9]+"} depending on \code{K}. In
#'   general, when using \code{\link{FactorHet_mbo}}, this argument is not used
#'   and rather set via the relevant options in
#'   \code{\link{FactorHet_mbo_control}} as this will ensure the same
#'   initialization for all runs of \code{FactorHet_mbo}.
#' @param return_data A logical value for whether the formatted data should be
#'   returned. The default is \code{FALSE}.
#' @param rare_threshold A numerical value setting the threshold for which
#'   interactions should be excluded. If an interaction of two factors has fewer
#'   than \code{rare_threshold} observations, the corresponding interaction term
#'   will not be included. This is a way to enforce randomization restrictions.
#'   The default is \code{5} but setting it to 0 will ensure that all
#'   interactions are included. The documentation of \code{\link{FactorHet}}
#'   provides more discussion.
#' @param rare_verbose A logical value as to whether to print information about
#'   the rare interactions. The default is \code{TRUE}.
#' @param beta_method A character value for the method by which \eqn{\beta} is
#'   updated. The default is \code{"cpp"}. An alternative that uses conjugate
#'   gradient (\code{"cg"}) is faster per-iteration but may introduce numerical
#'   differences across platforms.
#' @param beta_cg_it A numerical value of the number of conjugate gradient steps
#'   to use if \code{beta_method = "cg"}.
#' @param lambda_scale A function for internally rescaling lambda to be a
#'   function of \eqn{N}. Options are \code{"N"} (default; \code{lambda * N}),
#'   \code{"unity"} (i.e. no rescaling), or \code{"root_N"} (\code{lambda *
#'   sqrt(N)}).
#' @param weight_dlist A logical value for whether to weight additional
#'   penalties following Hastie and Lim (2015). The default is \code{FALSE}.
#' @param do_SQUAREM A logical value for whether to perform SQUAREM to
#'   accelerate convergence. The default is \code{TRUE}.
#' @param step_SQUAREM An argument specifying the step size to use for SQUAREM.
#'   The default is \code{NULL} which uses a data-driven step size. This
#'   generally performs well, but may introduce numerical differences across
#'   machines. See the documentation of \code{\link{FactorHet}} for more
#'   discussion.
#' @param backtrack_SQUAREM An integer that sets the number of backtracking
#'   steps to perform for SQUAREM. The default is 10.
#' @param df_method A character value specifying the method calculating degrees
#'   of freedom. Default of \code{"EM"} follows Goplerud et al. (2025) and
#'   calculates the degrees of freedom using the Polya-Gamma weights.
#'   \code{"IRLS"} uses \eqn{\zeta_{ik} (1 -
#'   \zeta_{ik})} as weights, 
#'   where \eqn{\zeta_{ik} = Pr(y_i = 1 | X_i, z_i = k)}.
#'   \code{"free_param"} counts the number of parameters after fusion and
#'   accounting for the sum-to-zero constraints. Use \code{"all"} to estimate
#'   all methods and compare.
#' @param forced_randomize A logical value that indicates, in the forced-choice
#'   setting, whether the "left" and "right" profiles should be randomized for
#'   each task. The default is \code{FALSE}.
#' @param tau_method A character value indicating the method for dealing with
#'   binding restrictions, i.e. numerically infinite \eqn{E[1/\tau^2]}. The two
#'   options are \code{"nullspace"} (i.e. perform inference assuming this
#'   restriction binds) or \code{"clip"} (set to a large value
#'   \code{tau_truncate}). The default is \code{"nullspace"}.
#' @param tau_truncate A numerical value to either truncate \eqn{E[1/\tau^2]}
#'   (i.e. set maximum \eqn{E[1/\tau^2]} in the E-Step for updating \eqn{\beta})
#'   if \code{tau_method = "clip"} or a threshold by which to declare that two
#'   levels are fused if 
#'   \code{tau_method="nullspace"}. The default is 1e6.
#' @param tau_stabilization An integer value of the number of steps to perform
#'   with \code{tau_method="clip"} before using the provided setting. The
#'   default is 5.
#' @param force_reset A logical argument about how the nullspace is computed. If
#'   \code{tau_method="nullspace"}, it forces nullspace to be estimated directly
#'   from all binding restrictions at each iteration versus the default method
#'   that updates the existing basis when possible. Default is \code{FALSE}.
#' @param calc_df A logical value for whether to calculate degrees of freedom of
#'   final model. The default is \code{TRUE}.
#' @param calc_se A logical value for whether standard errors of final model.
#'   The default is \code{TRUE}.
#' @param quiet_tictoc A logical value for whether to \emph{not} print
#'   information about the timing of the model. The default is \code{TRUE}.
#' @param override_BR A logical value for whether to ignore Bondell and Reich
#'   style-weights. The default is \code{FALSE}. If \code{TRUE} is provided,
#'   \code{sqrt(L) * (L + 1)} is used, where \code{L} is the number of factor
#'   levels.
#' @param debug A logical value for whether the algorithm should be debugged.
#'   The default is \code{FALSE}. In particular, it will verify that the
#'   log-posterior increases at each (intermediate) step and throw an exception
#'   otherwise.
#' @examples 
#' str(FactorHet_control())
#' 
#' @encoding UTF-8
#' @references 
#' 
#' Bondell, Howard D., and Brian J. Reich. 2009. "Simultaneous Factor Selection
#' and Collapsing Levels in ANOVA." Biometrics 65(1): 169-177.
#' 
#' Goplerud, Max, Kosuke Imai, and Nicole E. Pashley. 2025. "Estimating
#' Heterogeneous Causal Effects of High-Dimensional Treatments: Application to
#' Conjoint Analysis." arxiv preprint: \url{https://arxiv.org/abs/2201.01357}
#' 
#' Post, Justin B., and Howard D. Bondell. 2013. "Factor Selection and
#' Structural Identification in the Interaction ANOVA Model." \emph{Biometrics}
#' 69(1):70-79.
#' 
#' Lim, Michael, and Trevor Hastie. 2015. "Learning Interactions via Hierarchical
#' Group-Lasso Regularization." \emph{Journal of Computational and Graphical
#' Statistics} 24(3):627-654.
#' 
#' Städler, Nicolas, Peter Bühlmann, and Sara Van De Geer. 2010.
#' "l1-penalization for Mixture Regression Models." \emph{Test} 19(2):209-256.
#' 
#' Yan, Xiaohan and Jacob Bien. 2017. "Hierarchical Sparse Modeling: A Choice of
#' Two Group Lasso Formulations." \emph{Statistical Science} 32(4):531–560.
#' 
#' @return \code{FactorHet_control} returns a named list containing the elements
#'   listed in "Arguments".
#' @export
FactorHet_control <- function(
  iterations = 1000, 
  maxit_pi = NULL, 
  optim_phi_controls = list(method = 'lib_lbfgs'),
  prior_var_phi = 4, 
  prior_var_beta = Inf, gamma = 1, repeat_beta = 1, 
  adaptive_weight = 'B&R', init_method = 'short_EM',   
  return_data = FALSE, log_method = 'log_ginv',
  tolerance.parameters = 1e-5, tolerance.logposterior = 1e-5,
  rare_threshold = 5, rare_verbose = 1, 
  beta_method = 'cpp',
  beta_cg_it = 25, lambda_scale = 'N', 
  weight_dlist = FALSE,
  do_SQUAREM = TRUE, step_SQUAREM = NULL,
  backtrack_SQUAREM = 10, df_method = 'EM',
  forced_randomize = FALSE, single_intercept = NULL,
  tau_method = 'nullspace', tau_stabilization = 5, tau_truncate = 1e6, 
  debug = FALSE, force_reset = FALSE,
  calc_df = TRUE, 
  calc_se = TRUE,
  quiet_tictoc = TRUE, override_BR = FALSE){
  
  stopifnot(lambda_scale %in% c('N', 'unity', 'root_N'))
  if (lambda_scale == 'N'){
    lambda_scale <- function(lambda,N){lambda * N}
  }else if (lambda_scale == 'root_N'){
    lambda_scale <- function(lambda,N){lambda * sqrt(N)}
  }else if (lambda_scale == 'unity'){
    lambda_scale <- function(lambda, N){lambda}
  }
  stopifnot(df_method %in% c('all', 'EM', 'free_param', 'IRLS'))
  stopifnot(log_method %in% c('standard', 'log_0', 'log_ginv', 'log_random'))
  stopifnot(tau_method %in% c('clip', 'nullspace'))
  
  if (is.null(optim_phi_controls$method)){
    stop('optim_phi_controls must, at a mimum, contain list(method = "XX")')
  }else if (optim_phi_controls$method == 'lib_lbfgs'){
    check <- requireNamespace('lbfgs', quietly = TRUE)
    if (!check){
      stop('lbfgs must be installed for optim_phi_controls=list(method="lib_lbfgs")')
    }
    rm(check)
  }
  if (!is.null(step_SQUAREM)){
    if (step_SQUAREM > -1){
      stop('step_SQUAREM must be NULL or smaller than -1')
    }
  }
  
  output <- mget(ls())
  class(output) <- 'FactorHet_control'
  return(output)
}

#' Arguments for initializing FactorHet
#' 
#' A set of arguments that govern the initialization of \code{\link{FactorHet}}.
#' Use \code{\link{FactorHet_control}} to set arguments around estimation.
#' \code{\link{FactorHet_mbo}} ignores many of these arguments as it uses a
#' single fixed initialization set by \code{\link{FactorHet_mbo_control}}. All
#' arguments have default values.
#' 
#' @param nrep An integer value of the number of random iterations or runs of
#'   "short EM" should be used. The default value is 5.
#' @param short_EM A logical value indicating whether "short EM" should be used.
#'   The default value is \code{FALSE}. \code{TRUE} indicates a "short EM"
#'   should be followed. That is, run multiple short runs of EM with random
#'   initializations and then proceed with the best for full initialization.
#'   Biernacki et al. (2003) provides more discussion. If
#'   \code{\link{FactorHet_control}} has \code{init_method = "short_EM"}, this
#'   will override this setting.
#' @param short_EM_init An argument that sets the initialization procedure for
#'   "short EM". It must be some non-deterministic procedure that is valid in
#'   \code{\link{FactorHet_control}}. The default is \code{"random_member".}
#' @param short_EM_it A numerical value of the number of iterations to use for
#'   each "short" run of the EM algorithm. The default is 40.
#' @param short_EM_pi An argument for the maximum number of iterations for the
#'   moderator updates to use for each "short" run of the EM algorithm. The
#'   default is \code{NULL}.
#' @param short_EM_beta_method An argument for the update method for \eqn{\beta}
#'   to use for each "short" run of the EM algorithm. The default is
#'   \code{"cpp"}.
#' @param short_EM_cg_it An argument for the number of conjugate gradient
#'   iterations to use if \code{short_EM_beta_method = "cg"}.
#' @param force_rep A logical value for whether to repeat the algorithm if
#'   \code{K=1}. The default is \code{FALSE} and it should be used only for
#'   debugging.
#' @param verbose A logical value to print more information about the progress
#'   of each iteration. The default is \code{FALSE}.
#' @param debug_repeat A logical value for whether to debug the repeated runs.
#'   The default is \code{FALSE}.
#' @param plot_repeat A logical value for whether to plot the trajectory of the
#'   log-posterior for each run. The default is \code{FALSE}.
#' @param return_all A logical value for whether to return all repetitions of
#'   the model versus the one with the highest log-posterior. The default is
#'   \code{FALSE}.
#' 
#' @examples 
#' str(FactorHet_init())
#' 
#' @return \code{FactorHet_init} returns a named list containing the elements
#'   listed in "Arguments".
#' @encoding UTF-8
#' @references 
#' Biernacki, Christophe, Gilles Celeux, and Gérard Govaert. "Choosing Starting
#' Values for the EM algorithm for Getting the Highest Likelihood in
#' Multivariate Gaussian Mixture Models." 2003. \emph{Computational Statistics &
#' Data Analysis}. 41(3-4):561-575.
#' @export
FactorHet_init <- function(short_EM = FALSE, short_EM_it = 40, 
                           short_EM_init = 'random_member', short_EM_pi = NULL,
                           force_rep = FALSE, verbose = FALSE, 
                           short_EM_beta_method = 'cpp', 
                           short_EM_cg_it = 10, nrep = 5, 
                           debug_repeat = FALSE, 
                           plot_repeat = FALSE, return_all = FALSE){
  output <- mget(ls())
  class(output) <- 'FactorHet_init'
  return(output)
}
