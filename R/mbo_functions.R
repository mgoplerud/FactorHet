#' @rdname FactorHet
#' @param mbo_control A list of control parameters for MBO; see
#'   \link{FactorHet_mbo_control} for more information.
#' @export
FactorHet_mbo <- function(formula, design, K, moderator = NULL, weights = NULL, 
    group = NULL, task = NULL, choice_order = NULL, control = FactorHet_control(), 
    initialize = FactorHet_init(), mbo_control = FactorHet_mbo_control()){
  
  if (mbo_control$ic_method != control$df_method){
    message('mbo_control$ic_method does not equal control$df_method; over-riding to use mbo_control$ic_method.')
    control$df_method <- mbo_control$ic_method
  }
  
  # Prepare the regression data so it does not need to be
  # created for each model evaluation 
  message('Preparing Data')
  design <- prepare_regression_data(formula = formula, 
     design = design, weights = weights,
     forced_randomize = control$forced_randomize,
     moderator = moderator,
     group = group, task = task, choice_order = choice_order,
     single_intercept = control$single_intercept, 
     rare_threshold = control$rare_threshold,
     rare_verbose = control$rare_verbose
  )
  if (length(setdiff(names(mbo_control), names(FactorHet_mbo_control())))){
    stop('Some arguments from "FactorHet_mbo_control()" missing from "mbo_control".')
  }
  if (!is.null(mbo_control$mbo_initialize)){
    message('Computing Initalization')
    init <- prepare_deterministic_init(data = design, K = K, 
        method = mbo_control$mbo_initialize, 
        iterations = mbo_control$mm_init_iterations)
    init <- list('group_E.prob' = init$group_E.prob)
    control$init_method <- init
    rm(init)
  }else{
    message('Running MBO without deterministic initalization is not advised')
  }

  # Package all of the arguments needed to fit the model  
  mbo_control <- mbo_control[names(mbo_control) != 'mbo_initialize']
  mbo_control <- mbo_control[names(mbo_control) != 'mm_init_iterations']
  mbo_args <- c(
    list('FH_args' = mget(setdiff(ls(), c('mbo_control')))),
    mbo_control
  )
  
  
  message('Fitting MBO')
  fit_mbo <- do.call('internal_FH_mbo', mbo_args)
  # Repackage the data to be a FactorHet object 
  # plus some additional elements from the MBO
  fit_mbo$model$MBO_output <- fit_mbo[c('path', 'obj')]
  fit_mbo <- fit_mbo$model
  fit_mbo$MBO_output$control <- mbo_control
  class(fit_mbo) <- c('FactorHet_MBO', 'FactorHet')
  return(fit_mbo)
}

#' Control for Model-Based Optimization
#' 
#' \code{FactorHet_mbo_control} is used to adjust the settings for the MBO
#' (Model-Based Optimization). All arguments have default values. This relies
#' heavily on options from the \code{mlrMBO} package so please see this for more
#' detailed discussion.
#' 
#' @param mbo_type Options \code{"sparse"} (default) or \code{"ridge"}. The
#'   former (\code{"sparse"}) uses the structured sparse penalty discussed in Goplerud et al.
#'   (2022). \code{"ridge"} does a ridge regression.
#' @param mbo_initialize Initialization method for each MBO proposal. See
#'   "Details" for more discussion.
#' @param mm_init_iterations How many options to use if Murphy/Murphy
#'   initialization is required? See "Details" for discussion. Default of
#'   \code{NULL} uses default values of 100 if probabilistic and 50 if
#'   deterministic.
#' @param mbo_method Function used to propose new values of the regularization
#'   parameters. See information from \code{mlr3} for more details.
#' @param final_method How should the final regularization parameter be
#'   selected? Default is \code{"best_predicted"}. As documented in
#'   \code{makeMBOcontrol} for \code{final.method}, this uses the final model to
#'   predict all points visited and use the best one. Alternative options are
#'   those available in \code{makeMBOcontrol} and include \code{"last.proposed"}
#'   and \code{"best.true.y"}.
#' @param iters How many proposals to do after initialization?
#' @param criterion Which criterion to minimize? Options are \code{"BIC"}
#'   (default), \code{"AIC"}, \code{"GCV"}, or \code{"BIC_group"}.
#'   \code{"BIC_group"} counts the number of observations as the number of
#'   individuals (e.g. in the case of repeated observations per person).
#' @param ic_method Method for calculating degrees of freedom: \code{"EM"}
#'   (default), \code{"IRLS"}, \code{"free_param"}. See \link{FactorHet_control}
#'   for more information.
#' @param se_final Should standard errors be calculated for the final model?
#'   Default of \code{TRUE}.
#' @param mbo_design How to design the initial proposals for the MBO? See
#'   "Details".
#' @param mbo_range The range of values to consider on log10(lambda), before
#'   standardization. See "Details". Default is \code{c(-5,0)}.
#' @param fast_estimation Default \code{NULL} uses same convergence criterion
#'   for all models. See "Details" for more information.
#' @param verbose Show the progress of MBO at each iteration; default
#'   \code{FALSE}.
#' 
#' @details 
#' 
#' \bold{Initialization}: \code{FactorHet_mbo} relies on the same initialization for each
#' attempt. The default procedure (\code{"mm_mclust_prob"}) is discussed in detail in
#' the appendix of Goplerud et al. (2022) and builds on Murphy and Murphy (2020).
#' In brief, it deterministically initializes cluster memberships using only the
#' moderators (e.g. using \code{"mclust"}). Using those memberships, it uses an
#' EM algorithm (with probabilistic assignment, if \code{"prob"} is specified,
#' or hard assignment otherwise) for a few steps with only the main effects to
#' update the proposed cluster memberships. If the warning appears that
#' "Murphy/Murphy initialization did not fully converge" , this mean that this
#' initial step did not fully converge. The number of iterations could be
#' increased using \code{mm_init_iterations} if desired, although benefits are
#' usually modest beyond the default settings. These memberships are then used
#' to initialize the model at each proposed regularization value.
#' 
#' The options available are \code{"spectral"} and \code{"mclust"} that use
#' \code{"spectral"} or \code{"mclust"} on the moderators with no Murphy/Murphy style
#' tuning. Alternatively, \code{"mm_mclust"} and \code{"mm_spectral"} combine
#' the Murphy/Murphy tuning upon the corresponding initial deterministic
#' intialization (e.g. spectral or \code{"mclust"}). These use hard assignment
#' at each step and likely will converge more quickly although a hard initial
#' assignment may not be desirable. Adding the suffix \code{"_prob"} to the
#' \code{"mm_*"} options uses a standard (soft-assignment) EM algorithm during
#' the Murphy/Murphy tuning.
#' 
#' If one wishes to use a custom initialization for MBO, then set
#' \code{mbo_initialize=NULL} and provide an initialization via
#' \code{FactorHet_control}. It is strongly advised to use a deterministic
#' initialization if done manually, e.g. by providing a list of initial
#' assignment probabilities for each group.
#' 
#' \bold{Design of MBO Proposals}: The MBO procedure works as follows; there are
#' some initial proposals that evaluated in terms of the criterion. Given those
#' initial proposals, there are \code{iters} attempts to improve the criterion
#' through methods described in detail in \code{mlrMBO} (Bischl et al. 2018). A default
#' of 11 seems to work well, though one can examine \code{visualize_MBO} after
#' estimation to see how the criterion varied across the proposals.
#' 
#' By default, the regularization parameter is assumed to run from -5 to 0 on
#' the log10 scale, before standardizing by the size of the dataset. We found
#' this to be reasonable, but it can be adjusted using \code{mbo_range}.
#' 
#' It is possible to calibrate the initial proposals to help the algorithm find
#' a minimum of the criterion more quickly. This is controlled by
#' \code{mbo_design} which accepts the following options. Note that a manual
#' grid search can be provided using the \code{data.frame} option below.
#' 
#' \itemize{ 
#' \item{Scalar: }{By default, this is initialized with a scalar
#' (-1.5) that is the log10 of lambda, before standardization as discussed in
#' \code{FactorHet_control}. For a scalar value, four proposals are generated
#' that start with the scalar value and adjust it based on the level of sparsity
#' of the initial estimated model. This attempts to avoid initializations that
#' are too dense and thus are very slow to estimate, as well as ones that are
#' too sparse.} 
#' \item{"random": }{If the string "random" is provided, this
#' follows the default settings in \code{mlrMBO} and generates random
#' proposals.} 
#' \item{data.frame: }{A custom grid can be provided using a
#' data.frame that has two columns (\code{"l"} and \code{"y"}). \code{"l"}
#' provides the proposed values on the log10 lambda scale (before
#' standardization). If the corresponding BIC value is known, e.g. from a prior
#' run of the algorithm, the column \code{"y"} should contain this value. If it
#' is unknown, leave the value as \code{NA} and the value will be estimated.
#' Thus, if a manual grid search is desired, this can be done as follows. Create
#' a data.frame with the grid values \code{"l"} and all \code{"y"} as NA. Then,
#' set \code{iters = 0} to do no estimation \emph{after} the grid search. } 
#' }
#' 
#' \bold{Estimation}: Typically, estimation proceeds using the same settings for
#' each MBO proposal and the final model estimated given the best regularization
#' value (see option \code{final_method} for details). However, if one wishes to
#' use a lower convergence criterion for the MBO proposals to speed estimation,
#' this can be done using the \code{fast_estimation} option. This proceeds by
#' giving a named list with two members \code{"final"} and \code{"fast"}. Each
#' of these should be a list with two elements \code{"tolerance.logposterior"}
#' and \code{"tolerance.parameters"} with the corresponding convergence
#' thresholds. \code{"final"} is used for the final model and \code{"fast"} is
#' used for evaluating all of the MBO proposals.
#' 
#' @references 
#' 
#' Bischl, Bernd, Jakob Richter, Jakob Bossek, Daniel Horn, Janek Thomas and
#' Michel Lang. 2018. "mlrMBO: A Modular Framework for Model-Based Optimization
#' of Expensive Black-Box Functions." arxiv preprint:
#' \url{https://arxiv.org/abs/1703.03373}
#' 
#' Goplerud, Max, Kosuke Imai, and Nicole E. Pashley. 2022. "Estimating
#' Heterogeneous Causal Effects of High-Dimensional Treatments: Application to
#' Conjoint Analysis." arxiv preprint: \url{https://arxiv.org/abs/2201.01357}
#' 
#' Murphy, Keefe and Thomas Brendan Murphy. 2020. "Gaussian Parsimonious
#' Clustering Models with Covariates and a Noise Component." \emph{Advances in
#' Data Analysis and Classification} 14:293??? 325.
#' 
#' 
#' @examples 
#' str(FactorHet_mbo_control())
#' 
#' @export
FactorHet_mbo_control <- function(
  mbo_type = c('sparse', 'ridge'),
  mbo_initialize = 'mm_mclust_prob',
  mm_init_iterations = NULL, mbo_range = c(-5, 0),
  mbo_method = 'regr.bgp', final_method = 'best.predicted', iters = 11, 
  criterion = c('BIC', 'AIC', 'GCV', 'BIC_group'),
  ic_method = c('EM', 'IRLS', 'free_param'), 
  se_final = TRUE, mbo_design = -1.5,
  fast_estimation = NULL,
  verbose = FALSE){
  
  mbo_type <- match.arg(mbo_type)
  criterion <- match.arg(criterion)
  ic_method <- match.arg(ic_method)
  # Verify this option is specified correctly
  if (!is.null(fast_estimation)){
    check_valid <- sapply(c('final', 'fast'), FUN=function(s){
      (inherits(fast_estimation[[s]], 'list')) &
      (length(
        setdiff(names(fast_estimation[[s]]), 
                c('tolerance.logposterior', 'tolerance.parameters'))
      ) == 0) &
      (length(fast_estimation[[s]]) == 2) &
      (all(sapply(fast_estimation[[s]], is.numeric))) &
      (all(lengths(fast_estimation[[s]]) == 1)) 
    })
    if (length(check_valid) != 2 | !all(check_valid)){
      stop('fast_estimation must be NULL or a 2-length list described in the documentation.')
    }
    rm(check_valid)
  }
  
  if (mbo_method == 'regr.bgp'){
    if (!requireNamespace('tgp', quietly = FALSE)){
      warning('mbo_method="regr.bgp" needs "tgp" installed.')
    }
  }
  output <- mget(ls())
  return(output)
}
  
#' @import mlrMBO
#' @importFrom mlr makeLearner configureMlr
#' @importFrom smoof makeSingleObjectiveFunction
#' @importFrom ParamHelpers makeParamSet makeNumericParam generateDesign getParamSet
internal_FH_mbo <- function(FH_args,
                          mbo_method = 'regr.bgp', mbo_type,
                          fast_estimation = NULL, mbo_range,
                          final_method = 'best.predicted', iters = 11, 
                          criterion = c('BIC', 'AIC', 'GCV', 'BIC_group'),
                          ic_method = c('EM', 'IRLS', 'free_param'), 
                          se_final = TRUE, mbo_design = -1.5,
                          verbose = FALSE){

  ic_method <- match.arg(ic_method)
  criterion <- match.arg(criterion)
  if (is.null(FH_args$control)){
    FH_args$control <- FactorHet_control()
  }
  if (!is.null(fast_estimation)){
    FH_args$control[
      c('tolerance.logposterior', 'tolerance.parameters')
    ] <- fast_estimation$fast
    
  }
  FH_args$control$df_method <- ic_method
  
  MBO_func <- function(x, arg, full_output = FALSE){
    
    if (mbo_type == 'sparse'){
      arg$lambda <- 10^(x$l)
    }else if (mbo_type == 'ridge'){
      arg$lambda <- 0
      arg$control$prior_var_beta <- 10^(x$l)
    }else{stop('mbo_type must be sparse or ridge')}
    
    arg$control$calc_se <- FALSE
    
    output <- do.call('FactorHet', arg)
    IC <- output$information_criterion
    if (full_output == TRUE){
      return(output)
    }else{
      return(IC[[criterion]][IC$method == ic_method])
    }
  }
  
  copy_args <- FH_args

  if (inherits(mbo_design, 'data.frame')){
    design_type <- 'provided'
  }else if (inherits(mbo_design, 'numeric')){
    if (is.numeric(mbo_design) & length(mbo_design) == 1){
      design_type <- 'smart'
    }else{
      stop('mbo_design must be (i) data.frame with l and y columns, (ii) a single number log_10(lambda), or (iii) "random"')
    }
  }else if (identical(mbo_design, 'random')){
    design_type <- 'random'
  }else{
    stop('mbo_design must be (i) data.frame with l and y columns, (ii) a single number log_10(lambda), or (iii) "random"')
  }

  if (design_type == 'smart'){
    if (verbose){message('Generating Four Proposals for MBO')}
    initial_propose <- mbo_design
    #Propose a model that is likely to be fairly sparse
    time_init <- proc.time()
    initial_model <- MBO_func(list(l = initial_propose), 
      copy_args, full_output = TRUE)  
    time_init <- time_init - proc.time()
    if (verbose){message(paste0('Time for Propose 1: ', time_init[3]/60))}
    initial_sparse <- mean(abs(coef(initial_model)) < 1e-4)
    initial_IC <- initial_model$information_criterion
    if (verbose){print(initial_IC)}
    rm(initial_model); gc()
    initial_IC <- initial_IC$BIC[initial_IC$method == ic_method]
    if (initial_sparse > 0.90){
      second_propose <- initial_propose - 1
    }else{
      second_propose <- initial_propose - 1/2
    }
  
    time_second <- proc.time()
    second_model <- MBO_func(list(l = second_propose), copy_args, full_output = TRUE)  
    time_second <- time_second - proc.time()
    if (verbose){message(paste0('Time for Propose 2: ', time_second[3]/60))}
    
    second_sparse <- mean(abs(coef(second_model)) < 1e-4)
    second_IC <- second_model$information_criterion
    rm(second_model); gc()
    second_IC <- second_IC$BIC[second_IC$method == ic_method]
    percent_IC <- (initial_IC - second_IC)/initial_IC * 100
  
    #If still sparse, jump FARTHER
    if (second_sparse > 0.90){
      third_propose <- second_propose - 1
    }else{
      #If at least 1% improvement in the IC, 
      #decrease lambda by 50%
      if (percent_IC > 1){
        third_propose <- second_propose + log(0.5)/log(10)
      }else{
        #If going more denser is likely WORSE AND
        #initial model is not very sparse go SPARSER
        if (initial_sparse < 0.8){
          third_propose <- initial_propose + 1/2
        }else{
          #If going dense is likely worse and initial was sparse
          #try half way between
          third_propose <- (initial_propose + second_propose)/2
        }
      }
    }
    time_third <- proc.time()
    third_model <- MBO_func(list(l = third_propose), copy_args, full_output = TRUE)  
    time_third <- time_third - proc.time()
    if (verbose){message(paste0('Time for Propose 3: ', time_third[3]/60))}
    
    third_sparse <- mean(abs(coef(third_model)) < 1e-4)
    third_IC <- third_model$information_criterion
    rm(third_model); gc()
    third_IC <- third_IC$BIC[third_IC$method == ic_method]
    
    l_init <- c(initial_propose, second_propose, third_propose, NA)
    y_init <- c(initial_IC, second_IC, third_IC, NA)
    
    if (verbose){print(cbind(y_init, l_init))}
    
    if (third_sparse > 0.90){#if still sparse try YET AGAIN
      fourth_propose <- third_propose - 1
    }else{
      if (third_IC > second_IC){#If third is *worse* then do an average
        fourth_propose <- (third_propose + second_propose)/2
      }else{# If third is again better than second, make again 25% worse
        fourth_propose <- third_propose + log(0.25)/log(10)
      }
    }
    time_fourth <- proc.time()
    fourth_model <- MBO_func(list(l = fourth_propose), copy_args, full_output = TRUE)  
    time_fourth <- time_fourth - proc.time()
    if (verbose){message(paste0('Time for Propose 4: ', time_fourth[3]/60))}
    
    fourth_sparse <- mean(abs(coef(fourth_model)) < 1e-4)
    fourth_IC <- fourth_model$information_criterion
    rm(fourth_model); gc()
    fourth_IC <- fourth_IC$BIC[fourth_IC$method == ic_method]
    l_init[4] <- fourth_propose
    y_init[4] <- fourth_IC
    init_sparse <- c(initial_sparse, second_sparse, third_sparse, fourth_sparse)
    init_time <- c(time_init[3]/60, time_second[3]/60, time_third[3]/60, time_fourth[3]/60)
  }
  
  
  obj_FH <- makeSingleObjectiveFunction(
    name = "FactorHet (MBO)",
    fn = MBO_func,
    par.set = makeParamSet(
      makeNumericParam("l", lower = mbo_range[1], upper = mbo_range[2])
    ),
    has.simple.signature = FALSE,
    minimize = TRUE
  )
  
  surr.learner = makeLearner(mbo_method, predict.type = 'se')

  if (iters != 0){
    ctrl = makeMBOControl(final.method = final_method)
    ctrl = setMBOControlTermination(ctrl, iters = iters)
    ctrl = setMBOControlInfill(ctrl, crit = makeMBOInfillCritEI())
    configureMlr(show.learner.output = FALSE)
  }  
  
  if (design_type == 'random'){
    if (verbose){message('Generating Four (Random) Proposals')}
    design2 = generateDesign(n = 4, par.set = getParamSet(obj_FH))
    
    y_init <- sapply(design2$l, FUN=function(i){MBO_func(x = list(l = i), arg = FH_args)})
    design2$y <- y_init
    if (verbose){
      message('Summary of Initial Design')
      print(design2)
    }
  }else if (design_type == 'smart'){
    design2 <- data.frame(l = l_init, y = y_init, 
        sparse = init_sparse, time = init_time)
    if (verbose){
      message('Summary of Initial Design')
      print(design2)
    }
  }else if (design_type == 'provided'){
    design2 <- mbo_design
    # Fill in the missing "y"
    if (!('y' %in% names(design2))){
      design2$y <- NA
    }
    fill_y <- which(is.na(design2$y))
    if (length(fill_y) != 0){
      design2$y[fill_y] <- sapply(design2$l[fill_y], FUN=function(i){MBO_func(x = list(l = i), arg = FH_args)})
    }
  }else{stop('invalid type for mbo_design')}
  
  if (iters > 0){
    run2 <- mbo(obj_FH, design = design2[,c('l', 'y')], more.args = list(arg = copy_args),
                learner = surr.learner, control = ctrl, show.info = TRUE)
    path_2 <- as.data.frame(run2$opt.path)
  }else{
    run2 <- list(x = data.frame(l = design2$l[which.min(design2$y)]))
    path_2 <- design2
  }

  final_args <- FH_args
  if (!is.null(fast_estimation)){
    
    final_args$control[
      c('tolerance.logposterior', 'tolerance.parameters')
    ] <- fast_estimation$final
    
    FH_args$control[
      c('tolerance.logposterior', 'tolerance.parameters')
    ] <- fast_estimation$fast
    
  }
  final_args$control$calc_se <- se_final
  
  if (mbo_type == 'sparse'){
    if (final_method == 'best.predicted'){
      final_args$lambda <- 10^run2$x$l  
    }else{
      final_args$lambda <- 10^path_2$l[which.min(path_2$y)]
    }
  }else{
    final_args$lambda <- 0
    if (final_method == 'best.predicted'){
      final_args$control$prior_var_beta <- 10^(run2$x$l)
    }else{
      final_args$control$prior_var_beta <- 10^(path_2$l[which.min(path_2$y)])
    }
  }
  final_model <- do.call('FactorHet', final_args)
  
  if (iters == 0){
    run2 <- NULL
  }
  
  return(list(model = final_model, path = path_2, obj = run2))
}

