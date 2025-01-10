#' Predict after using FactorHet
#' 
#' Predicted values for choosing particular profiles based on
#' \code{\link{FactorHet}}.
#' 
#' @param object A model estimated using \code{\link{FactorHet}} or
#'   \code{\link{FactorHet_mbo}}.
#' @param newdata A dataset on which to generate predictions. The default,
#'   \code{NULL}, uses the estimation data.
#' @param type An argument that specifies how to deal with group-membership
#'   probabilities when making predictions. The default is \code{"posterior"}
#'   which use the posterior probabilities for each observation in the training
#'   data for weighting the groups. If \code{"posterior_predictive"} is used,
#'   the group membership probabilities implied by the moderators, i.e.
#'   \eqn{\hat{\pi}_k(X_i)}, will be used. If an observation in \code{newdata}
#'   is not in the estimation data (i.e., its value of \code{group}) is not
#'   found, then \code{"posterior_predictive"}, i.e. \eqn{\hat{\pi}_k(X_i)}, is
#'   used.
#' @param by_group A logical value as to whether the predictions should be
#'   returned for each group or whether a weighted averaged based on the group
#'   membership probabilities (specified by \code{type}) should be reported. The
#'   default is \code{FALSE}.
#' @param return A character value that determines the type of prediction
#'   return. The default is \code{"prediction"} that returns the predicted
#'   probability. The option \code{"detailed"} returns a variety of additional
#'   information. This is mostly called internally for other functions such as
#'   \code{\link{AME}} or \code{\link{margeff_moderators}}.
#' @param ... Miscellaneous options used internally and not documented.
#' @examples 
#' data(immigration)
#' set.seed(1)
#' # Fit a model once for simplicity
#' fit <- FactorHet(Chosen_Immigrant ~ Plans + Ed + Country,
#'  design = immigration, lambda = 1e-4,
#'  # Randomly initialize, do only one iteration for speed
#'  init = FactorHet_init(nrep = 1),
#'  control = FactorHet_control(init = 'random_member'),
#'  K = 2, group = ~ CaseID, task = ~ contest_no, 
#'  choice_order = ~ choice_id)
#' immigration$pred_FH <- predict(fit)  
#' @importFrom stats predict coef formula plogis delete.response update.formula
#'   na.omit qlogis
#' @return Returns an estimate of the predicted probability of choosing a
#'   profile for each observation. "Arguments" outlines different behavior if
#'   certain options are chosen.
#' @export 
predict.FactorHet <- function(object, newdata = NULL, type = 'posterior', 
                              by_group = FALSE, 
                              return = 'prediction', ...){
  options <- list(...)
  
  type <- match.arg(type, c('posterior', 'posterior_predictive', 'both'))
  #Only relevant ... is "require_response" but not to be ever adjusted by user.
  if (length(setdiff(names(options), c('return_task', 'quick_exit',
                                       'override_weights', 
                                       'ignore_moderator', 'ignore_treatment',
                                       'calc_gradient', 'require_response'))) > 0){
    stop('... is ignored for predicting with FactorHet')
  }
  
  if (is.null(newdata)){newdata <- object$internal_parameters$data$design}
  if (is.null(newdata)){stop('newdata must be provided if return_data is NA')}
  
  if (!is.null(options$calc_gradient)){
    calc_gradient <- options$calc_gradient
  }else{
    calc_gradient <- FALSE
  }

  if (!is.null(options$override_weights)){
    override_weights <- options$override_weights
  }else{
    override_weights <- TRUE
  }
  if (!is.null(options$quick_exit)){
    quick_exit <- options$quick_exit
  }else{
    quick_exit <- FALSE
  }
  if (!is.null(options$return_task)){
    return_task <- options$return_task
  }else{
    return_task <- FALSE
  }
  if (!is.null(options$require_response)){
    delete_response <- !options$require_response
    type <- 'both'
  }else{
    delete_response <- TRUE
  }

  formula_object <- formula(object)
  
  if (override_weights){
    formula_object$weights <- NULL
  }

  if (isTRUE(options$ignore_moderator)){
    formula_object$mod <- ~ 1
  }
  if (isTRUE(options$ignore_treatment)){
    formula_object$het <- ~ 1
    formula_object$other_parameters <- NULL
  }
  
  prep_newdata <- prepare_formula(fmla_main = formula_object$het, 
    fmla_moderator = formula_object$mod, design = newdata, 
    group = formula_object$other_parameters$name_group, 
    task = formula_object$other_parameters$name_task, 
    weights = formula_object$weights,
    choice_order = formula_object$other_parameters$name_choice_order ,
    delete_response = delete_response)
  
  if (!isTRUE(options$ignore_treatment)){
    new_X <- create_data(design = prep_newdata$design, 
     penalty_for_regression = object$internal_parameters$penalty,
     verif_row = FALSE,
     remove_cols = names(object$internal_parameters$rare$rare_col))
    new_group <- data.frame(prep_newdata$design[prep_newdata$name_group], stringsAsFactors = F)
    new_task <- data.frame(prep_newdata$design[prep_newdata$name_task], stringsAsFactors = F)
    new_choice_order <- data.frame(prep_newdata$design[prep_newdata$name_choice_order], stringsAsFactors = F)
    
    if (ncol(new_group) == 0){
      new_group <- NULL
    }else{
      new_group <- as.character(as.vector(new_group[,1]))
    }
    if (ncol(new_task) == 0){
      null_task <- TRUE
    }else{
      null_task <- FALSE
      new_task <- as.character(as.vector(new_task[,1]))
    }
    if (ncol(new_choice_order) == 0){
      new_choice_order <- NULL
    }else{
      new_choice_order <- as.character(as.vector(new_choice_order[,1]))
    }
    
    conjoint_names <- prep_newdata[c('name_group', 'name_task', 'name_choice_order')]
    factor_names <- prep_newdata$factor_names
    #Verify same factors
    stopifnot(identical(factor_names, names(object$internal_parameters$factor_levels)))
    
  }else{
    new_X <- NULL
    new_group <- new_task <- new_choice_order <- NULL
  }
  
  if (!delete_response){
    new_y <- prep_newdata$design[[prep_newdata$outcome]]
  }else{
    new_y <- NULL
  }
  
  
  if (is.null(new_group)){
    new_group <- 1:nrow(prep_newdata$design)
    null_group <- TRUE
  }else{
    null_group <- FALSE
  }
  unique_new_group <- unique(new_group)
  
  
  if (any(is.na(new_group))){stop('No Missingness Allowed for Group')}
  
  if (!by_group){#Get predictions averaged across, need info on moderators
    
    #Build the moderator variables
    if (object$K == 1){
      formula_object$mod <- ~ 1
    }
    if (isTRUE(options$ignore_moderator)){
      stop('Set up option')
    }
    new_concom <- suppressMessages(create_moderator(design = prep_newdata$design, 
        args = object$internal_parameters$W$args, check_rank = FALSE,
        moderator = formula_object$mod, group = new_group, 
        unique_group = unique_new_group))
    new_W <- new_concom$W
    new_colMeans_W <- new_concom$colMeans
    ncol_new_W <- new_concom$ncol_W
    
    
    checksum_names <- object$internal_parameters$W$args$names
    if (length(checksum_names) == 0){checksum_names <- NULL}
    if (!identical(colnames(new_W), checksum_names)){
      easy_message('Non Identical Names in OOS W and W')
      easy_message('OOS NAMES')
      easy_message(colnames(new_W))
      easy_message('MODEL NAMES')
      easy_message(object$internal_parameters$W$args$names[-1])
      stop('error: inspect prior messages')
    }
    
  }else{
    new_W <- matrix(1, nrow =length(unique_new_group))
  }
  
  new_weights <- prep_newdata$weights
  
  if (quick_exit){return(new_W)}
  
  if (!is.null(new_choice_order)){
    if (null_group){stop('Group must be provided for forced choice')}
    
    forced_output <- make_forced_choice(y = new_y, X = new_X, group = new_group,
        task = new_task, choice_order = new_choice_order, 
        unique_choice = object$internal_parameters$unique_choice,
        estimate_intercept = TRUE, weights = new_weights,
        randomize = object$internal_parameters$control$forced_randomize)
    new_X <- forced_output$X
    new_y <- forced_output$y
    new_group <- forced_output$group
    new_task <- forced_output$task
    new_weights <- forced_output$weights
    
    rm(forced_output)
  }
  
  if (isTRUE(options$ignore_treatment)){
    new_group_mapping <- new_weights_W <- new_weights_sq_W <- NULL
  }else{
    new_group_mapping <- sparseMatrix(i = 1:nrow(new_X), j = match(new_group, unique_new_group), x = 1)
    new_weights_W <- Diagonal(x = 1/colSums(new_group_mapping)) %*% t(new_group_mapping) %*% new_weights
    new_weights_sq_W <- Diagonal(x = 1/colSums(new_group_mapping)) %*% t(new_group_mapping) %*% new_weights^2
    
    if (any(abs(as.vector(new_weights_sq_W - new_weights_W^2)) > sqrt(.Machine$double.eps))){
      stop('weights are not constant within group (e.g. person) across tasks.')
    }
    std_weights <- 1/sum(new_weights_W) * nrow(new_W)
    new_weights_W <- new_weights_W * std_weights
    new_weights <- new_weights * std_weights
    if (abs(sum(new_weights_W) - nrow(new_W)) > sqrt(.Machine$double.eps)){
      warning('Weird Standardization Issue')
    }
    names(new_weights) <- NULL
    #Linear predictor and probability of voting "yes"
    new_xb <- new_X %*% coef(object)
    new_prob <- apply(new_xb, MARGIN = 2, plogis)
  }
  
  #Get the posterior predictive probability for each group 
  if (!by_group & !isTRUE(options$ignore_moderator)){
    if (inherits(object, 'FactorHet_refit')){
      coef_phi <- object$internal_parameters$phi
    }else{
      coef_phi <- coef(object, 'phi')
    }
    new_group_postpred <- softmax_matrix(new_W %*% t(coef_phi))
    #Project to each observation (probability of each observations being each group)
    if (!isTRUE(options$ignore_treatment)){
      new_postpred <- apply(new_group_postpred, MARGIN = 2, FUN=function(i){as.vector(new_group_mapping %*% i)})
    }else{
      new_postpred <- NULL
    }
  }else{
    new_group_postpred <- NA
  }
  if (return == 'postpred_only'){
    if (calc_gradient){
      attributes(new_group_postpred)$W <- new_W
      attributes(new_group_postpred)$unique_group <- unique_new_group
    }
    
    norm_weights <- as.vector(new_weights_W/sum(new_weights_W))
    
    attributes(new_group_postpred)$norm_weights <- norm_weights
    return(new_group_postpred)
  }

  if (!by_group & !isTRUE(options$ignore_moderator)){
    # Are the new groups in the training data? If so, get the posterior probabilities
    if (null_group){
      new_group_posterior <- matrix(NA, nrow = length(unique_new_group), ncol = ncol(coef(object)))    
    }else{
      if (!is.null(object$posterior)){
        new_group_posterior <- object$posterior$posterior[match(unique_new_group, object$posterior$posterior[,1]),-1,drop=F]
      }else{
        new_group_posterior <- matrix(NA, nrow = length(unique_new_group), ncol = ncol(coef(object)))    
      }
    }
    new_posterior <- apply(new_group_posterior, MARGIN = 2, FUN=function(i){as.vector(new_group_mapping %*% i)})
    
    if (nrow(new_group_mapping) == 1){
      new_prob_postpred <- sum(new_prob * new_postpred)
      new_prob_posterior <- sum(new_prob * new_posterior)
    }else{
      new_prob_postpred <- rowSums(new_prob * new_postpred)
      new_prob_posterior <- rowSums(new_prob * new_posterior)
    }
    
    if (type %in% c('posterior', 'both')){
      new_prob_agg <- new_prob_posterior
      #Replace with posterior probability when group is included
      new_prob_agg[which(is.na(new_prob_agg))] <- new_prob_postpred[which(is.na(new_prob_agg))]
    }else if (type == 'posterior_predictive'){
      new_prob_agg <- new_prob_postpred
    }
    
    if (type == 'both'){
      new_prob_agg <- cbind(new_prob_agg, new_prob_postpred)
    }
    
  }else{
    new_prob_agg <- NA
    new_postpred <- NA
  }
  
  
  if (return == 'prediction'){
    
    if (by_group){
      if (return_task){
        attributes(new_prob)$task <- new_task
        attributes(new_prob)$group <- new_group
      }
      if (calc_gradient){
        attributes(new_prob)$X <- new_X
      }
      return(new_prob)
    }else{
      
      new_prob_agg <- data.frame(prediction = new_prob_agg, 
                 task = new_task,
                 group = new_group, stringsAsFactors = F)
        
      if (!null_group){
        # Get unique groups in the data
        orig_data <- data.frame(
          group = prep_newdata$design[[prep_newdata$name_group]],
          stringsAsFactors = FALSE
        )
      }else{
        orig_data <- data.frame(
          group = new_group,
          stringsAsFactors = F
        )
      }

      if (!is.null(new_choice_order)){
        orig_data$choice_order <- prep_newdata$design[[prep_newdata$name_choice_order]]
      }
      if (!null_task){
        orig_data$task <- prep_newdata$design[[prep_newdata$name_task]]
        orig_data <- merge(orig_data, new_prob_agg, by = c('group', 'task'), sort = F)
      }else{
        orig_data <- merge(orig_data, new_prob_agg, by = c('group'), sort = F)
      }
      
      if (length(prep_newdata$name_group) != 0){
        stopifnot(all(orig_data$group == prep_newdata$design[[prep_newdata$name_group]]))
      }
      
      if (!null_task){
        stopifnot(all(orig_data$task == prep_newdata$design[[prep_newdata$name_task]]))
      }
      
      if (!is.null(new_choice_order)){
        orig_data$mod_prediction <- ifelse(orig_data$choice_order == object$internal_parameters$unique_choice[2], orig_data$prediction, 1 - orig_data$prediction)
      }else{
        orig_data$mod_prediction <- orig_data$prediction
      }
      
      output <- orig_data$mod_prediction
      return(output)
    }
  }else if (return == 'detailed'){
    return(list(
      X = new_X,
      prediction = new_prob_agg,
      posterior_predictive_group = new_group_postpred,
      posterior_predictive = new_postpred,
      prediction_by_group = new_prob,
      actual_posterior = new_posterior,
      group = new_group,
      unique_new_group = unique_new_group,
      outcome = new_y
    ))
  }else{stop('return must be prediction or detailed.')}
  
}