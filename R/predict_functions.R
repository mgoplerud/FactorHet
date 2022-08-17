#' Predict after using FactorHet
#' 
#' Get predictions of chosing particular profiles for either the data used
#' during estimation or some out-of-sample dataset.
#' 
#' @param object A model estimated using \code{FactorHet} or
#'   \code{FactorHet_mbo}.
#' @param newdata Data to estimate predictions on; default \code{NULL} uses
#'   estimation data.
#' @param type How to generate the prediction? If \code{"posterior"} (default),
#'   use the posterior probabilities for each observation in the training data
#'   for weighting the clusters. If \code{"posterior_predictive"}, use the
#'   cluster memberships implied by the moderator parameters, e.g. pi_{ik}. If
#'   using this for evaluating out-of-sample accuracy, using
#'   \code{"posterior_predictive"} may be sensible.
#' @param by_cluster Default of \code{FALSE}. If \code{TRUE}, return predicted
#'   probabilities by cluster.
#' @param return Default of \code{"prediction"} returns the prediction. Use
#'   \code{"detailed"} to return a variety of additional terms. This is mostly
#'   used internally for estimating marginal effects and other quantities of
#'   interest.
#' @param ... Miscellanous options used internally and not documented.
#' @examples 
#' data(immigration)
#' # Fit a model once for simplicity
#' fit_MBO <- FactorHet(Chosen_Immigrant ~ Plans + Ed + Country,
#'  design = immigration, lambda = 1e-4,
#'  initialize = FactorHet_init(nrep = 1),
#'  K = 2, group = ~ CaseID, task = ~ contest_no, 
#'  choice_order = ~ choice_id)
#' immigration$pred_FH <- predict(fit_MBO)  
#' @importFrom stats predict coef formula plogis delete.response update.formula
#'   na.omit qlogis
#' @export 
predict.FactorHet <- function(object, newdata = NULL, type = 'posterior', 
                              by_cluster = FALSE, 
                              return = 'prediction', ...){
  options <- list(...)
  
  type <- match.arg(type, c('posterior', 'posterior_predictive', 'both'))
  #Only relevant ... is "require_response" but not to be ever adjusted by user.
  if (length(setdiff(names(options), c('return_task', 'quick_exit', 
                                       'override_weights',
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
  
  prep_newdata <- prepare_formula(fmla_main = formula_object$het, 
    fmla_moderator = formula_object$mod, design = newdata, 
    group = formula_object$other_parameters$name_group, 
    task = formula_object$other_parameters$name_task, weights = formula_object$weights,
    choice_order = formula_object$other_parameters$name_choice_order ,
    delete_response = delete_response)
  
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
  if (!delete_response){
    new_y <- prep_newdata$design[[prep_newdata$outcome]]
  }else{
    new_y <- NULL
  }
  
  conjoint_names <- prep_newdata[c('name_group', 'name_task', 'name_choice_order')]
  factor_names <- prep_newdata$factor_names
  #Verify same factors
  stopifnot(identical(factor_names, names(object$internal_parameters$factor_levels)))
  
  if (is.null(new_group)){
    new_group <- 1:nrow(prep_newdata$design)
    null_group <- TRUE
  }else{
    null_group <- FALSE
  }
  unique_new_group <- unique(new_group)
  
  
  if (any(is.na(new_group))){stop('No Missingness Allowed for Group')}
  
  if (!by_cluster){#Get predictions averaged across, need info on moderators
    
    #Build the moderator variables
    if (object$K == 1){
      formula_object$mod <- ~ 1
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
      print('Non Identical Names in OOS W and W')
      print('OOS NAMES')
      print(colnames(new_W))
      print('MODEL NAMES')
      print(object$internal_parameters$W$args$names[-1])
      stop()
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
  new_group_mapping <- sparseMatrix(i = 1:nrow(new_X), j = match(new_group, unique_new_group), x = 1)
  
  new_group_mapping <- new_group_mapping
  new_weights <- new_weights
  
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
  
  
  #Get the posterior predictive probability for each group 
  if (!by_cluster){
    new_group_postpred <- softmax_matrix(new_W %*% t(coef(object, 'phi')))
    #Project to each observation (probability of each observations being each cluster)
    new_postpred <- apply(new_group_postpred, MARGIN = 2, FUN=function(i){as.vector(new_group_mapping %*% i)})
  }else{
    new_group_postpred <- NA
  }
  if (return == 'postpred_only'){
    if (calc_gradient){
      attributes(new_group_postpred)$W <- new_W
    }
    
    norm_weights <- as.vector(new_weights_W/sum(new_weights_W))
    
    attributes(new_group_postpred)$norm_weights <- norm_weights
    return(new_group_postpred)
  }

  if (!by_cluster){
    # Are the new groups in the training data? If so, get the posterior probabilities
    if (null_group){
      new_group_posterior <- matrix(NA, nrow = length(unique_new_group), ncol = ncol(coef(object)))    
    }else{
      new_group_posterior <- object$posterior$posterior[match(unique_new_group, object$posterior$posterior[,1]),-1,drop=F]
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
    
    if (by_cluster){
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
      prediction = new_prob_agg,
      posterior_predictive = new_postpred,
      prediction_by_cluster = new_prob,
      group = new_group,
      outcome = new_y
    ))
  }else{stop('return must be prediction or detailed.')}
  
}