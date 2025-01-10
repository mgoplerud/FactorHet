#' Estimate heterogeneous treatment effects by individual or moderator
#' @name HTE
#' @description These functions estimate heterogeneous effects from
#'   \code{\link{FactorHet}} at the individual level or by the average value of
#'   a moderator. They can be used to produce individual-level estimates that
#'   can be compared against other methods for estimating heterogeneous effects.
#' @details The functions here allow for, first, estimation of conditional
#'   average marginal effects for each individual given their pre-treatment
#'   moderators (\code{HTE_by_individual}). This is a weighted average of the
#'   AME for each group by the individual's group membership probabilities, i.e.
#'   \eqn{\hat{\pi}(X_i)} (Goplerud et al. 2025). These are also averaged
#'   together to return an estimate to produce a population-level effect.
#'   
#'   Second, one can estimate conditional average marginal effects using
#'   \code{HTE_by_moderator}. This takes a moderator such as party
#'   identification and counterfactually sets each observation to some level
#'   (e.g., "Strong Democrat"). It then reports the average of the
#'   individual-level conditional effects across the sample population as the
#'   "conditional" average marginal effect. If \code{overall_AME} is true, it
#'   also returns the average of the individual heterogeneous effects given the
#'   observed distribution of pre-treatment moderators. It and the
#'   \code{population} element of the list produced by \code{HTE_by_individual}
#'   coincide exactly.
#'   
#'   Both functions can be used with split-sample or refit, i.e.
#'   \code{\link{FactorHet_refit}}, and the computed AME, although this will not
#'   take into account uncertainty in the moderator estimation as they are
#'   assumed fixed when refitting the model.
#'   
#'   To use these functions, first estimate the AMEs by group, i.e., using
#'   \code{\link{AME}} and then pass this object and the original
#'   \code{\link{FactorHet}} model to the functions for computing heterogeneous
#'   effects by moderator or individual.
#' 
#' @param object An object from \code{\link{FactorHet}},
#'   \code{\link{FactorHet_mbo}}.
#' @param AME An estimate of the average marginal effects by group from
#'   \code{\link{AME}}.
#' @param design An optional data.frame of moderator values on which to produce
#'   the individual-level or average conditional average marginal effects.
#'   \bold{Note:} There should be one row per observation if this function is
#'   used. The default is \code{NULL} which uses the estimation data.
#' 
#' @return \code{HTE_by_individual} returns a list with two data.frames. The
#'   first \code{individual} contains the estimated individual conditional
#'   average marginal effects. The second \code{population} contains the average
#'   of those individual effects. Standard errors (via the column \code{var})
#'   are also included.
#'   
#'   \code{HTE_by_population} returns a list for each moderator that consists
#'   itself of a list of each value of the moderator used. The value
#'   \code{"out"} contains the conditional average marginal effects.
#' @examples 
#' \donttest{
#' data(immigration)
#' # Estimate model with arbitrary choice of lambda
#' fit <- FactorHet(Chosen_Immigrant ~ Plans + Ed + Country,
#'   design = immigration, lambda = 1e-2,
#'   moderator = ~ party_ID,
#'   K = 2, group = ~ CaseID,
#'   control = FactorHet_control(init = 'mclust'),
#'   task = ~ contest_no, choice_order = ~ choice_id)
#' # Estimate AME 
#' est_AME <- AME(fit)
#' # Get individual AME; note only seven distinct
#' # patterns will exist as partyID is the only (discrete)
#' # moderator
#' iAME <- HTE_by_individual(
#'   object = fit,
#'   AME = est_AME)
#' # Get conditional AME by level of party ID
#' cAME_pID <- HTE_by_moderator(
#'   object = fit,
#'   AME = est_AME, overall_AME = TRUE)
#' 
#' AME_1 <- cAME_pID$`Overall AME`$out[[1]][,c('factor', 'level', 'est', 'var')]
#' AME_2 <- iAME$population[,c('factor', 'level', 'est', 'var')]
#' rownames(AME_1) <- rownames(AME_2) <- NULL
#' stopifnot(isTRUE(all.equal(AME_1, AME_2)))
#' }
#' @importFrom stats vcov
#' @export
HTE_by_individual <- function(object, AME, design = NULL){
  
  K <- object$K
  if (inherits(object, 'FactorHet_refit')){
    warning('Uncertainty in phi ignored for refit model')
    n_W <- length(object$internal_parameters$W$args$names)
    full_vcov_object <- bdiag(vcov(object, phi = FALSE), Diagonal(x = rep(0,n_W * (K-1))))
  }else{
    full_vcov_object <- vcov(object)
  }
  AME_estimates <- AME$data[AME$data$baseline == FALSE,]
  AME_gradients <- AME$gradient
  unique_combinations <- unique(AME_estimates[, c('factor', 'level')])
  seq_K <- 1:K
  format_AME <- lapply(1:nrow(unique_combinations), FUN=function(v){
    combo_v <- unique_combinations[v,]
    pos_v <- which((AME_estimates$factor == combo_v$factor) & (AME_estimates$level == combo_v$level))
    mfx_v <- AME_estimates[pos_v, ]  
    est_v <- mfx_v$marginal_effect
    var_v <- mfx_v$var
    group_v <- as.numeric(mfx_v$group)
    grad_v <- AME_gradients[pos_v]
    grad_v <- do.call('cbind', grad_v)
    stopifnot(isTRUE(all.equal(group_v, seq_K)))
    return(list(est = est_v, grad = grad_v))
  })
  if (is.null(design)){
    # Get the model frame for the original set of moderators
    design <- object$internal_parameters$data$design
    design <- clean_design_for_HTE(object, design)
  }
  
  pred_design <- predict(object,
     newdata = design,
     calc_gradient = TRUE,
     return = 'postpred_only',
     ignore_treatment = TRUE,
     override_weights = FALSE)
  W <- attributes(pred_design)$W
  if (nrow(W) != nrow(design)){stop('W and design sizes do not line up.')}
  
  grad_pi_ik <- est_grad_pi_ik(K = K, 
     prob = pred_design, W = W, 
     weights = NULL, individual = TRUE)

  counter <- 0
  individual_HTE <- mapply(format_AME, unique_combinations$level, 
   unique_combinations$factor, SIMPLIFY = FALSE, FUN=function(data_jl, level_jl, factor_jl){
     counter <<- counter + 1
     term_1 <- Reduce('+', mapply(data_jl$est, 
        grad_pi_ik, SIMPLIFY = FALSE, 
        FUN=function(i,j){i * j}))
     term_2 <- do.call('cbind', lapply(1:K, FUN=function(k){
       outer(pred_design[,k], data_jl$grad[,k], '*')
     }))
     grad_all <- cbind(term_2, term_1)
     var_all <- rowSums( (grad_all %*% full_vcov_object) * grad_all )
     est_all <- as.vector(pred_design %*% data_jl$est)
     out_all <- data.frame(raw_id = 1:length(est_all), est = est_all, var = var_all, factor = factor_jl, level = level_jl, stringsAsFactors = FALSE)
     out_all <- out_all[order(est_all),]
     out_all$order_id <- 1:nrow(out_all)
     out_all$counter <- counter
     return(list(out = out_all, grad = grad_all))
  })
  individual_grad <- sapply(individual_HTE, `[[`, 'grad')
  # Get the average of the gradients
  grad_pHTE <- lapply(individual_grad, colMeans)
  grad_pHTE <- do.call('rbind', grad_pHTE)
  est_pHTE <- sapply(individual_HTE, FUN=function(i){mean(i$out$est)})
  var_pHTE <- rowSums( (grad_pHTE %*% full_vcov_object) * grad_pHTE)
  population_HTE <- cbind(unique_combinations, data.frame(est = est_pHTE, var = var_pHTE))

  individual_HTE <- do.call('rbind', lapply(individual_HTE, `[[`, 'out'))
  if ('...group' %in% names(design)){
    individual_HTE$group <- design[['...group']][individual_HTE$raw_id]
  }else{
    individual_HTE$group <- attributes(pred_design)$unique_group[individual_HTE$raw_id]
  }
  
  return(list(
    individual = individual_HTE,
    population = population_HTE
  ))
  
}

#' @rdname HTE
#' @param moderators An argument that contains a list of moderators to evaluate.
#'   The default is \code{NULL} and considers all moderators.
#' @param points_continuous A positive integer value that indicates the number
#'   of equally spaced points to evaluate a continuous moderator over.
#' @param overall_AME A logical value that indicates whether to compute the AME
#'   over the entire \code{design} without modification. The default is
#'   \code{FALSE}.
#' @param verbose A logical value that indicates whether progress should be
#'   reported. The default is \code{FALSE}.
#' @importFrom stats vcov
#' @export
HTE_by_moderator <- function(object, AME,
    moderators = NULL, design = NULL,
    points_continuous = 10, overall_AME = FALSE,
    verbose = FALSE){
  
  K <- object$K
  if (inherits(object, 'FactorHet_refit')){
    warning('Uncertainty in phi ignored for refit model')
    n_W <- length(object$internal_parameters$W$args$names)
    full_vcov_object <- bdiag(vcov(object, phi = FALSE), Diagonal(x = rep(0, n_W * (K-1))))
  }else{
    full_vcov_object <- vcov(object)
  }
  
  AME_estimates <- AME$data[AME$data$baseline == FALSE,]
  AME_gradients <- AME$gradient
  unique_combinations <- unique(AME_estimates[, c('factor', 'level')])
  seq_K <- 1:K
  format_AME <- lapply(1:nrow(unique_combinations), FUN=function(v){
    combo_v <- unique_combinations[v,]
    pos_v <- which((AME_estimates$factor == combo_v$factor) & (AME_estimates$level == combo_v$level))
    mfx_v <- AME_estimates[pos_v, ]  
    est_v <- mfx_v$marginal_effect
    var_v <- mfx_v$var
    group_v <- as.numeric(mfx_v$group)
    grad_v <- AME_gradients[pos_v]
    grad_v <- do.call('cbind', grad_v)
    stopifnot(isTRUE(all.equal(group_v, seq_K)))
    return(list(est = est_v, grad = grad_v))
  })
  
  if (is.null(design)){
    # Get the model frame for the original set of moderators
    design <- object$internal_parameters$data$design
    design <- clean_design_for_HTE(object, design)
  }
  raw_mf <- model.frame(formula(object)$mod, design)
  
  factors <- object$internal_parameters$factor_levels

  if (!is.null(moderators)){
    # Select only those of interest
    mf <- raw_mf[, moderators, drop = F]
  }else{
    mf <- raw_mf
  }
  
  unique_mf <- lapply(1:ncol(mf), FUN=function(i){
    mf_i <- mf[[i]]
    if (class(mf_i) %in% c('factor', 'character')){
      u_i <- unique(mf_i)
    }else{
      u_i <- tryCatch(seq(min(mf_i), max(mf_i), length.out=points_continuous), error = function(e){NULL})
      if (is.null(u_i)){
        warning(paste0('Moderator ', colnames(mf)[i], ' seems to be neither factor nor character. Using unique values.'))
        message(paste0('Moderator ', colnames(mf)[i], ' seems to be neither factor nor character. Using unique values.'))
        u_i <- unique(mf_i)
      }
    }
    return(u_i)
  })
  names(unique_mf) <- colnames(mf)
  if (overall_AME){
    unique_mf <- c(unique_mf, list(`Overall AME` = NA))
  }
  store_moderator <- list()
  if (verbose){
    message('Looping over moderators')
  }
  for (v in names(unique_mf)){
    if (verbose){message(v)}
    copy_design <- design
    # Loop over value of moderators
    u_mf <- unique_mf[[v]]
    counter_size <- length(u_mf)
    copy_all <- copy_all_prob <- copy_all_grad <- as.list(rep(NA, counter_size))
    for (counter in 1:counter_size){
      value <- u_mf[counter]
      if (is.na(value) & v == 'Overall AME'){
        copy_design <- design
      }else{
        copy_design <- design
        copy_design[,v] <- value
      }
      pred_copy <- predict(object,
         newdata = copy_design,
         calc_gradient = TRUE,
         return = 'postpred_only', 
         ignore_treatment = TRUE,
         override_weights = FALSE)
      copy_W <- attributes(pred_copy)$W
      norm_weights <- attributes(pred_copy)$norm_weights
      if (is.null(norm_weights) || length(norm_weights) == 0){
        mean_copy <- colMeans(pred_copy)
        norm_weights <- rep(1/nrow(pred_copy), nrow(pred_copy))
        mean_copy_2 <- colSums(Diagonal(x = norm_weights) %*% pred_copy)
        stopifnot(isTRUE(all.equal(mean_copy, mean_copy_2)))
      }else{
        mean_copy <- colSums(Diagonal(x = norm_weights) %*% pred_copy)
      }
      grad_pi_ik <- est_grad_pi_ik(K = K, 
       prob = pred_copy, W = copy_W, 
       weights = norm_weights)
      copy_AME <- mapply(format_AME, unique_combinations$level, 
        unique_combinations$factor, SIMPLIFY = FALSE, FUN=function(data_jl, level_jl, factor_jl){
        term_1 <- Reduce('+', mapply(data_jl$est, 
                     grad_pi_ik, SIMPLIFY = FALSE, 
                     FUN=function(i,j){i * j}))
        term_2 <- data_jl$grad %*% Diagonal(x = mean_copy)
        grad_all <- c(as.vector(term_2), as.vector(term_1))
        var_all <- as.numeric(t(grad_all) %*% full_vcov_object %*% grad_all)
        est_all <- sum(data_jl$est * mean_copy)
        out_all <- data.frame(est = est_all, var = var_all, factor = factor_jl, level = level_jl, stringsAsFactors = FALSE)
        return(list(out = out_all, grad = grad_all))
      })
      copy_grad <- sapply(copy_AME, `[[`, 'grad')
      copy_AME <- do.call('rbind', lapply(copy_AME, `[[`, 'out'))
      copy_AME$moderator <- v
      copy_AME$value <- value
      copy_prob <- data.frame(group = 1:K, prob = mean_copy, moderator = v, value = value, stringsAsFactors = FALSE)
      copy_all[[counter]] <- copy_AME
      copy_all_grad[[counter]] <- copy_grad
      copy_all_prob[[counter]] <- copy_prob
    }
    store_moderator[[v]] <- list(out = copy_all, grad = copy_all_grad, prob = copy_all_prob, value = u_mf)
  } 
  return(store_moderator)
}


est_grad_pi_ik <- function(K, prob, W, vcov, weights, individual=FALSE){
  
  if (individual){
    
    delta_out <- lapply(1:K, FUN=function(k){
      grad <- do.call('cbind', lapply(2:K, FUN=function(l){
        g <- Diagonal(x = prob[,k] * ((l == k) - prob[,l])) %*% W
        return(g)
      }))
      return(grad)
    })
    
  }else{

    delta_out <- lapply(1:K, FUN=function(k){
      grad <- sapply(2:K, FUN=function(l){
        g <- colSums(Diagonal(x = weights) %*% Diagonal(x = prob[,k] * ((l == k) - prob[,l])) %*% W)
        return(g)
      })
      return(grad)
    })
    
  }
  return(delta_out)
}

clean_design_for_HTE <- function(object, design){

  formula_object <- formula(object)
  design <- prepare_formula(fmla_main = formula_object$het, 
                            fmla_moderator = formula_object$mod, design = design, 
                            group = formula_object$other_parameters$name_group, 
                            task = formula_object$other_parameters$name_task, 
                            weights = formula_object$weights,
                            choice_order = formula_object$other_parameters$name_choice_order ,
                            delete_response = TRUE)
  
  new_group <- data.frame(design$design[design$name_group], stringsAsFactors = F)
  
  if (ncol(new_group) == 0){
    new_group <- NULL
  }else{
    new_group <- as.character(as.vector(new_group[,1]))
  }
  
  if (is.null(new_group)){
    new_group <- 1:nrow(design$design)
    null_group <- TRUE
  }else{
    null_group <- FALSE
  }
  unique_new_group <- unique(new_group)
  
  design <- model.frame(formula_object$mod, design$design)
  design[['...group']] <- new_group
  design <- unique(design)
  if (anyDuplicated(design[['...group']]) != 0){
    stop('moderators are not unique by ~ group')
  }
  return(design)  
}
