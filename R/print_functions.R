#' Generic methods for FactorHet models
#' 
#' Brief descriptions of generic methods (e.g. print, summary) for FactorHet as
#' well as a way to visualize the progress of the model-based optimization.
#' 
#' @name FactorHet-class 
#' @param object Object fit using \code{\link{FactorHet}} or \code{\link{FactorHet_mbo}}.
#' @param x Object fit using \code{\link{FactorHet}} or \code{\link{FactorHet_mbo}}.
#' @param y Not used; required to maintain compatibility.
#' @param ... Optional arguments; only used by \code{plot.FactorHet} to pass
#'   arguments to \code{\link{cjoint_plot}}.
#' @param show_interactions Used by \code{summary.FactorHet}; indicates whether the
#'   interaction terms be shown. Default \code{FALSE}. See "Details" for more
#'   discussion.
#' @details The following methods with the arguments given above exist. All
#'   methods work on models with using \code{\link{FactorHet}} and
#'   \code{\link{FactorHet_mbo}}. 
#'   \describe{
#'   \item{plot: }{A shorthand for \code{\link{cjoint_plot}} on a fitted object.}
#'   \item{formula: }{Get the underlying formula for the treatment effects and
#'   moderators. This also returns the values used for \code{group},
#'   \code{task}, and \code{choice_order} if provided.}
#'   \item{print: }{Two print methods. For \code{\link{FactorHet}}, it summarizes the
#'   model and fusion of the factor levels. \code{fusion.tolerance} sets the
#'   threshold at which levels are reported as fused. For outputs of
#'   \code{\link{AME}} (and similar), this plots the corresponding plot. See
#'   that documentation for more details.}
#'   \item{summary: }{Summarizes the main effects by group with standard errors.
#'   It is typically more common to visualize this with \code{\link{cjoint_plot}} (and
#'   the accompanying data.frame) or \code{\link{AME}}.
#'   \code{show_interactions = TRUE} shows the interactions in addition to the
#'   main effects.}
#'   \item{coef: }{Returns the coefficient matrix on the original scale (i.e. with
#'   the sum-to-zero constraints). \code{code_type = "phi"} returns the
#'   moderator coefficients instead of the treatment effect coefficients.}
#'   \item{AIC and BIC: }{Reports the AIC or BIC. If multiple degrees of freedom
#'   options specified, returns a matrix. Otherwise, it returns a number.}
#'   \item{logLik: }{Returns the log-likelihood, log-posterior or sequence of
#'   log-posterior values at each iteration of the algorithm. See argument
#'   \code{type} above for details.}
#'   \item{visualize_MBO: }{For a model fit with \code{\link{FactorHet_mbo}}, show
#'   information about the MBO, i.e. proposed values and objectives.}
#'   \item{posterior_FactorHet: }{For a model with \code{K > 1}, visualize the
#'   posterior for each observation and the posterior predictive implied by the
#'   moderators.}
#'   \item{vcov.FactorHet}{Extracts the estimated variance-covariance 
#'   matrix of the parameters.}
#'  }
#' @export
plot.FactorHet <- function(x, y = NULL, ...){
  if (!is.null(y)){warning('y not used in plot.FactorHet')}
  suppressMessages(cjoint_plot(object = x, ...))
}

#' @rdname FactorHet-class
#' @importFrom stats as.formula
#' @export
formula.FactorHet <- function(x, ...){
  fmla <- lapply(x$internal_parameters$formula[c('het', 'mod')], as.formula, globalenv())
  
  if (x$internal_parameters$formula$weight == ""){
    fmla$weight <- NULL
  }else{
    fmla$weight <- as.formula(x$internal_parameters$formula$weight, globalenv())
  }
    
  fmla <- c(fmla, list('other_parameters' = x$internal_parameters$formula$other_parameters))
  return(fmla)
}

#' @rdname FactorHet-class 
#' @param x Model from FactorHet
#' @param fusion.tolerance Threshold at which to declare levels fused
#' @method print FactorHet
#' @export
print.FactorHet <- function(x, fusion.tolerance = 1e-3, ...){
  if (length(list(...)) > 0){stop('... not used for FactorHet print')}
  factor_levels <- x$internal_parameters$factor_levels
  beta <- x$parameters$beta
  K <- ncol(beta)
  J <- length(factor_levels)
  L_J <- lengths(factor_levels)
  all_combn <- sapply(L_J, choose, k = 2)
  ordered_factor <- x$internal_parameters$ordered_factor
  fused_levels <- x$internal_parameters$fusion
  fused_levels <- subset(fused_levels, fused_levels$main < fusion.tolerance & fused_levels$inter <= fusion.tolerance)
  
  if (nrow(fused_levels) > 0){
    
    any_fusion <- TRUE
    
    fusion_report <- list()
    for (k in 1:K){
      sub_k <- fused_levels[fused_levels$group == k, ]
      split_k <- with(sub_k, split(pair, factor))
      split_k <- mapply(split_k, names(split_k), SIMPLIFY = FALSE, FUN=function(f,j){
         o <- ordered_factor[j]
         if (!o & (length(f) == all_combn[j])){
           return('All levels fused.')
         }else if (o & (length(f) == (L_J[j] - 1))){
           return('All levels fused.')
         }else{
           return(paste(f, collapse = ', '))
         }
      })
      fusion_report[[k]] <- split_k
    }
    
  }else{
    any_fusion <- FALSE
  }
  cat('\n---Summary of FactorHet---\n')
  cat(paste0('The model was estimated with ', K, ' groups and ', J, ' factors with each group having ', nrow(beta), ' parameters.\n'))
  cat('\n')
  if (any(x$internal_parameters$ordered_factors)){
    ordered_factors <- x$internal_parameters$ordered_factors
    cat(paste0('The following factors are treated as ordered: ', 
      paste0(names(ordered_factors[which(ordered_factors)]), collapse = ', ')))  
    cat('\n\n')
  }
  
  cat('---Summary of Fusion---\n')
  if (any_fusion){
    for (k in 1:K){
      cat(paste0('Group ', k, ':'))
      fusion_k <- fusion_report[[k]]
      if (length(fusion_k) == 0){
        cat(' No levels fused.\n')
      }else{
        cat('\n')
        for (j in names(fusion_k)){
          cat(paste0(j, ' -- ', fusion_k[[j]]))
          cat('\n')
        }
          
      }
    }
  }else{
    cat('\n')
    cat('No levels were fused together.')
  }
  cat('\n----\n')
  cat('Use AME or AMIE to examine the effects in detail.\n')
}

#' @rdname FactorHet-class
#' @param digits Number of digits to include
#' @importFrom Matrix drop0
#' @method summary FactorHet
#' @export
summary.FactorHet <- function(object, show_interactions = FALSE, digits = 3, ...){
  if (length(list(...)) > 0){stop('... not used in FactorHet summary.')}  
  beta <- object$parameters$beta
  J <- length(object$internal_parameters$factor_levels)
  K <- ncol(beta)
  n_p <- nrow(beta)
  
  main_terms <- which(apply(object$internal_parameters$penalty$term_position, MARGIN = 1, FUN=function(i){sum(is.na(i))}) == J)
  all_se <- sqrt(diag(vcov.FactorHet(object, phi = FALSE)))
  if (!show_interactions){
    sum_beta <- (drop0(as.matrix(round(beta[main_terms,], digits))))
    sum_se <- (round(matrix(all_se, nrow = n_p)[main_terms,], digits))
  }else{
    sum_beta <- as.character(drop0(as.matrix(round(beta, digits))))
    sum_se <- (round(matrix(all_se, nrow = n_p), digits))
  }
  sum_names <- rownames(sum_beta)
  #https://stackoverflow.com/questions/36136808/make-all-elemants-of-a-character-vector-the-same-length
  sum_beta <- gsub("\\s", " ", format(as.vector(sum_beta), width=max(nchar(as.character(as.vector(sum_beta))))))
  sum_se <- paste0('(', sprintf(paste0('%.', digits, 'f'), sum_se), ')')
  sum_se <- gsub("\\s", " ", format(sum_se, width=max(nchar(sum_se))))
  
  sum_beta <- matrix(sum_beta, ncol = K)
  sum_se <- matrix(sum_se, ncol = K)
  sum_names <- gsub("\\s", " ", format(sum_names, width = max(nchar(sum_names))))

  sum_coef <- paste(sum_names, apply(sum_beta, MARGIN = 1, paste, collapse = ' '))
  sum_se <- paste(format(' ', width = max(nchar(sum_names))), apply(sum_se, MARGIN = 1, paste, collapse = ' '))
  
  sum_out <- c(rbind(sum_coef, sum_se))
  
  sum_out <- c(
    paste(format(' ', width = max(nchar(sum_names))), paste(format(paste0('Clus ', 1:K), width = max(nchar(sum_beta))), collapse = ' ')),
    sum_out
  )
  if (show_interactions){
    cat('Summary of All Effects\n')
  }else{
    cat('Summary of Main Effects\n')
  }
  cat(paste(sum_out, collapse='\n'))
  cat('\n')  
  
  data_output <- data.frame(est = as.vector(beta), se = all_se)
  data_output$variable <- rownames(data_output)
  rownames(data_output) <- NULL
  invisible(data_output)
}

#' @rdname FactorHet-class
#' @param coef_type Type of coefficient (beta for treatment effects; phi for moderators)
#' @export
coef.FactorHet <- function(object, coef_type = 'beta', ...){
  if (length(list(...)) > 0){stop('... not used in coef.FactorHet')}
  if (coef_type == 'beta'){
    as.matrix(object$parameters$beta)
  }else{
    if (inherits(object, 'FactorHet_refit')){
      stop('option not available for FactorHet_refit.')
    }
    as.matrix(object$parameters$phi)
  }
}

#' @rdname FactorHet-class
#' @importFrom stats logLik
#' @param type For "logLik", should the log-likelihood (\code{"loglik"}),
#'   log-posterior (\code{"log_posterior"}), or sequence of log-posterior values
#'   at each iteration (\code{"log_posterior_seq"}) be returned?
#' @export
logLik.FactorHet <- function(object, type = 'loglik', ...){
  
  if (inherits(object, 'FactorHet_refit')){
    stop('function not available for FactorHet_refit.')
  }
  
  if (length(list(...)) > 0){
    stop('... not used in logLik after FactorHet')
  }
  if (type == 'loglik'){
    lp_sequence <- object$internal_parameters$trajectory$ll[,2]
  }else if (type == 'log_posterior'){
    lp_sequence <- object$internal_parameters$trajectory$ll[,1]
  }else if (type == 'log_posterior_seq'){
    lp_sequence <- object$internal_parameters$trajectory$ll[,1]
    return(lp_sequence[!is.na(lp_sequence)])
  }else{stop('type for logLik must be "loglik" or "log_posterior" or "log_posterior_seq".')}
  lp_sequence <- lp_sequence[!is.na(lp_sequence)]
  return(tail(lp_sequence, 1))
}

#' @rdname FactorHet-class
#' @importFrom stats BIC
#' @export
BIC.FactorHet <- function(object, ...){
  ic <- object$information_criterion
  if (length(ic$BIC) == 1){
    return(ic$BIC)
  }else{
    return(ic[,c('method', 'BIC')])
  }
}
#' @rdname FactorHet-class
#' @export
AIC.FactorHet <- function(object, ...){
  ic <- object$information_criterion
  if (length(ic$AIC) == 1){
    return(ic$AIC)
  }else{
    return(ic[,c('method', 'AIC')])
  }
}

#' @rdname FactorHet-class
#' @method print FactorHet_vis
#' @export
print.FactorHet_vis <- function(x, ...){print(x$plot)}

#' @rdname FactorHet-class
#' @export
visualize_MBO <- function(object){
  
  if (inherits(object, 'FactorHet_refit')){
    stop('function not available for FactorHet_refit.')
  }
  if (!inherits(object, 'FactorHet_MBO')){
    stop('visualize_MBO can only be fit on an object fit FactorHet_mbo.')
  }
  
  .data <- NULL
  object$MBO_output$path$prop.type <- ifelse(object$MBO_output$path$prop.type == 'initdesign', 'Initial', 'MBO')
  object$MBO_output$path$power_ten <- 10^object$MBO_output$path$l
  show_curve <- ggplot(object$MBO_output$path) +
    geom_point(aes(x=.data[['power_ten']],y=.data[['y']],
                   col = .data[['prop.type']], 
                   pch = .data[['prop.type']])) +
    ylab('Criterion') +
    theme_bw() + scale_x_log10() +
    scale_color_manual(values = c('red', 'black')) +
    theme(legend.position = 'bottom') +
    labs(col = 'Type: ', pch = 'Type: ')
  
  final_criterion <- object$information_criterion[[object$MBO_output$control$criterion[1]]]
  if (object$MBO_output$control$mbo_type == 'sparse'){
    final_lambda <- object$parameters$orig_lambda
    lab_x <- bquote('(Unscaled) Regularization Parameter (' * lambda * ')')
  }else{
    final_lambda <- object$internal_parameters$control$prior_var_beta
    lab_x <- bquote('Variance of Normal Prior on ' * beta)
    
  }
  .data <- NULL
  show_curve <- show_curve +
    # Show the *final* estimated values
    geom_hline(aes(yintercept=final_criterion), linetype = 'dashed') +
    geom_vline(aes(xintercept=final_lambda), linetype = 'dashed') +
    xlab(lab_x)
  
  return(show_curve)
}


#' @rdname FactorHet-class
#' @export
posterior_FactorHet <- function(object){
  
  if (inherits(object, 'FactorHet_refit')){
    stop('function not available for FactorHet_refit.')
  }
  K <- ncol(coef.FactorHet(object))
  if (K == 1){stop('No posterior to visualize when K = 1.')}
  post.predict <- object$posterior$posterior_predictive
  posterior <- object$posterior$posterior
  
  compare_moderator <- data.frame()
  #If two groups, only show first group.
  if (K == 2){K <- 1}
  for (k in 1:K){
    compare_moderator <- rbind(compare_moderator, data.frame(group = post.predict[,1], post.predict = post.predict[,k+1], posterior = posterior[,k+1], K = k))
  }
  if (K == 1){
    compare_baseline <- data.frame(K = 1:K, baseline = object$parameters$pi[1])
    compare_posterior_mode <- data.frame(K = 1:K, baseline = posterior$group_1)
  }else{
    compare_baseline <- data.frame(K = 1:K, baseline = object$parameters$pi)
    compare_posterior_mode <- data.frame(K = 1:K, baseline = colMeans(posterior[,-1]))
  }
  
  vis_moderator <- do.call('rbind', lapply(setdiff(names(compare_moderator), c('group', 'K')),
                                           FUN=function(i){
                                             out <- compare_moderator[,c('group', 'K', i)]
                                             names(out)[ncol(out)] <- 'value'
                                             out$variable <- as.character(i)
                                             return(out)
                                           }))
  vis_moderator <- vis_moderator[, c('group', 'K', 'variable', 'value')]
  
  vis_moderator$variable <- ifelse(vis_moderator$variable == 'posterior', 'Posterior', 'Posterior Predictive')
  vis_moderator$variable <- factor(vis_moderator$variable, levels = c('Posterior Predictive', 'Posterior'))
  vis_moderator <- ggplot(vis_moderator) + 
    geom_density(aes(x=.data[['value']], 
                     group=.data[['variable']], 
                     fill = .data[['variable']], 
                     col = .data[['variable']]), alpha = 0.10) + 
    facet_wrap(~paste0('Group ', K)) + 
    geom_vline(aes(xintercept=.data[['baseline']]), data = compare_baseline, linetype = 'dashed') +
    theme_bw() + xlab('Probability') + ylab('Density') +
    scale_color_manual(values = c('black', 'red')) +
    scale_fill_manual(values = c('grey', 'red')) +
    theme(legend.position = 'bottom') +
    labs(col = '', fill = '')
  print(vis_moderator)
  output <- list(plot = vis_moderator, compare = compare_moderator, baseline = compare_baseline)
  invisible(output)
}

