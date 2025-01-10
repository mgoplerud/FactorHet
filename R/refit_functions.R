#' @importFrom stats glm binomial vcov
refit_model_FH <- function(object, tolerance, data, hard_assign,
                           iter_refit = 50){
  
  recons.beta <- coef(object)
  K <- ncol(recons.beta)
  single_intercept <- object$internal_parameters$single_intercept
  factor_levels <- object$internal_parameters$factor_levels
  term_position <- object$internal_parameters$refit$term_position
  coef_names <- object$internal_parameters$refit$coef_names
  make_X_refit <- object$internal_parameters$refit$make_X_refit
  clip_tiny <- object$internal_parameters$misc$clip_tiny
  ordered_factors <- object$internal_parameters$ordered_factors
  # Get the restrictions that are nearly binding
  fusion_analysis <- mapply(list_from_cols(recons.beta), 1:ncol(recons.beta), SIMPLIFY = FALSE, FUN=function(b, cl){
    out <- prepare_fusion(factor_levels = factor_levels, term_position = term_position, 
                          coef_names = coef_names, ordered_factors = ordered_factors, beta = b, simplify = FALSE)
    out <- lapply(out, FUN=function(i){
      sapply(i, FUN=function(j){max(unlist(j))})
    })
    out <- lapply(out, FUN=function(j){which(j <= tolerance)})
    return(out)
  })

  if (object$internal_parameters$control$log_method != 'standard'){
    Fmatrix <- object$internal_parameters$refit$Fmatrix_orig
    basis_M <- object$internal_parameters$refit$make_X_refit$add_args$M
    Fmatrix <- lapply(Fmatrix, FUN=function(j){lapply(j, FUN=function(l){t(basis_M) %*% l %*% basis_M})})
    p_X <- ncol(basis_M)
  }else{
    Fmatrix <- object$internal_parameters$Fmatrix
    p_X <- object$internal_parameters$refit$p_X
    basis_M <- object$internal_parameters$basis_M

  }
  # Calculate the appropriate corresponding nullspace
  if (single_intercept){
    binding_null_basis <- lapply(fusion_analysis, FUN=function(r_k){
      binding_k <- do.call('rbind', mapply(r_k, Fmatrix, SIMPLIFY = FALSE, FUN=function(i, F_j){
        subset_Fj <- F_j[i]
        if (!identical(names(subset_Fj), names(i))){stop('name misalignment')}
        return(do.call('rbind', subset_Fj))
      }))
      binding_k <- binding_k[,-p_X]
      if (is.null(binding_k)){
        return(sparse_diag(rep(1, p_X - 1)))
      }else{
        drop0(zapsmall(calculate_nullspace_basis(binding_k), clip_tiny))
      }
    })
    binding_null_basis <- bdiag(c(1, bdiag(binding_null_basis)))
  }else{
    binding_null_basis <- lapply(fusion_analysis, FUN=function(r_k){
      binding_k <- do.call('rbind', mapply(r_k, Fmatrix, SIMPLIFY = FALSE, FUN=function(i, F_j){
        subset_Fj <- F_j[i]
        if (!identical(names(subset_Fj), names(i))){stop('name misalignment')}
        return(do.call('rbind', subset_Fj))
      }))
      if (is.null(binding_k)){
        return(sparse_diag(rep(1, p_X)))
      }else{
        drop0(zapsmall(calculate_nullspace_basis(binding_k), clip_tiny))
      }
    })
    binding_null_basis <- bdiag(binding_null_basis)
  }
  # Generate the probabilities of group membership
  
  pred_new <- predict(object, newdata = data, return = 'detailed', require_response = TRUE)
  
  y <- pred_new$outcome
  oX <- pred_new$X
  X <- do.call(make_X_refit$add_col, c(list(X = oX), make_X_refit$add_args))
  obs.E.prob <- pred_new$posterior_predictive 
  
  
  stopifnot(nrow(obs.E.prob) == nrow(X))
  if (hard_assign){
    obs.E.prob <- as.matrix(sparseMatrix(i = 1:nrow(obs.E.prob), 
      j = apply(obs.E.prob, MARGIN = 1, which.max), 
      x = 1, dims = c(nrow(obs.E.prob), K)))
  }
  # Generate the data
  if (single_intercept){
    blocked_X <- cbind(1, kronecker(diag(K), X[,-ncol(X)]))
  }else{
    blocked_X <- kronecker(diag(K), X)
  }
  rX <- as.matrix(blocked_X %*% binding_null_basis)
  if (iter_refit > 0){
    G <- sparseMatrix(i = 1:nrow(X), 
                      j = match(pred_new$group, pred_new$unique_new_group), 
                      x = 1, dims = c(nrow(X), length(pred_new$unique_new_group)))
    refit_object <- simple_EM_refit(y = y, X = rX, K = K, G = G,
        prior = pred_new$posterior_predictive_group, 
        iter = iter_refit)
    refit_model <- refit_object$model
    refit_posterior <- refit_object$posterior
  }else{
    refit_posterior <- NULL
    # refit_model <- suppressWarnings(
    #   arm::bayesglm(rep(y, K) ~ 0 + rX,
    #            weights = as.vector(obs.E.prob), 
    #            prior.scale = 1, prior.df = Inf,
    #            scaled = FALSE,
    #            family = binomial())
    # )
    refit_model <- suppressWarnings(
      glm(rep(y, K) ~ 0 + rX,
          weights = as.vector(obs.E.prob), family = binomial())
    )
  }
  
  if (iter_refit > 0){
    vcov_refit_raw <- refit_object$vcov
  }else{
    vcov_refit_raw <- vcov(refit_model)
  }
  # Estimate the refit model
  # Adjust the coefficients
  if (single_intercept){
    refit_beta <- raw_beta <- as.vector(binding_null_basis %*% coef(refit_model))
    refit_mu <- coef(refit_model)[1]
    
    refit_beta <- matrix(refit_beta[-1], ncol = K)
    refit_beta <- rbind(refit_beta, refit_mu)
    
    vcov_refit_beta <- as.matrix(
      binding_null_basis %*% vcov_refit_raw %*% t(binding_null_basis)
    )  
    # Permutation matrix
    P <- do.call('rbind', lapply(1:K, FUN=function(k){
      rbind(
        cbind((k-1) * p_X + 1:(p_X - 1), 1 + (k-1) * (p_X - 1) + 1:(p_X - 1)),
        cbind(k * p_X, 1)
      )
    }))
    P <- sparseMatrix(i = P[,1], j = P[,2], x = 1)
    stopifnot(all.equal(as.vector(P %*% raw_beta), as.vector(refit_beta))) 
    vcov_refit_beta <- P %*% vcov_refit_beta %*% t(P)
  }else{
    refit_beta <- as.vector(binding_null_basis %*% coef(refit_model))
    refit_beta <- matrix(refit_beta, ncol = K)
    
    vcov_refit_beta <- as.matrix(
      binding_null_basis %*% vcov_refit_raw %*% t(binding_null_basis)
    )
  }
  
  Pout <- kronecker(Diagonal(n = K), basis_M[1:length(coef_names),]) 
  vcov_refit_beta <- Pout %*% vcov_refit_beta %*% t(Pout)
  # Project back to the original space
  recons.refit_beta <- (basis_M %*% refit_beta)[1:length(coef_names),,drop=F]
  rownames(recons.refit_beta) <- rownames(recons.beta)
  
  names_vcov_refit <- do.call('c', lapply(1:K, FUN=function(k){paste0('beta', k, '_', rownames(recons.beta))}))
  rownames(vcov_refit_beta) <- colnames(vcov_refit_beta) <- names_vcov_refit

  # Get the fusion analysis on the refit model
  fusion_analysis <- mapply(list_from_cols(recons.beta), 1:ncol(recons.beta), SIMPLIFY = FALSE, FUN=function(b, cl){
    out <- prepare_fusion(factor_levels = factor_levels, term_position = term_position, 
                          coef_names = coef_names, ordered_factors = ordered_factors, beta = b, simplify = FALSE)
    out <- lapply(out, FUN=function(i){
      sapply(i, FUN=function(j){max(unlist(j))})
    })
    out <- lapply(out, FUN=function(j){which(j <= tolerance)})
    return(out)
  })
  
  out <- list(est = recons.refit_beta, 
              vcov = vcov_refit_beta, 
              refit_posterior = refit_posterior,
              fusion = fusion_analysis)
  return(out)
}

#' @title Refit model using estimated sparsity patterns
#' @description Using a previously estimated model, this function takes the
#'   estimated sparsity patterns (e.g., which levels are fused together) and the
#'   estimates of the moderator parameters, \eqn{\hat{\phi}}, and re-estimates the
#'   regression parameters \eqn{\beta}.
#' @details The main use of this function is to enable sample-splitting as
#'   discussed in Goplerud et al. (2025) to improve coverage and remove bias
#'   from the initial estimates. An example is provided below.
#' @param object An object from \code{\link{FactorHet}} or
#'   \code{\link{FactorHet_mbo}}.
#' @param newdata A data.frame containing the data to be estimated in the refit
#'   model.
#' @param tolerance A numerical value that sets the threshold at which to
#'   declare two levels as "fused"; the default is 1e-3. Two levels meet this
#'   threshold if the maximum difference between the main effects and any
#'   interactions is \code{tolerance}.
#' @param hard_assign A logical value that sets whether observations should be
#'   be assigned to the most probable cluster given \eqn{\phi} from the original
#'   model or whether they should they be weighted according to their estimated
#'   group membership probabilities, \eqn{\hat{\pi}(X_i)}. The default is
#'   \code{FALSE} which uses the weighted method.
#' @param iter_refit An integer value that sets the number of iterations used in
#'   fitting the refit model. The default is 200. A warning will be produced if
#'   it does not converge in this many iterations.
#' @examples
#' data(immigration)
#' set.seed(1)
#' # Split the data into two parts for sample-splitting
#' train_data <- subset(immigration, CaseID < 900)
#' refit_data <- subset(immigration, CaseID >= 900)
#' # Fit using fixed lambda for demonstration
#' # only
#' fit <- FactorHet(Chosen_Immigrant ~ Plans + Ed + Country,
#'   design = train_data, lambda = 1e-2,
#'   moderator = ~ party_ID + census_div,
#'   control = FactorHet_control(init = 'mclust'),
#'   K = 2, group = ~ CaseID, task = ~ contest_no, choice_order = ~ choice_id)
#' # Refit using the other half of data
#' \donttest{
#' refit <- FactorHet_refit(fit, newdat = refit_data)
#' # AME (etc.) for treatment effects can be computed as normal
#' AME_refit <- AME(refit)
#' # As can heterogeneous effects, although uncertainty in 
#' # phi is ignored
#' HTE_refit <- HTE_by_individual(refit, AME_refit, design = immigration)
#' }
#' @importFrom stats vcov
#' @return An object of class \code{\link{FactorHet}} that contains the output
#'   described the linked documentation.
#' @export
FactorHet_refit <- function(object, newdata, tolerance = 1e-3, 
    hard_assign = FALSE, iter_refit = 200){
  
  if (missing(newdata)){stop('newdata cannot be missing for FactorHet_refit')}
  # Refit the model
  refit_coef <- refit_model_FH(object = object, 
    tolerance = tolerance, data = newdata, 
    hard_assign = hard_assign, iter_refit = iter_refit)
  refit_object <- object
  refit_object$parameters <- list(beta = refit_coef$est)
  stopifnot(all.equal(coef(refit_object), as.matrix(refit_coef$est)))
  refit_object$vcov <- list(vcov = refit_coef$vcov)
  stopifnot(all.equal(vcov(refit_object), as.matrix(refit_coef$vcov)))
  
  refit_internal <- object$internal_parameters[c('factor_levels', 
                                                 'ordered_factors',
      'formula', 'penalty', 'rare', 'W', 'unique_choice', 'fusion')]
  refit_object <- refit_object[c('parameters', 'vcov', 'K')]
  refit_object$internal_parameters <- refit_internal
  # assign original phi parameters for use in other predictions
  refit_object$internal_parameters$phi <- object$parameters$phi
  refit_object$internal_parameters$control <- list(
    forced_randomize = object$internal_parameters$control$forced_randomize
  )
  refit_object$internal_parameters$data <- list(design = newdata)
  class(refit_object) <- c('FactorHet_refit', 'FactorHet')
  return(refit_object)
}

internal_align <- function(m1, m2, f = function(x){mean(abs(x))}, return_error = FALSE){
  # https://stackoverflow.com/questions/11095992/generating-all-distinct-permutations-of-a-list-in-r  
  getPerms <- function(x) {
    if (length(x) == 1) {
      return(x)
    }
    else {
      res <- matrix(nrow = 0, ncol = length(x))
      for (i in seq_along(x)) {
        res <- rbind(res, cbind(x[i], Recall(x[-i])))
      }
      return(res)
    }
  }
  orderings <- getPerms(1:ncol(m1))
  error_ordering <- apply(orderings, MARGIN = 1, FUN=function(i){
    f(m1 - m2[,i])
  })
  if (return_error){
    return(error_ordering)
  }else{
    return(orderings[which.min(error_ordering),])
  }
}

#' @importFrom stats glm binomial vcov
simple_EM_refit <- function(y, X, K, G, prior, iter){
  
  rep_y <- rep(y, K)
  old_posterior <- posterior <- prior
  
  prior_coef <- NULL
  for (it in 1:iter){
    # print(it)
    posterior_obs <- G %*% posterior
    # fit <- suppressWarnings(
    #   arm::bayesglm(rep_y ~ 0 + as.matrix(X), 
    #                 family = binomial(),
    #                 start = prior_coef,
    #                 prior.scale = 1, prior.df = Inf,
    #                 scaled = FALSE,
    #                 # control = list(trace = TRUE),
    #                 weights = as.vector(posterior_obs))
    # )
    fit <- suppressWarnings(glm(rep_y ~ 0 + as.matrix(X), 
      family = binomial(),
      # control = list(trace = TRUE),
      weights = as.vector(posterior_obs))
    )
    prior_coef <- coef(fit)
    ll <- plogis(as.vector(X %*% coef(fit)) * (2 * rep_y - 1), log.p = TRUE)
    # ll <- plogis(as.vector(X %*% coef(fit)))
    # ll <- ifelse(rep(y,K) == 1, log(ll), log(1-ll))
    ll <- matrix(ll, ncol = K)
    ll <- apply(ll, MARGIN = 2, FUN=function(i){as.vector(i %*% G)})
    posterior <- softmax_matrix(ll + log(prior))
    change <- max(abs(old_posterior - posterior))
    if (change < 1e-5){
      # print(paste0('Converged at ', it))
      break
    }
    old_posterior <- posterior
  }
  if (it == iter){message('Refitting mixture did not converge.')}
  
  # Calculate the correct standard errors
  # using Louis (1982)
  split_id <- split(1:nrow(X), rep(1:K, each = nrow(G)))
  beta <- coef(fit)
  p <- plogis(as.vector(X %*% beta))
  HC <- -t(X) %*% Diagonal(x = as.vector(G %*% posterior) * p * (1 - p)) %*% X
  
  # if (regularize){
  #   HC <- HC + - Diagonal(x = rep(1, ncol(X)))
  # }
  
  weight_X <- Diagonal(x = (rep_y - p)) %*% X
  weight_X <- lapply(split_id, FUN=function(i){
    t(G) %*% weight_X[i,]
  })
  
  var_score <- Reduce("+", lapply(1:K, FUN=function(k){
    Reduce("+", lapply(1:K, FUN=function(l){
      mod_weight <- posterior[,k] * ( (k == l) - posterior[,l])
      out <- t(weight_X[[k]]) %*% Diagonal(x = mod_weight) %*% weight_X[[l]]
      return(out)
    }))
  }))
  
  var_full <- solve(-HC - var_score)
  # if (max(abs(- solve(HC) - vcov(fit))) > 1e-5){
  #   browser()
  #   diag_HC <- sqrt(diag(-solve(HC)))
  #   diag_vcov <- sqrt(diag(vcov(fit)))
  #   warning('Difference between Louis and vcov.glm over 1e-5.')
  #   print(range(diag_HC - diag_vcov))
  # }
  return(list(posterior = posterior, model = fit, vcov = var_full))
}
