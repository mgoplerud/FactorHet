
# Not exported:
# Function for creating moderator data
#' @importFrom stats model.frame terms .getXlevels model.matrix
create_moderator <- function(moderator, group, unique_group, design,
  args = NULL, check_rank = TRUE){
  
  if (is.null(args)){
    #A simplified version of what lm does
    mf_mod <- model.frame(moderator, design)
    tt_mod <- terms(mf_mod)
    tt_xlevels <- .getXlevels(tt_mod, mf_mod)
    
    W <- model.matrix(tt_mod, mf_mod)
    if (ncol(W) == 1 & all(colnames(W) == '(Intercept)')){
      W_full <- W
    }else{
      W_full <- do.call('cbind', lapply(1:ncol(mf_mod), FUN=function(i){
        model.matrix(~ 0 + ., data = mf_mod[i])
      }))
    }

  }else{
    #A simplified version of what predict.lm does
    tt_mod <- args$terms
    mf_mod <- model.frame(args$terms, design, xlev = args$xlev)
    tt_xlevels <- args$xlev
    W <- model.matrix(object = args$terms, mf_mod, contrasts.arg = args$contrasts)
    W_full <- NULL
  }
  contrasts_mod <- attr(W, 'contrasts')
  names_mod <- colnames(W)
  
  if (nrow(W) != nrow(design)){
    stop('No missing data allowed for moderator variables')
  }
  #Get unique for each group
  if (any(sapply(lapply(split(as.data.frame(W), group), unique), nrow) != 1)){
    stop('Moderators are not unique by ~ group.')
  }
  
  # W <- unique(group_by(data.frame(group, as.data.frame(W), stringsAsFactors = F, check.names = FALSE)))
  W <- unique(data.frame(group, as.data.frame(W), 
    stringsAsFactors = F, check.names = FALSE))
  
  if (!is.null(W_full)){
    # W_full <- unique(group_by(data.frame(group, as.data.frame(W_full), stringsAsFactors = F, check.names = FALSE)))
    W_full <- unique(data.frame(group, as.data.frame(W_full),
        stringsAsFactors = F, check.names = FALSE))
    W_full <- as.matrix(W_full[,-1,drop=F])
  }
  W_group <- as.vector(W[,1,drop=T])
  W <- as.matrix(W[,-1,drop=F])
  stopifnot(all.equal(W_group, unique_group))
  
  if (check_rank){
    if (ncol(W) != rank_via_null(W)){
      stop('moderator formula must make full rank design.')
    }
  }
  ncol_W <- ncol(W)

  return(list(W = W, ncol_W = ncol_W, W_full = W_full,
    args_W = list(terms = tt_mod, xlev = tt_xlevels, contrasts = contrasts_mod, names = names_mod))
  )
}

# Create ridge penalty following Tutz and XX
make_TMatrix <- function(K){
  TMatrix <- matrix(-1/K, K-1, K-1)
  diag(TMatrix) <- (K-1)/K
  return(TMatrix)  
}


update_phi <- function(phi, W, group_E.prob, K, 
                       ridge_phi, weights_W,
                       maxit_mod, extra_opt_args){
  
  ridge_penalty <- make_TMatrix(K) * ridge_phi
  
  norm_weights <- as.vector(weights_W/sum(weights_W))
  
  pi_bar <- colSums(Diagonal(x = norm_weights) %*% softmax_matrix(W %*% t(phi)))
  
  vec_phi <- as.vector(phi[-1,,drop=F])
  
  optim_method <- extra_opt_args$method
  
  if (is.null(optim_method)){
    stop('"method" must be provided in FactorHet_control(optim_phi_controls = list(method = "..."), ...)')
  }else{
    extra_opt_args <- extra_opt_args[
      !(names(extra_opt_args) %in% c('method', 'force_nocpp'))
    ]
  }
  
  if (optim_method == 'lib_lbfgs'){
    
    if (!is.null(maxit_mod)){
      it_arg <- list(max_iterations = maxit_mod)
    }else{
      it_arg <- NULL
    }
    
    args <- c(list(
      call_eval = function(...){-1 * objective_phi(...)},
      call_grad = function(...){-1 * gradient_phi(...)},
      vars = vec_phi, invisible = 1,
      W = W, weights_W = weights_W,
      group_E.prob = group_E.prob, ridge_penalty =ridge_penalty, K = K),
      extra_opt_args,
      it_arg
    )
    
    optim_moderator <- do.call('lbfgs', args)
    
  }else{
    
    if (!is.null(maxit_mod)){
      required_args <- list(fnscale = -1, maxit = maxit_mod)
    }else{
      required_args <- list(fnscale = -1)
    }
    
    optim_moderator <- tryCatch(optim(
      par = vec_phi, fn = objective_phi, gr = gradient_phi,
      W = W, weights_W = weights_W,
      group_E.prob = group_E.prob, ridge_penalty =ridge_penalty, K = K,
      method = optim_method,
      control = c(required_args, extra_opt_args)
    ), error = function(e){NULL})
    if (is.null(optim_moderator)){
      warning('BFGS used instead of optimization method provided to optim', 
              immediate. = TRUE)
      optim_moderator <- optim(
        par = vec_phi, fn = objective_phi, gr = gradient_phi,
        W = W, weights_W = weights_W,
        group_E.prob = group_E.prob, ridge_penalty =ridge_penalty, K = K,
        method = optim_method,
        control = c(list(fnscale = -1, maxit = maxit_mod), extra_opt_args)
      )
      
    }
  }
  
  prior_obj <- objective_phi(vec_phi = vec_phi, 
      W = W, K = K, weights_W = weights_W,
      group_E.prob = group_E.prob, ridge_penalty =ridge_penalty)
  
  if (optim_moderator$value < prior_obj){
    warning('Optim Failed for Phi')
  }
  phi <- rbind(0, matrix(optim_moderator$par, nrow = K- 1))
  
  pi_bar <- colSums(Diagonal(x = norm_weights) %*% softmax_matrix(W %*% t(phi)))
  
  return(list(phi = phi, pi_bar = pi_bar))
}

objective_phi <- function(vec_phi, group_E.prob, W, K, weights_W, ridge_penalty){
  phi <- matrix(vec_phi, nrow = K - 1)
  
  objective <- sum( Diagonal(x = weights_W) %*% group_E.prob * log(softmax_matrix(cbind(0, W %*% t(phi)))) )
  regularization <- -1/2 * sum(diag(t(phi) %*% ridge_penalty %*% phi))
  
  return(objective + regularization)
}

gradient_phi <-  function(vec_phi, group_E.prob, W, K, weights_W, ridge_penalty){
  
  phi <- matrix(vec_phi, nrow = K - 1)
  #Get pi_ik
  pi_ik <- softmax_matrix(cbind(0, W %*% t(phi)))
  
  weight_grad_ll <- Diagonal(x = weights_W) %*% (group_E.prob - pi_ik)[,-1, drop = F]
  
  grad_ll_phi <- apply(weight_grad_ll, MARGIN = 2, FUN=function(w_k){
    colSums(Diagonal(x = w_k) %*% W)
  })
  
  if (ncol(W) == 1){
    grad_ll_phi <- matrix(as.vector(grad_ll_phi))
  }else{
    grad_ll_phi <- t(grad_ll_phi)
  }
  grad_reg_phi <- apply(phi, MARGIN = 2, FUN=function(i){-1 * ridge_penalty %*% i})
  
  all_grad_vecphi <- as.vector(grad_ll_phi + grad_reg_phi)
  return(all_grad_vecphi)
}

reverse_softmax <- function(pi, baseline = NULL){
  if (is.null(baseline)){
    baseline <- 1
  }
  log_pi <- log(pi)
  log_pi <- log_pi - log_pi[baseline]
  return(log_pi)
}

#' @importFrom lbfgs lbfgs
#' @importFrom stats optim
update_moderator <- function(phi, beta, W, Fmatrix, group_E.prob, K, 
  weights_W, 
  ridge_phi, ridge_beta, lambda, gamma, rank_F, adaptive_weight, single_intercept,
  maxit_mod = 10, 
  use_grad = TRUE, extra_opt_args = NULL){
  
  #If using ridge, must use manual optim.
  if (ridge_phi > 0){
    ridge_penalty <- make_TMatrix(K) * ridge_phi
  }else{
    ridge_penalty <- diag(rep(0, K-1))
  }
  
  
  if (lambda == 0){
    b_r <- sapply(list_from_cols(beta), FUN=function(b){
      as.numeric(t(b) %*% ridge_beta %*% b)
    })
    b_r <- 1/2 * b_r
    power_pi <- gamma
    
  }else{
    b_r <- mapply(list_from_cols(beta), adaptive_weight, FUN=function(b, aw){
      F_prior_kernel(beta = b, Fmat = Fmatrix, aw = aw)
    })
    b_r <- -b_r
    power_pi <- gamma * (lambda > 0)
  }
  n_W <- nrow(W)
  
  vec_phi <- as.vector(phi[-1,,drop=F])
  
  if (use_grad){
    optim_gr <- cpp_gradient_phi
  }else{
    optim_gr <- NULL
  }
  
  optim_fn <- objective_moderator
  optim_method <- extra_opt_args$method
  if (is.null(optim_method)){
    stop('"method" must be provided in FactorHet_control(optim_phi_controls = list(method = "..."), ...)')
  }else{
    if (!is.null(extra_opt_args$force_nocpp)){
      if (extra_opt_args$force_nocpp == TRUE){
        optim_gr <- gradient_moderator
        optim_fn <- objective_moderator
      }
    }
    extra_opt_args <- extra_opt_args[
      !(names(extra_opt_args) %in% c('method', 'force_nocpp'))
    ]
  }
  
  norm_weights <- as.vector(weights_W/sum(weights_W))
  weights_W <- as.vector(weights_W)
  
  if (optim_method == 'lib_lbfgs'){
    if (!is.null(maxit_mod)){
      it_arg <- list(max_iterations = maxit_mod)
    }else{
      it_arg <- NULL
    }
    args <- c(list(
      call_eval = function(...){-1 * optim_fn(...)},
      call_grad = function(...){-1 * optim_gr(...)},
      vars = vec_phi, invisible = 1,
      norm_weights = norm_weights, weights_W = weights_W,
      K = K, W = W, group_E_prob = group_E.prob, ridge_penalty = ridge_penalty,
      gamma = gamma, rank_F = rank_F, power_pi = power_pi, b_r = b_r, 
      lambda = lambda),
      extra_opt_args,
      it_arg
    )
    
    optim_moderator <- do.call('lbfgs', args)
    
  }else{
    
    if (!is.null(maxit_mod)){
      required_args <- list(fnscale = -1, maxit = maxit_mod)
    }else{
      required_args <- list(fnscale = -1)
    }
    
    optim_moderator <- tryCatch(optim(par = vec_phi,
        fn = optim_fn, gr = optim_gr,
        method = optim_method,
        norm_weights = norm_weights, weights_W = weights_W,
        K = K, W = W, group_E_prob = group_E.prob, ridge_penalty = ridge_penalty,
        gamma = gamma, rank_F = rank_F, power_pi = power_pi, b_r = b_r, lambda = lambda,
        control = c(required_args, extra_opt_args)
    ), error = function(e){NULL})
    if (is.null(optim_moderator)){
      warning('BFGS used instead of optimization method provided to optim', 
        immediate. = TRUE)
      optim_moderator <- optim(par = vec_phi,
         fn = optim_fn, gr = optim_gr,
         method = 'BFGS',
         norm_weights = norm_weights, weights_W = weights_W,
         K = K, W = W, group_E_prob = group_E.prob, ridge_penalty = ridge_penalty,
         gamma = gamma, rank_F = rank_F, power_pi = power_pi, b_r = b_r, lambda = lambda,
         control = c(list(fnscale = -1, maxit = maxit_mod), extra_opt_args)
      )
      
    }
  }
  
  prior_moderator <- objective_moderator(
    par = vec_phi, weights_W = weights_W, norm_weights = norm_weights,
    K = K, n_W = n_W, W = W, group_E_prob = group_E.prob, ridge_penalty = ridge_penalty, 
    gamma = gamma, rank_F = rank_F, power_pi = power_pi, b_r = b_r, lambda = lambda
  )
  if (optim_moderator$value - prior_moderator < -1e-7){
    warning('Optim Failed for Phi')
  }
  phi <- rbind(0, matrix(optim_moderator$par, nrow = K-1))
  
  #Average of Posterior Predictive Probs
  pi_bar <- colSums((Diagonal(x = norm_weights) %*% softmax_matrix(W %*% t(phi))))
  
  return(list(phi = phi, pi_bar = pi_bar, b_r = b_r, convergence = optim_moderator$convergence))
}


objective_moderator <- function(par, K, n_W, W, 
      weights_W, norm_weights,
      group_E_prob, ridge_penalty, gamma, rank_F, power_pi, b_r, lambda){

  phi <- matrix(par, nrow = K - 1)
  #Get multinomial loglik with weighted outcomes.
  pi_ik <- softmax_matrix(cbind(0, W %*% t(phi)))
  loglik_zik <- sum(Diagonal(x = weights_W) %*% (group_E_prob * log(pi_ik)))
  #Set pi_bar to the average of all posterior predictive probabilities
  pi_bar <- colSums(Diagonal(x = norm_weights) %*% pi_ik)
  
  #Get prior / regularization on phi
  regularization_phi <- -1/2 * sum(diag(t(phi) %*% ridge_penalty %*% phi))
  
  #Prior on beta, ignoring K m \ln(\lambda) if lambda >0 and constant if lambda = 0
  if (lambda == 0){
    regularization_beta <- sum(log(pi_bar) * power_pi * rank_F - pi_bar^power_pi * b_r)
  }else{
    regularization_beta <- sum(log(pi_bar) * power_pi * rank_F - lambda * pi_bar^power_pi * b_r)
  }

  objective <- loglik_zik + regularization_phi + regularization_beta
  
  return(objective)
}

#' @importFrom Matrix Diagonal
gradient_moderator <-  function(par, K, W,  weights_W, norm_weights,
    group_E_prob, ridge_penalty, gamma, rank_F, power_pi, b_r, lambda){
  
  n_W <- nrow(W)
  
  phi <- matrix(par, nrow = K - 1)

  #Get pi_ik
  pi_ik <- softmax_matrix(cbind(0, W %*% t(phi)))
  pi_bar <- colSums(Diagonal(x = norm_weights) %*% pi_ik)
  
  weight_grad_ll <- (group_E_prob - pi_ik)[,-1, drop = F]
  
  grad_ll_phi <- apply(weight_grad_ll, MARGIN = 2, FUN=function(w_k){
    colSums(Diagonal(x = weights_W * w_k) %*% W)
  })
  
  if (ncol(W) == 1){
    grad_ll_phi <- matrix(as.vector(grad_ll_phi))
  }else{
    grad_ll_phi <- t(grad_ll_phi)
  }

  grad_reg_phi <- apply(phi, MARGIN = 2, FUN=function(i){-1 * ridge_penalty %*% i})
  
  if (gamma != 0){
    if (lambda == 0){
      weight_k <- (power_pi * rank_F)/pi_bar - gamma * pi_bar^(gamma - 1) * b_r
    }else{
      weight_k <- (power_pi * rank_F)/pi_bar - lambda * gamma * pi_bar^(gamma - 1) * b_r
    }
    
    branch_ik <- list_from_cols(pi_ik)
    w_pi_ik <- sapply(branch_ik, FUN=function(i){colMeans(Diagonal(x = i) %*% W)})
    if (ncol(W) == 1){
      w_pi_ik <- as.matrix(t(w_pi_ik))
    }
    array_weights <- array(NA, dim = c(K, K, ncol(W)))
    for (i in 2:K){
      for (j in 1:K){
        cross_ij <- colMeans(Diagonal(x = -branch_ik[[i]] * branch_ik[[j]]) %*% W)
        if (i == j){
          cross_ij <- cross_ij + w_pi_ik[,i]
        }
        array_weights[i,j,] <- weight_k[j]  * cross_ij
      }
    }
    
    array_weights <- apply(array_weights, MARGIN = 1, colSums)
    if (ncol(W) == 1){
      array_weights <- as.matrix(t(array_weights))
    }
    all_grad_reg <- grad_reg_phi + t(array_weights[,-1, drop = F])
    # glb <- gamma * lambda * b_r * pi_bar^gamma
    # grad_reg_xi <- gamma * rank_F - glb - pi_bar * (K * gamma * rank_F - sum(glb))    
    
  }else{
    all_grad_reg <- grad_reg_phi
  }
  all_grad <- as.vector(grad_ll_phi + all_grad_reg)
  return(all_grad)
}
