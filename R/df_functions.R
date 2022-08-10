###
# Calculate degrees of freedom for use in information criterion
###

# Use the Khalili and Chen
calculate_df_kc <- function(X, ridge, beta, binding_null_basis = NULL){
  
  # mm <<- mget(ls())
  
  K <- ncol(beta)
  if (!is.null(binding_null_basis)){
    
    ridge <- t(binding_null_basis) %*% ridge %*% binding_null_basis
    X <- X %*% binding_null_basis
    
    proj_beta <- matrix(as.vector(solve(crossprod(binding_null_basis), t(binding_null_basis) %*% beta)))
  }else{
    proj_beta <- beta
  }

  #Do tr[(X^T diag[p * 1-p] X + R)^{-1} X^T diag(...)X] in C++
  
  # meat <- plogis(as.vector(X %*% proj_beta))
  # meat <- Diagonal(x = meat * (1 - meat))
  # sum(diag(solve(t(X) %*% meat %*% X + ridge, t(X) %*% meat %*% X)))
  
  df_kc <- trace_df_cpp(X = drop0(X), beta = proj_beta, ridge = ridge)
  
  return(df_kc)
}

# Use the number of free parameters
calculate_df_freeparam <- function(blocked_X, binding_null_basis, outer = TRUE){

  rank_projX <- rank_via_null(blocked_X %*% binding_null_basis, outer = outer)
  
  df_model <- rank_projX
  return(df_model)
}

calculate_df_EM <- function(X, ridge, omega_weight, binding_null_basis = NULL){
  

  if (length(omega_weight) != nrow(X)){stop('omega misaligned')}
  K <- ncol(beta)
  if (!is.null(binding_null_basis)){
    
    ridge <- t(binding_null_basis) %*% ridge %*% binding_null_basis
    X <- X %*% binding_null_basis
  }
  
  df_EM <- sum(diag(solve(t(X) %*% Diagonal(x = omega_weight) %*% X + ridge, 
                          t(X) %*% Diagonal(x = omega_weight) %*% X)))
  return(df_EM)
}

calculate_df <- function(beta, X, binding_restrictions, group_mapping, 
     phi, W, y, lambda, gamma, 
     ridge_phi, rank_F, Fmatrix, K, 
     tol_epsilon, adaptive_weight, 
     ridge_beta, single_intercept, clip_tiny){
  
  args <- as.list(match.call())
  args <- args[names(args) != '']
  
  if (lambda == 0){
    if (single_intercept){
      binding_null_basis <- lapply(1:K, FUN=function(k){sparse_diag(rep(1, ncol(X) - 1))})
      X_list <- lapply(1:K, FUN=function(i){X[, -ncol(X), drop = F]})
      rank_F <- (ncol(X)-1)/2
    }else{
      binding_null_basis <- lapply(1:K, FUN=function(k){sparse_diag(rep(1, ncol(X)))})
      X_list <- lapply(1:K, FUN=function(i){X})
      rank_F <- ncol(X)/2
    }
  }else{
    if (single_intercept){
      binding_null_basis <- lapply(binding_restrictions, FUN=function(r_k){
        binding_k <- do.call('rbind', mapply(r_k, Fmatrix, SIMPLIFY = FALSE, FUN=function(i, F_j){
          do.call('rbind', F_j[i])
        }))
        binding_k <- binding_k[,-ncol(X)]
        if (is.null(binding_k)){
          return(sparse_diag(rep(1, ncol(X)-1)))
        }else{
          drop0(zapsmall(calculate_nullspace_basis(binding_k), clip_tiny))
        }
      })
      X_list <- lapply(binding_null_basis, FUN=function(i){X[, -ncol(X)] %*% i})
      
    }else{
      binding_null_basis <- lapply(binding_restrictions, FUN=function(r_k){
        binding_k <- do.call('rbind', mapply(r_k, Fmatrix, SIMPLIFY = FALSE, FUN=function(i, F_j){
          do.call('rbind', F_j[i])
        }))
        if (is.null(binding_k)){
          return(sparse_diag(rep(1, ncol(X))))
        }else{
          drop0(zapsmall(calculate_nullspace_basis(binding_k), clip_tiny))
        }
      })
      
      X_list <- lapply(binding_null_basis, FUN=function(i){X %*% i})
    }
  }

  if (single_intercept){
    beta_list <- mapply(binding_null_basis, 1:K, FUN=function(i,k){
      solve(crossprod(i), t(i) %*% beta[-nrow(beta),k,drop=F])
    })
  }else{
    beta_list <- mapply(binding_null_basis, 1:K, FUN=function(i,k){
      solve(crossprod(i), t(i) %*% beta[,k,drop=F])
    })
    
  }

  loglik.k <- plogis(as.matrix(X %*% beta) * (2 * y - 1), log.p = TRUE) 
  
  E.z <- calculate_posterior_zi(loglik.k = loglik.k, 
      group_mapping = group_mapping, K = K, W = W, ncol_W = ncol(W), 
      pi = NA, phi = phi)
  
  args$E.z <- E.z
  args$X_list <- X_list
  args$binding_null_basis <- binding_null_basis
  args$beta_list <- beta_list
  args$rank_F <- rank_F
  if (single_intercept){
    args$mu <- as.numeric(beta[nrow(beta),1])
  }else{
    args$mu <- NULL
  }

  args$ll_type <- 'no_prior'
  args <- args[!(names(args) %in% c('beta', 'X', 'binding_restrictions', 'clip_tiny', 'single_intercept'))]
  ll <- do.call('standard_error_louis', args = args)
  args$ll_type <- 'no_ll'
  prior <- do.call('standard_error_louis', args = args)

  DfHessianLogPrior <- drop0(zapsmall(-1 * with(prior, -HC - VScore), 15))
  DfHessianLogLik <- drop0(zapsmall(-1 * with(ll, -HC - VScore), 15))
  
  # DfHessianLogPrior <<- DfHessianLogPrior
  # DfHessianLogLik <<- DfHessianLogLik
  
  # #ADD AS UNIT TEST
  # all.equal(DfHessianLogLik + DfHessianLogPrior, DfHessianLogPosterior)
  gen_df <- sum(diag(solve(DfHessianLogLik + DfHessianLogPrior, DfHessianLogLik)))
  # gen_df <- fast_sparse_df(LL = DfHessianLogLik, PR = DfHessianLogPrior)
  return(gen_df)

}

calculate_df_kc_new <- function(X, ridge, beta, E.z, group_mapping, binding_null_basis = NULL){
  
  K <- ncol(beta)
  if (!is.null(binding_null_basis)){
    
    ridge <- t(binding_null_basis) %*% ridge %*% binding_null_basis
    X <- X %*% binding_null_basis
    
    proj_beta <- matrix(as.vector(solve(crossprod(binding_null_basis), t(binding_null_basis) %*% beta)))
  }else{
    proj_beta <- beta
  }
  
  p <- plogis(as.vector(X %*% proj_beta))
  meat <- t(X) %*% Diagonal(x = as.vector(group_mapping %*% E.z) * p * (1-p)) %*% X
  df <- diag(solve(meat + ridge, meat))
  df <- sum(df)
  return(df)
}
