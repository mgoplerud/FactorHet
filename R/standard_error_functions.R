
#Gives the gradient of \bar{pi}_{k'} w.r.t. each \phi_k
#in list form where the list is for \bar{pi}_{k'} and the
#matrix is for each \bm{\phi}_k
gradient_pibar <- function(phi, W, K, norm_weights){
  pi_ik <- softmax_matrix(W %*% t(phi))
  all_out <- lapply(1:K, FUN=function(i){matrix(NA, nrow = ncol(W), ncol = K)})
  for (k in setdiff(1:K, 1)){
    for (kprime in 1:K){
      all_out[[kprime]][,k] <- (k == kprime) * colSums(Diagonal(x = norm_weights * pi_ik[,k]) %*% W) - colSums(Diagonal(x = norm_weights * pi_ik[,k] * pi_ik[,kprime]) %*% W)
    }
  }
  return(all_out)
}
#Gives the Hessian of \bar{pi}_{k'} w.r.t. each \phi_l, \phi_m
#in list where the list is for each \bar{pi}_{k'} and the matrix
#of lists
hessian_pibar <- function(phi, W, K, norm_weights){
  pi_ik <- softmax_matrix(W %*% t(phi))
  all_out <- lapply(1:K, FUN=function(i){l <- lapply(1:K^2, FUN=function(i){matrix(NA, nrow = 0, ncol = 0)}); dim(l) <- c(K, K); return(l)})
  for (kprime in 1:K){
    for (l in setdiff(1:K, 1)){
      for (m in setdiff(1:K, 1)){
        weight <- ((kprime == l) - pi_ik[,l]) * ((kprime == m) - pi_ik[,m]) + 
          -((l == m) - pi_ik[,m]) * pi_ik[,l]
        weight <- pi_ik[,kprime] * weight
        all_out[[kprime]][l,m][[1]] <- as.matrix(t(W) %*% Diagonal(x = norm_weights * weight) %*% W)
      }
    }
  }
  return(all_out)
}
#Gradient as above but for log of pi-bar
gradient_ln_pibar <- function(phi, W, K, norm_weights){
  pi_ik <- softmax_matrix(W %*% t(phi))
  pi_bar <- colSums(Diagonal(x = norm_weights) %*% pi_ik)
  all_out <- lapply(1:K, FUN=function(i){matrix(NA, nrow = ncol(W), ncol = K)})
  for (k in setdiff(1:K, 1)){
    for (kprime in 1:K){
      all_out[[kprime]][,k] <- 1/pi_bar[kprime] * (
        (k == kprime) * colSums(Diagonal(x = norm_weights * pi_ik[,k]) %*% W) - 
          colSums(Diagonal(x = norm_weights * pi_ik[,k] * pi_ik[,kprime]) %*% W)
      )
    }
  }
  return(all_out)
}
#Hessian as above but for log of pi-bar
hessian_ln_pibar <- function(phi, W, K, norm_weights){
  pi_ik <- softmax_matrix(W %*% t(phi))
  pi_bar <- colSums(Diagonal(x = norm_weights) %*% pi_ik)
  
  grad_bar <- gradient_pibar(phi, W, K, norm_weights)
  hess_bar <- hessian_pibar(phi, W, K, norm_weights)
  
  all_out <- lapply(1:K, FUN=function(i){l <- lapply(1:K^2, FUN=function(i){matrix(NA, nrow = 0, ncol = 0)}); dim(l) <- c(K, K); return(l)})
  
  for (kprime in 1:K){
    for (l in setdiff(1:K, 1)){
      for (m in setdiff(1:K, 1)){
        all_out[[kprime]][l,m][[1]] <- -(pi_bar[kprime])^(-2) * matrix(grad_bar[[kprime]][,l]) %*% t(matrix(grad_bar[[kprime]][,m])) +
          hess_bar[[kprime]][l,m][[1]] * 1/pi_bar[kprime]
      }
    }
  }
  return(all_out)
}

hessian_complete <- function(y, X_list, W, lambda, gamma, mu, beta_list, phi, ridge_phi,
   K, Fmatrix, rank_F, adaptive_weight, group_mapping, binding_null_basis,
   weights = weights, weights_W = weights_W,
   E.z, tol_epsilon, ll_type, ridge_beta, fisher = FALSE, oelker = FALSE){
  
  k_range <- setdiff(1:K, 1)
  
  E.z <- group_mapping %*% E.z
  
  if (fisher){
    E.z <- group_mapping %*% softmax_matrix(W %*% t(phi))
  }
  
  single_intercept <- !is.null(mu)
  if (single_intercept){
    split_p <- mapply(X_list, beta_list, SIMPLIFY = FALSE, FUN=function(x_k, b_k){plogis(mu + as.vector(x_k %*% b_k))})
  }else{
    split_p <- mapply(X_list, beta_list, SIMPLIFY = FALSE, FUN=function(x_k, b_k){plogis(as.vector(x_k %*% b_k))})
  }

  norm_weights <- weights_W/sum(weights_W)
  pi_ik <- softmax_matrix(W %*% t(phi))
  pi_bar <- colSums((Diagonal(x = norm_weights) %*% pi_ik))
  
  ######
  # Get the terms for \beta_k, \beta_l
  ######
  
  h_beta <- mapply(1:K, beta_list, X_list, binding_null_basis,
   list_from_cols(E.z), SIMPLIFY = FALSE, FUN=function(k, b.k, X.k, basis.k, z.ik){
     p.ik <- split_p[[k]]
     
     if (ll_type == 'no_ll'){
       ll.meat <- Diagonal(x = rep(0, ncol(X.k)))
     }else{
       ll.meat <- -1 * t(X.k) %*% Diagonal(x = weights * p.ik * (1-p.ik) * z.ik) %*% X.k
     }
     
     prior.beta.meat <- 0
     if (ll_type != 'no_prior'){
       
       if (lambda > 0){
         
         if (single_intercept){
           Fmatrix_k <- lapply(Fmatrix, FUN=function(F_l){lapply(F_l, FUN=function(F_lj){
             t(basis.k) %*% F_lj[-nrow(F_lj),-nrow(F_lj)] %*% basis.k
           })})
         }else{
           Fmatrix_k <- lapply(Fmatrix, FUN=function(F_l){lapply(F_l, FUN=function(F_lj){
             t(basis.k) %*% F_lj %*% basis.k
           })})
         }

         prior.beta.meat <- mapply(Fmatrix_k, adaptive_weight[[k]], SIMPLIFY = FALSE, FUN=function(F_l,exo_j){
           out_l <- mapply(F_l, exo_j, SIMPLIFY = FALSE, FUN=function(F_g, exo_jg){
             
             row_g <- F_g %*% b.k
             meat <- as.numeric(t(b.k) %*% F_g %*% b.k) + tol_epsilon
             if (oelker){
               out_g <- F_g * meat^(-1/2)
             }else{
               out_g <- -(meat)^(-3/2) * tcrossprod(row_g) +
                 meat^(-1/2) * F_g
             }
             return(out_g * exo_jg)
           })
           return(Reduce('+', out_l))
         })
         prior.beta.meat <- Reduce('+', prior.beta.meat)
         prior.beta.meat <-  -lambda * pi_bar[k]^gamma * prior.beta.meat
       }else{
         ridge.k <- ridge_beta[[k]]
         if (!single_intercept){
           prior.beta.meat <- - pi_bar[k]^gamma * ridge.k
         }else{
           prior.beta.meat <- - pi_bar[k]^gamma * ridge.k[-ncol(ridge.k), -ncol(ridge.k)]
         }
       }
     }
     
     out.beta <- ll.meat + prior.beta.meat
     
     return(out.beta)
   })
  
  # Legacy
  # if (K == 1){
  #   
  #   h_mu <- -sum(split_p[[1]] * (1-split_p[[1]]))
  #   h_complete <- bdiag(h_mu, bdiag(h_beta))
  #   
  #   hcross_betamu <- -colSums(Diagonal(x = E.z[,1] * split_p[[1]] * (1-split_p[[1]])) %*% X_list[[1]])
  #   
  #   h_complete[1,-1] <- hcross_betamu
  #   h_complete[-1,1] <- hcross_betamu
  #   return(h_complete)
  # }
  
  #############
  # Get the terms for phi_k, phi_l
  ############

  hess_lnpi_bar <- hessian_ln_pibar(phi = phi, W = W, K = K, norm_weights = norm_weights)
  hess_pi_bar <- hessian_pibar(phi = phi, W = W, K = K, norm_weights = norm_weights)
  grad_pi_bar <- gradient_pibar(phi = phi, W = W, K = K, norm_weights = norm_weights)
  
  h_phi <- lapply(1:K^2, FUN=function(i){matrix(0, nrow=ncol(W),ncol=ncol(W))})
  dim(h_phi) <- c(K,K)
  
  for (k in k_range){
    for (l in k_range){
      if (ll_type != 'no_ll'){
        
        h_phi[k,l][[1]] <- h_phi[k,l][[1]] + -1 * as.matrix(t(W) %*% Diagonal(x = weights_W * ((k == l) - pi_ik[,k]) * pi_ik[,l]) %*% W)
      }
      
    }
  }
  
  if (ll_type != 'no_prior'){#Add in the structured sparse prior
    
    for (kprime in 1:K){
      
      b.kprime <- beta_list[[kprime]]
      
      if (lambda > 0){
        basis.k <- binding_null_basis[[kprime]]
        Fmatrix_k <- lapply(Fmatrix, FUN=function(F_l){lapply(F_l, FUN=function(F_lj){
          if (single_intercept){
            out.k <- t(basis.k) %*% F_lj[-nrow(F_lj),-nrow(F_lj)] %*% basis.k
          }else{
            out.k <- t(basis.k) %*% F_lj %*% basis.k
          }
          return(out.k)
        })})
        
        F_meat <- mapply(Fmatrix_k, adaptive_weight[[kprime]], SIMPLIFY = FALSE, FUN=function(F_l,exo_j){
          out_l <- mapply(F_l, exo_j, SIMPLIFY = FALSE, FUN=function(F_g, exo_jg){
            meat <- as.numeric(t(b.kprime) %*% F_g %*% b.kprime) + tol_epsilon
            return(sqrt(meat) * exo_jg)
          })
          return(Reduce('+', out_l))
        })
        F_meat <- -lambda * gamma * Reduce('+', F_meat)
      }else{
        ridge.kprime <- ridge_beta[[kprime]]
        if (single_intercept){
          ridge.kprime <- ridge.kprime[-ncol(ridge.kprime), -ncol(ridge.kprime)]
        }
        
        F_meat <- -1/2 * as.numeric(t(b.kprime) %*% ridge.kprime %*% b.kprime)
      }
      
      for (k in k_range){
        for (l in k_range){
          term_3 <- rank_F * gamma * hess_lnpi_bar[[kprime]][k,l][[1]]
          term_4 <- F_meat * pi_bar[kprime]^(gamma - 1) * hess_pi_bar[[kprime]][k,l][[1]]
          if (!(gamma %in% c(0,1))){
            term_4 <- term_4 + F_meat * (gamma - 1) * pi_bar[kprime]^(gamma - 2) * matrix(grad_pi_bar[[kprime]][,k]) %*% t(matrix(grad_pi_bar[[kprime]][,l]))
          }
          h_phi[k,l][[1]] <- h_phi[k,l][[1]] + as.matrix(term_3 + term_4)
        }
      }
    }
  }
  #Ridge prior on phi_k,phi_l
  if (K != 1){
    prec_prior_phi <- bdiag(lapply(seq_len(ncol(W)), FUN=function(i){ridge_phi * make_TMatrix(K)}))
    
    order_w <- as.vector(matrix(seq_len(ncol(W) * (K-1)), byrow = T, nrow = K - 1))
    Perm <- sparseMatrix(i=seq_len(ncol(W) * (K-1)), j =order_w, x = 1)
    
    prec_prior_phi <- t(Perm) %*% prec_prior_phi %*% Perm
    if (ll_type == 'no_prior'){
      prec_prior_phi[,] <- 0
      prec_prior_phi <- drop0(prec_prior_phi)
    }
  }else{
    prec_prior_phi <- matrix(nrow = 0, ncol = 0)
  }
  
  ############
  # Get the terms for \beta_k, phi_l (only when lambda > 0)
  ############
  h_cross <- lapply(1:K^2, FUN=function(i){matrix(NA, nrow=0, ncol = 0)})
  dim(h_cross) <- c(K, K)
  
  if (ll_type != 'no_prior'){
    for (k in 1:K){
      
      if (lambda > 0){
        
        b.k <- beta_list[[k]]
        basis.k <- binding_null_basis[[k]]
        Fmatrix_k <- lapply(Fmatrix, FUN=function(F_l){lapply(F_l, FUN=function(F_lj){
          if (single_intercept){
            out.k <- t(basis.k) %*% F_lj[-nrow(F_lj),-nrow(F_lj)] %*% basis.k
          }else{
            out.k <- t(basis.k) %*% F_lj %*% basis.k
          }
          return(out.k)
        })})
        
        weight_k <- -lambda * gamma * pi_bar[k]^(gamma-1)
        beta.weight <- mapply(Fmatrix_k, adaptive_weight[[k]], SIMPLIFY = FALSE, FUN=function(F_l,exo_j){
          out_l <- mapply(F_l, exo_j, SIMPLIFY = FALSE, FUN=function(F_g, exo_jg){
            row_g <- F_g %*% b.k
            meat <- as.numeric(t(b.k) %*% F_g %*% b.k) + tol_epsilon
            return(meat^(-1/2) * exo_jg * row_g)
          })
        })
        beta.weight <- weight_k * Reduce('+', unlist(beta.weight))
        
      }else{
        
        weight_k <- -gamma * pi_bar[k]^(gamma-1)

        b.k <- beta_list[[k]]
        ridge.k <- ridge_beta[[k]]
        if (single_intercept){
          ridge.k <- ridge.k[-ncol(ridge.k), -ncol(ridge.k)]
        }
        
        beta.weight <- weight_k * (ridge.k %*% b.k)

      }
      
      for (l in k_range){#Loop for phi only over relevant
        h_cross[k,l][[1]] <- as.matrix(tcrossprod(as.vector(beta.weight), grad_pi_bar[[k]][,l]))
      }
    }
    
  }
  
  ##Build together all terms
  h_complete <- do.call('cbind', lapply(k_range, FUN=function(k){
    do.call('rbind', h_phi[-1,k])
  }))
  h_complete <- h_complete - prec_prior_phi
  
  X_size <- sapply(X_list, ncol)
  X_cumsize <- c(0, cumsum(X_size))
  
  if (single_intercept){
    
    h_complete <- bdiag(bdiag(h_beta), h_complete)
    if (ll_type == 'no_ll'){
      h_mu <- 0
    }else{
      h_mu <- sum(sapply(1:K, FUN=function(i){-sum(weights * E.z[,i] * split_p[[i]] * (1 - split_p[[i]]))}))
    }
    h_complete <- bdiag(h_mu, h_complete)

    for (k in 1:K){
      X_k <- X_list[[k]]
      if (ll_type == 'no_ll'){
        hcross_betamu <- rep(0, ncol(X_k))
      }else{
        hcross_betamu <- -colSums(Diagonal(x = weights * E.z[,k] * split_p[[k]] * (1-split_p[[k]])) %*% X_k)
      }
      if (length(hcross_betamu) != 0){
        h_complete[1, 1 + X_cumsize[k]  + 1:(X_size[k])] <- hcross_betamu
        h_complete[1 + X_cumsize[k]  + 1:(X_size[k]), 1] <- hcross_betamu
      }
    }
    if (ll_type != 'no_prior'){
      for (k in k_range){
        hcross_k <- do.call('rbind', lapply(h_cross[,k], FUN=function(i){i}))
        h_complete[1 + seq_len(X_cumsize[K + 1]), 1 + X_cumsize[K + 1] + ncol(W) * (k-2) + seq_len(ncol(W))] <- hcross_k
        h_complete[1 + X_cumsize[K + 1] + ncol(W) * (k-2) + seq_len(ncol(W)), 1 + seq_len(X_cumsize[K + 1])] <- t(hcross_k)
      }
    }
  }else{
    h_complete <- bdiag(bdiag(h_beta), h_complete)
    if (ll_type != 'no_prior'){
      for (k in k_range){
        hcross_k <- do.call('rbind', h_cross[,k])
        h_complete[seq_len(X_cumsize[K + 1]), X_cumsize[K + 1] + ncol(W) * (k-2) + seq_len(ncol(W))] <- hcross_k
        h_complete[X_cumsize[K + 1] + ncol(W) * (k-2) + seq_len(ncol(W)), seq_len(X_cumsize[K + 1])] <- t(hcross_k)
      }
      
    }
  }
  
  return(h_complete)
}

variance_score <- function(y, X_list, W, lambda, gamma, mu, beta_list, phi, 
       binding_null_basis, ridge_beta, weights, weights_W,
       group_mapping, K, Fmatrix, adaptive_weight, rank_F, ridge_phi,
       E.z, tol_epsilon, ll_type, fisher = FALSE, oelker = FALSE){
  
  if (ll_type == 'no_ll'){
    size_out <- sum(sapply(X_list, ncol)) + 1 + (K-1) * ncol(phi)
    var_score <- sparseMatrix(i = 1, j = 1, x = 0, dims = rep(size_out, 2))
    var_score <- drop0(var_score)
    return(var_score)
  }
  
  single_intercept <- !is.null(mu)
  k_range <- setdiff(1:K, 1)
  
  pi_ik <- softmax_matrix(W %*% t(phi))
  
  if (single_intercept){
    split_p <- mapply(X_list, beta_list, SIMPLIFY = FALSE, FUN=function(x_k, b_k){plogis(mu + as.vector(x_k %*% b_k))})
  }else{
    split_p <- mapply(X_list, beta_list, SIMPLIFY = FALSE, FUN=function(x_k, b_k){plogis(as.vector(x_k %*% b_k))})
  }
  var_beta <- var_cross <- var_phi <- lapply(1:K^2, FUN=function(i){matrix(nrow=0,ncol=0)})
  dim(var_phi) <- dim(var_cross) <- dim(var_beta) <- c(K, K)
  
  for (k in 1:K){
    for (l in 1:K){
      X_k <- X_list[[k]]
      X_l <- X_list[[l]]
      
      hX_k <- t(group_mapping) %*% Diagonal(x = (y - split_p[[k]])) %*% X_k
      hX_l <- t(group_mapping) %*% Diagonal(x = (y - split_p[[l]])) %*% X_l
      
      mod_weight <- E.z[,k] * ( (k == l) - E.z[,l])
      
      var_beta[k,l][[1]] <- t(hX_k) %*% Diagonal(x = weights_W * mod_weight) %*% hX_l
      
      if (l != 1){
        var_cross[k,l][[1]] <- as.matrix(t(hX_k) %*% Diagonal(x = weights_W * mod_weight) %*% W)
        if  (k != 1){
          var_phi[k,l][[1]] <- as.matrix(t(W) %*% Diagonal(x = weights_W * mod_weight) %*% W)
        }
      }
    }
  }
  
  var_beta <- do.call('cbind', lapply(1:K, FUN=function(i){do.call('rbind', var_beta[,i])}))
  
  if (K == 1){
    var_phi <- matrix(nrow=0,ncol=0)
  }else{
    var_phi <- do.call('cbind', lapply(k_range, FUN=function(i){do.call('rbind', var_phi[-1,i])}))
  }
  
  
  
  var_score <- bdiag(var_beta, var_phi)      
  
  X_size <- sapply(X_list, ncol)
  X_cumsize <- c(0, cumsum(X_size))
  
  if (single_intercept){
    
    var_mu <- sum(sapply(1:K, FUN=function(k){
      sapply(1:K, FUN=function(l){
        
        hX_k <- rowSums(t(group_mapping) %*% Diagonal(x = (y - split_p[[k]])))
        hX_l <- rowSums(t(group_mapping) %*% Diagonal(x = (y - split_p[[l]])))
        
        sum( weights_W * E.z[,k] * ((k==l) - E.z[,l]) * hX_l * hX_k)
      })
    }))  
    
    var_score <- bdiag(var_mu, var_score)
    var_score <- as.matrix(var_score)
    for (k in k_range){
      cross_phimu_k <- lapply(1:K, FUN=function(kprime){
        hX_kprime <- rowSums(t(group_mapping) %*% Diagonal(x = (y - split_p[[kprime]])))
        t(hX_kprime) %*% Diagonal(x = weights_W * E.z[,k] * ((k==kprime)- E.z[,kprime])) %*% W
      })
      cross_phimu_k <- Reduce('+', cross_phimu_k)
      var_score[1,1 + X_cumsize[K+1] + (k-2) * ncol(W) +  seq_len(ncol(W))] <- as.matrix(cross_phimu_k)
      var_score[1 + X_cumsize[K+1] + (k-2) * ncol(W) +  seq_len(ncol(W)), 1] <- as.matrix(cross_phimu_k)
    }
    
    for (k in 1:K){
      X_k <- X_list[[k]]
      if (ncol(X_k) == 0){next}
      
      cross_betamu_k <- lapply(1:K, FUN=function(kprime){
        hX_k <- t(group_mapping) %*% Diagonal(x = (y - split_p[[k]])) %*% X_k
        hX_kprime <- rowSums(t(group_mapping) %*% Diagonal(x = (y - split_p[[kprime]])))
        t(hX_k) %*% Diagonal(x = weights_W * E.z[,k] * ((k==kprime) - E.z[,kprime])) %*% hX_kprime
      })
      cross_betamu_k <- Reduce('+', cross_betamu_k)
      var_score[1, 1 + X_cumsize[k] + seq_len(X_size[k])] <- as.matrix(cross_betamu_k)
      var_score[1 + X_cumsize[k] + seq_len(X_size[k]), 1] <- as.matrix(cross_betamu_k)
    }

    for (k in k_range){
      var_cross_k <- as.matrix(do.call('rbind', var_cross[,k]))
      var_score[1 + seq_len(X_cumsize[K+1]), 1 + X_cumsize[K+1] + seq_len(ncol(W)) + ncol(W) * (k-2)] <- var_cross_k
      var_score[1+ X_cumsize[K+1] + seq_len(ncol(W)) + ncol(W) * (k-2), 1 + seq_len(X_cumsize[K+1])] <- t(var_cross_k)
    }
  }else{
    for (k in k_range){
      var_cross_k <- do.call('rbind', var_cross[,k])
      var_score[seq_len(X_cumsize[K+1]), X_cumsize[K+1] + seq_len(ncol(W)) + ncol(W) * (k-2)] <- var_cross_k
      var_score[X_cumsize[K+1] + seq_len(ncol(W)) + ncol(W) * (k-2), seq_len(X_cumsize[K+1])] <- t(var_cross_k)
    }
  }
  var_score <- drop0(var_score)
  return(var_score)
}

standard_error_louis <- function(y, X_list, Z, W, lambda, gamma, mu, beta_list, phi, ridge_phi,
                                 binding_null_basis, ridge_beta,
                                 K, Fmatrix, adaptive_weight, rank_F, group_mapping,
                                 weights = weights, weights_W = weights_W,
                                 E.z, tol_epsilon, ll_type = 'all', fisher = FALSE, oelker = FALSE){
  if (missing(mu)){
    mu <- NULL
  }
  args <- as.list(match.call())
  args <- args[names(args) != '']

  HC <- do.call('hessian_complete', args = args)
  VScore <- do.call('variance_score', args = args)
  
  I_obs <- -HC - VScore
  I_complete <- -HC
  
  if (ll_type == 'all'){
    
    vcov_obs <- solve(I_obs)
    vcov_complete <- solve(I_complete)
  }else{
    vcov_obs <- vcov_complete <- NULL
  }
  
  return(list('vcov' = vcov_obs, 'var_score' = VScore, 'I_C' = I_complete,
              'HC' = HC, 'VScore' = VScore))
}

standard_error_sandwich <- function(y, X_list, Z, W, lambda, gamma, mu, beta_list, phi, 
   ridge_phi, binding_null_basis, ridge_beta, K, Fmatrix, adaptive_weight, weights,
   weights_W,
   rank_F, group_mapping, E.z, tol_epsilon, ll_type, fisher = FALSE, oelker = FALSE){
  
  
  args <- as.list(match.call())
  args <- args[names(args) != '']
  args$ll_type <- NULL
  args$ll_type <- 'all'
  #Delta^2 L(\theta | y) or J
  louis_SE <- do.call('standard_error_louis', args)$vcov
  
  k_range <- setdiff(1:K, 1)
  
  pi_ik <- softmax_matrix(W %*% t(phi))
  norm_weights <- weights_W/sum(weights_W)
  pi_bar <- colSums((Diagonal(x = norm_weights) %*% pi_ik))
  
  single_intercept <- !is.null(mu)
  
  if (single_intercept){
    split_p <- mapply(X_list, beta_list, SIMPLIFY = FALSE, FUN=function(x_k, b_k){plogis(mu + as.vector(x_k %*% b_k))})
  }else{
    split_p <- mapply(X_list, beta_list, SIMPLIFY = FALSE, FUN=function(x_k, b_k){plogis(as.vector(x_k %*% b_k))})
  }
  
  #Get the outer product of the score
  score_all_beta <- do.call('cbind', lapply(1:K, FUN=function(k){
    Diagonal(x = E.z[,k]) %*% (t(group_mapping) %*% Diagonal(x = weights * (y - split_p[[k]])) %*% X_list[[k]])
  }))
  score_all_phi <- do.call('cbind', lapply(k_range, FUN=function(k){
    Diagonal(x = weights_W * (E.z[,k] - pi_ik[,k])) %*% W
  }))
  
  if (single_intercept){
    score_all_mu <- rowSums(sapply(1:K, FUN=function(k){
      as.vector(Diagonal(x = E.z[,k]) %*% t(group_mapping) %*% (weights * (y - split_p[[k]])))
    }))
    score_all <- cbind(score_all_mu, score_all_beta, score_all_phi)  
  }else{
    score_all <- cbind(score_all_beta, score_all_phi)  
  }
  
  #Outer Product of the Score or "K"
  meilijson_hessian <- crossprod(score_all) - tcrossprod(colSums(score_all)) * 1/nrow(score_all)
  # meilijson_hessian <- var(as.matrix(score_all))

  vcov_sandwich <- louis_SE %*% meilijson_hessian %*% louis_SE
  
  return(
    list('vcov' = vcov_sandwich,
         'louis_SE' = louis_SE,
         'meilijson_hessian' = meilijson_hessian)
  )
}


project_standard_errors <- function(vcov_obj, flat_estimates, binding_null_basis, basis_M, single_intercept){
  
  K <- length(binding_null_basis)
  beta_size <- sapply(binding_null_basis, ncol)
  beta_cumsize <- c(0, cumsum(beta_size))
  

  
  if (single_intercept){
    all_proj <- matrix(nrow = 0, ncol = ncol(vcov_obj))
    for (k in 1:K){
      basis_k <- binding_null_basis[[k]]
      empty_matrix <- matrix(0, nrow = nrow(basis_k) + 1, ncol = ncol(vcov_obj))  
      empty_matrix[nrow(empty_matrix), 1] <- 1
      if (ncol(basis_k) == 0){
      }else if (ncol(basis_k) == 1){
        empty_matrix[-nrow(empty_matrix),1 + beta_cumsize[k] + 1:ncol(basis_k)] <- as.vector(basis_k)
      }else{
        empty_matrix[-nrow(empty_matrix),1 + beta_cumsize[k] + 1:ncol(basis_k)] <- as.matrix(basis_k)
      }
      all_proj <- rbind(all_proj, empty_matrix)
    }
    phi_size <- ncol(vcov_obj) - (1 + beta_cumsize[K+1])
    all_proj <- rbind(all_proj,
      cbind(matrix(0, nrow = phi_size, ncol = 1 + beta_cumsize[K+1]), 
            Diagonal(n = phi_size))
    )
    
  }else{
    all_proj <- bdiag(binding_null_basis)
    if (K == 1){
      phi_size <- 0
    }else{
      phi_size <- ncol(vcov_obj) - (beta_cumsize[K+1])
      all_proj <- cbind(all_proj, sparseMatrix(i=1,j=1,x=0,dims=c(nrow(all_proj), phi_size)))
      all_proj <- rbind(all_proj,
                        cbind(matrix(0, nrow = phi_size, ncol = beta_cumsize[K+1]), Diagonal(n = phi_size))
      )
    }
    all_proj <- drop0(all_proj)  
  }
  
  proj_vcov <- all_proj %*% vcov_obj %*% t(all_proj)
  proj_to_original <- bdiag(bdiag(lapply(1:K, FUN=function(i){basis_M})), Diagonal(n = phi_size))
  
  
  proj_original <- proj_to_original %*% proj_vcov %*% t(proj_to_original)
  return(list('original' = proj_original, 'estimates' = proj_to_original %*% all_proj %*% flat_estimates, 'unconstrained' = proj_vcov))
}

estimate.vcov <- function(object, tol_epsilon = 1e-4, type = 'louis', format_se = FALSE,
                          return_raw_se = FALSE){

  if (inherits(object, 'FactorHet')){
    X <- object$internal_parameters$data$X
    W <- object$internal_parameters$data$W
    y <- object$internal_parameters$data$y
    phi <- coef(object, 'phi')
    basis_M <- object$internal_parameters$data$basis_M
    beta <- object$parameters$nullspace_beta
    K <- ncol(beta)
    lambda <- object$parameters$eff_lambda
    gamma <- object$parameters$gamma
    Fmatrix <- object$internal_parameters$data$Fmatrix
    rank_F <- object$internal_parameters$data$rank_F
    adaptive_weight <- object$internal_parameters$data$adaptive_weight
    ridge_phi <- 1/object$internal_parameters$control$prior_var_phi
    group_mapping <- object$internal_parameters$group$group_mapping
    group_mapping <- sparseMatrix(i = 1:nrow(X), 
      j = match(object$internal_parameters$data$group, 
      object$internal_parameters$group$unique_groups), x = 1)
    factor_levels <- object$internal_parameters$factor_levels
    tau_truncate <- object$internal_parameters$control$tau_truncate
    # tau_equality_threshold <- object$internal_parameters$control$tau_equality_threshold
    clip_tiny <- object$internal_parameters$misc$clip_tiny
    recons_beta <- coef(object)
    term_position <- object$internal_parameters$penalty$term_position
    coef_names <- object$internal_parameters$penalty$coef
    single_intercept <- object$internal_parameters$single_intercept
    
    weights <- object$internal_parameters$weights$weights
    weights_W <- object$internal_parameters$weights$weights_W
    
    log_method <- object$internal_parameters$control$log_method
  }else{stop('')}
  
  if (lambda == 0){
    ridge_beta <- 1/object$internal_parameters$control$prior_var_beta
    ridge_beta <- c(0, rep(ridge_beta, nrow(basis_M) - 1))
    
    ridge_beta <- t(basis_M) %*% sparse_diag(ridge_beta) %*% basis_M
    ridge_beta <- drop0(zapsmall(ridge_beta, clip_tiny))
    ridge_beta <- as(ridge_beta, 'dgCMatrix')
    ridge_beta <- lapply(1:K, FUN=function(i){ridge_beta})
  }else{
    ridge_beta <- NULL
  }
  
  pi_bar <- colSums(Diagonal(x = weights_W/sum(weights_W)) %*% softmax_matrix(W %*% t(phi)))
  
  E.tau <- mapply(list_from_cols(beta), pi_bar^gamma, 
    adaptive_weight, SIMPLIFY = FALSE, FUN=function(b, pi.g, aw){
      F_EM_update(beta = b, Fmat = Fmatrix, lambda = lambda, pi.gamma = pi.g, exo_weights = aw)
  })
  
  if (lambda > 0){
    
    binding_restrictions <- lapply(E.tau, FUN=function(etau.k){
      lapply(etau.k, FUN=function(k){which(k[,1] >= tau_truncate)})
    })

  }else{
    binding_restrictions <- NULL
  }
  
  # binding_lookup <- build_binding_lookup(Fmatrix = Fmatrix, factor_levels = factor_levels, term_position = term_position, coef_names = coef_names)
  # 
  # binding_restrictions <- apply(recons_beta, MARGIN = 2, FUN=function(b){
  #   br <- calc_distance_binding(b, binding_lookup)
  #   # br <- lapply(prepare_fusion(factor_levels = factor_levels, term_position = term_position, coef_names = coef_names, beta = b, simplify = FALSE), FUN=function(i){
  #   #   sapply(i, FUN=function(j){max(unlist(j))})
  #   # })
  #   br <- lapply(br, FUN=function(j){which(j <= tau_equality_threshold)})
  #   return(br)
  # })
  
  K <- ncol(beta)

  loglik.k <- plogis(as.matrix(X %*% beta) * (2 * y - 1), log.p = TRUE) 
  
  E.z <- calculate_posterior_zi(loglik.k = loglik.k, 
    group_mapping = group_mapping, K = K, W = W, ncol_W = ncol(W), 
    pi = NA, phi = phi)
  
  copy_raw_beta <- beta
  copy_phi <- phi
  copy_pi <- softmax_matrix(W %*% t(copy_phi))
  rm(beta, phi)
  if (single_intercept){
    copy_mu <- as.numeric(copy_raw_beta[nrow(copy_raw_beta),1])
    names(copy_mu) <- NULL
    copy_beta <- copy_raw_beta[-nrow(copy_raw_beta),,drop=F]
  }else{
    copy_beta <- copy_raw_beta
    copy_mu <- NULL
  }
  
  if (lambda == 0){
    if (single_intercept){
      binding_null_basis <- lapply(1:K, FUN=function(k){sparse_diag(rep(1, ncol(X) - 1))})
      X_list <- lapply(binding_null_basis, FUN=function(i){X[, -ncol(X)] %*% i})
    }else{
      binding_null_basis <- lapply(1:K, FUN=function(k){sparse_diag(rep(1, ncol(X)))})
      X_list <- lapply(binding_null_basis, FUN=function(i){X %*% i})
    }
  }else{
    if (single_intercept){
      binding_null_basis <- lapply(binding_restrictions, FUN=function(r_k){
        binding_k <- do.call('rbind', mapply(r_k, Fmatrix, SIMPLIFY = FALSE, FUN=function(i, F_j){
          do.call('rbind', F_j[i])
        }))
        if (is.null(binding_k)){
          return(sparse_diag(rep(1, ncol(X)-1)))
        }else{
          binding_k <- binding_k[,-ncol(binding_k)]
          out_k <- drop0(zapsmall(calculate_nullspace_basis(binding_k), clip_tiny))
          return(out_k)
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
  
  beta_list <- mapply(binding_null_basis, 1:K, FUN=function(i,k){
    solve(crossprod(i), t(i) %*% copy_beta[,k,drop=F])
  })
  
  beta_size <- sapply(beta_list, nrow)
  
  flat_beta <- do.call('c', lapply(beta_list, as.vector))
  flat_par <- c(copy_mu, flat_beta, as.vector(t(copy_phi[-1,])))
  mu <- copy_mu
  phi <- copy_phi
  
  
  
  if (type == 'louis'){
    
    analytical <- standard_error_louis(beta_list = beta_list, mu = mu, 
     group_mapping = group_mapping, binding_null_basis = binding_null_basis,
     phi = phi, X_list = X_list, W = W, y = y, lambda = lambda, gamma = gamma, 
     ridge_phi = ridge_phi, rank_F = rank_F, Fmatrix = Fmatrix, K = K, ll_type = 'all',
     tol_epsilon = tol_epsilon, adaptive_weight = adaptive_weight, E.z = E.z,
     ridge_beta = ridge_beta, weights_W = weights_W, weights = weights)
    
  }else if (type == 'sandwich'){
    
    analytical <- standard_error_sandwich(beta_list = beta_list, mu = mu, 
      group_mapping = group_mapping, binding_null_basis = binding_null_basis,
      phi = phi, X_list = X_list, W = W, y = y, lambda = lambda, gamma = gamma, 
      ridge_phi = ridge_phi, rank_F = rank_F, Fmatrix = Fmatrix, K = K, ll_type = 'all',
      tol_epsilon = tol_epsilon, adaptive_weight = adaptive_weight, E.z = E.z,
      ridge_beta = ridge_beta, weights_W = weights_W, weights = weights)
    
  }else{
    stop('...')
  }
  if (return_raw_se){
    return(analytical)
  }

  SE <- project_standard_errors(vcov_obj = analytical$vcov, flat_estimates = flat_par, 
    binding_null_basis = binding_null_basis, basis_M = basis_M, single_intercept = !is.null(mu))
  k_range <- setdiff(1:K, 1)
  
  #
  # Get back to ONLY the main effects (i.e. ignoring over-parameterization)
  #
  
  if (grepl(log_method, pattern='^log')){
    
    naug <- nrow(SE$estimates)
    if (K > 1){
      extract_phi <- seq(1 + naug - ncol(copy_phi) * (K-1), naug)
    }else{
      extract_phi <- c()
    }
    extract_beta <- do.call('c', lapply(1:K, FUN=function(k){
      nrow(basis_M) * (k-1) + 1:length(coef_names)
    }))
    
    extract_main <- c(extract_beta, extract_phi)
    SE$estimates <- SE$estimates[extract_main,,drop=F]
    SE$original <- SE$original[extract_main, extract_main]
  }else if (log_method == 'standard'){
    
  }else{stop('Check SE for log_method')}
  
  checksum_SE <- abs(max(c(as.vector(coef(object)), as.vector(t(coef(object, 'phi')[-1,]))) +
                           -as.vector(SE$estimates)))
  if (checksum_SE > 1e-4){
    print(checksum_SE)
    warning('Alignment error after re-projection')
  }
  #Format or not SE
  if (format_se){
    proj_SE <- data.frame(est = as.vector(SE$estimates), 
      se = sqrt(diag(SE$original)),
      manual_estimates = c(as.vector(coef(object)), 
                           as.vector(t(coef(object, 'phi')[-1,]))))
    proj_SE$id <- 1:nrow(proj_SE)
    proj_SE$type <- c(as.vector(sapply(1:K, FUN=function(i){paste0('beta', i, '_', rownames(coef(object)))})),
                      as.vector(sapply(k_range, FUN=function(i){paste0('phi', i, '_', colnames(coef(object,'phi')))})))
    
    
    return(list('formatted_SE' = proj_SE))
  }else{
    vcov_names <- c(as.vector(sapply(1:K, FUN=function(i){paste0('beta', i, '_', rownames(coef(object)))})),
      as.vector(sapply(k_range, FUN=function(i){paste0('phi', i, '_', colnames(coef(object,'phi')))})))
    rownames(SE$original) <- colnames(SE$original) <- vcov_names
    return(list('vcov' = SE$original, 'dim' = dim(coef(object)), 
                'est' = SE$estimates))
  }
}

#' @export
vcov.FactorHet <- function(object, phi = TRUE, se.method = NULL, ...){
  
  options <- list(...)
  if(length(options) != 0){
    stop('Only argumnets to vcov for FactorHet are "phi" and "se.method".')
  }
  if (is.null(se.method)){
    est_vcov <- object$vcov
    if (is.null(est_vcov)){
      est_vcov <- tryCatch(estimate.vcov(object), error = function(e){NULL})
    }
  }else{
    est_vcov <- estimate.vcov(object, type = se.method)
  }
  if (is.null(est_vcov)){return(NULL)}
  dim_beta <- est_vcov$dim
  dim_beta <- prod(dim_beta)
  est_vcov <- est_vcov$vcov
  if (!phi){#Return *only* the standard errors on beta
    est_vcov <- est_vcov[1:dim_beta, 1:dim_beta]
  }
  return(as.matrix(est_vcov))
}