context('Conjugate Gradient')

test_that('Conjugate Gradient methods agree for simple case.', {
  #Simple base R for conjugate gradient
  N <- 5000
  p <- 20
  
  X <- rsparsematrix(nrow = N, ncol = p, density = 0.01)
  #Ensure ridge is PD to make problem well-defined.
  ridge <- rsparsematrix(nrow = p, ncol = p, density = 0.01, rand.x = runif, symmetric = TRUE)
  ridge <- t(ridge) %*% ridge + Diagonal(n = p)
  ridge <- as(ridge, 'dgCMatrix')
  omega <- rexp(N)
  y <- rnorm(N)
  
  direct_solve <- solve(t(X) %*% Diagonal(x = omega) %*% X + ridge, t(X) %*% y)
  
  base_R_CG <- function(X, y, omega, ridge, old_beta, tol, it_max = 0){
    
    if (it_max == 0){
      it_max <- ncol(X)
    }
    adj_y <- y / sqrt(omega)
    adj_X <- Diagonal(x = sqrt(omega)) %*% X
    
    precond_diag <- colSums(adj_X^2) + diag(ridge)
    precond_diag <- Diagonal(x = 1/precond_diag)
    
    t_adj_X <- t(adj_X)
    
    beta <- old_beta
    
    residual <- adj_y - adj_X %*% beta
    rhsNorm2 <- as.numeric(crossprod(t_adj_X %*% adj_y))
    
    normal_residual <- t_adj_X %*% residual - ridge %*% beta
    
    p <- precond_diag %*% normal_residual
    
    absNew <- t(normal_residual) %*% p
    
    it <- 0
    convg <- FALSE
    
    tol <- sqrt(.Machine$double.eps)
    
    threshold = tol^2 * rhsNorm2
    
    adenom <- anew <- step_seq <- alpha_seq <- rep(NA, it_max)
    
    while (it < it_max & !convg){
      tmp <- adj_X %*% p
      alpha_denom <- (crossprod(tmp) + t(p) %*% ridge %*% p)
      
      alpha <- as.numeric(absNew / alpha_denom)
      beta <- beta + alpha * p
      residual <- residual - alpha * tmp
      normal_residual <- t_adj_X %*% residual - ridge %*% beta
      
      adenom[it] <- alpha_denom
      anew[it] <- absNew
      
      residualNorm2 <- sum(normal_residual^2)
      tol_error <- sqrt(residualNorm2 / rhsNorm2)
      
      if (residualNorm2 < threshold){
        break
        
      }
      
      z <- precond_diag %*% normal_residual
      absOld <- absNew
      absNew <- t(normal_residual) %*% z
      step_beta <- absNew/absOld
      p <- z + as.numeric(step_beta) * p
      
      alpha_seq[it] <- alpha
      step_seq[it] <- step_beta
      
      it <- it  + 1
      
    }
    
    return(list(it = it, tol_error = tol_error, beta = beta, step_seq = step_seq, alpha_seq = alpha_seq, anew = anew, adenom = adenom))
  }
  
  base_CG <- base_R_CG(X = X, y = y, omega = omega, ridge = ridge, old_beta = matrix(0, ncol(X)), tol = .sqrt(.Machine$double.eps))
  custom_CG <- cg_custom(X = X, omega = omega, list_ridge = list(ridge), s = y, old_beta = matrix(0, ncol(X)), 
                         weights = matrix(0, nrow = 1, ncol = 0),
                         tol = sqrt(.Machine$double.eps), K = 1)

  expect_equal(as.vector(custom_CG$beta), as.vector(base_CG$beta))
  expect_equal(custom_CG$error, base_CG$tol_error)
  expect_equal(log(custom_CG$error), log(base_CG$tol_error), tolerance = 0.1)
  
  # Should be very close to direct solve
  expect_equal(as.vector(direct_solve), as.vector(base_CG$beta), tolerance = 1e-6)
  expect_equal(as.vector(direct_solve), as.vector(custom_CG$beta), tolerance = 1e-6)
  
})

test_that('Conjugate Gradient (Custom) agree for harder case where LeastSquaresConjugateGradient fails', {
  #Simple base R for conjugate gradient
  N <- 5000
  p <- 1000
  
  X <- rsparsematrix(nrow = N, ncol = p, density = 0.01)
  #Ensure ridge is PD to make problem well-defined.
  ridge <- rsparsematrix(nrow = p, ncol = p, density = 0.01, rand.x = runif, symmetric = TRUE)
  ridge <- t(ridge) %*% ridge + Diagonal(n = p)
  ridge <- as(ridge, 'dgCMatrix')
  omega <- rexp(N)
  y <- rnorm(N)
  
  direct_solve <- solve(t(X) %*% Diagonal(x = omega) %*% X + ridge, t(X) %*% y)
  
  base_R_CG <- function(X, y, omega, ridge, old_beta, tol, it_max = 0){
    
    if (it_max == 0){
      it_max <- ncol(X)
    }
    adj_y <- y / sqrt(omega)
    adj_X <- Diagonal(x = sqrt(omega)) %*% X
    
    precond_diag <- colSums(adj_X^2) + diag(ridge)
    precond_diag <- Diagonal(x = 1/precond_diag)
    
    t_adj_X <- t(adj_X)
    
    beta <- old_beta
    
    residual <- adj_y - adj_X %*% beta
    rhsNorm2 <- as.numeric(crossprod(t_adj_X %*% adj_y))
    
    normal_residual <- t_adj_X %*% residual - ridge %*% beta
    
    p <- precond_diag %*% normal_residual
    
    absNew <- t(normal_residual) %*% p
    
    it <- 0
    convg <- FALSE
    
    tol <- sqrt(.Machine$double.eps)
    
    threshold = tol^2 * rhsNorm2
    
    adenom <- anew <- step_seq <- alpha_seq <- rep(NA, it_max)
    
    while (it < it_max & !convg){
      tmp <- adj_X %*% p
      alpha_denom <- (crossprod(tmp) + t(p) %*% ridge %*% p)
      
      alpha <- as.numeric(absNew / alpha_denom)
      beta <- beta + alpha * p
      residual <- residual - alpha * tmp
      normal_residual <- t_adj_X %*% residual - ridge %*% beta
      
      adenom[it] <- alpha_denom
      anew[it] <- absNew
      
      residualNorm2 <- sum(normal_residual^2)
      tol_error <- sqrt(residualNorm2 / rhsNorm2)
      
      if (residualNorm2 < threshold){
        break
        
      }
      
      z <- precond_diag %*% normal_residual
      absOld <- absNew
      absNew <- t(normal_residual) %*% z
      step_beta <- absNew/absOld
      p <- z + as.numeric(step_beta) * p
      
      alpha_seq[it] <- alpha
      step_seq[it] <- step_beta
      
      it <- it  + 1
      
    }
    
    return(list(it = it, tol_error = tol_error, beta = beta, step_seq = step_seq, alpha_seq = alpha_seq, anew = anew, adenom = adenom))
  }
  
  base_CG <- base_R_CG(X = X, y = y, omega = omega, ridge = ridge, old_beta = matrix(0, ncol(X)), tol = .sqrt(.Machine$double.eps))
  custom_CG <- cg_custom(X = X, omega = omega, list_ridge = list(ridge), s = y, old_beta = matrix(0, ncol(X)), weights = matrix(0, nrow = 1, ncol = 0),
                         tol = sqrt(.Machine$double.eps), K = 1)
  # CG methods should agree:
  expect_equal(as.vector(base_CG$beta), as.vector(custom_CG$beta))
  expect_equal(base_CG$tol_error, custom_CG$error)
  expect_equal(log(base_CG$tol_error), log(custom_CG$error), tolerance = 0.1)
  

  # Should be very close to direct solve
  expect_equal(as.vector(direct_solve), as.vector(base_CG$beta), tolerance = 1e-6)
  expect_equal(as.vector(direct_solve), as.vector(custom_CG$beta), tolerance = 1e-6)

})


test_that('Conjugate Gradient methods agree for K > 1', {
  #Simple base R for conjugate gradient
  N <- 5000
  p <- 20
  K <- 2
  
  X <- rsparsematrix(nrow = N, ncol = p, density = 0.01)
  #Ensure ridge is PD to make problem well-defined.
  list_ridge <- list()
  omega <- matrix(NA, nrow = N, ncol = K)
  for (k in 1:K){
    ridge <- rsparsematrix(nrow = p, ncol = p, density = 0.01, rand.x = runif, symmetric = TRUE) + Diagonal(n = p)
    ridge <- t(ridge) %*% ridge + Diagonal(n = p)
    ridge <- as(ridge, 'dgCMatrix')
    list_ridge[[k]] <- ridge  
    omega[,k] <- rexp(N)
  }
  y <- rnorm(N)
  
  direct_solve <- matrix(NA, nrow = p, ncol = K)
  for (k in 1:K){
    direct_solve[,k] <- as.vector(solve(t(X) %*% Diagonal(x = omega[,k]) %*% X + list_ridge[[k]], t(X) %*% y))
  }
  custom_CG <- cg_custom(K = 2, X = X, omega = omega, list_ridge = list_ridge, s = y, weights = matrix(0, nrow = 1, ncol = 0),
                         old_beta = matrix(0, nrow = p, ncol = K), tol = sqrt(.Machine$double.eps))
  
  expect_equal(as.matrix(custom_CG$beta), as.matrix(direct_solve))
  
  #Doing it one-by-one should also agree
  
  loop_CG <- matrix(NA, nrow = p, ncol = K)
  for (k in 1:K){
    a <- cg_custom(K = 1, X = X, omega = omega[,k],
                   list_ridge = list_ridge[k], s = y,
                   weights = matrix(0, nrow = 1, ncol = 0),
                   old_beta = matrix(0, nrow = p, ncol = 1),
                   tol = sqrt(.Machine$double.eps))
    loop_CG[,k] <- a$beta
  }
  expect_equal(as.matrix(custom_CG$beta), as.matrix(loop_CG))
  expect_equal(as.matrix(loop_CG), as.matrix(direct_solve))
  
})

