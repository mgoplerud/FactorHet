if (isTRUE(as.logical(Sys.getenv("CI")))){
  # If on CI
  NITER <- 2
  env_test <- "CI"
}else if (!identical(Sys.getenv("NOT_CRAN"), "true")){
  # If on CRAN
  NITER <- 2
  env_test <- "CRAN"
  set.seed(5)
}else{
  # If on local machine
  NITER <- 2000
  env_test <- 'local'
}

test_that('K = 1 and lambda = 0 agrees with ridge ', {
  N <- 1000
  K <- 1
  dta <- data.frame(
    state = sample(state.name[1:4], N, replace = T),
    letter = sample(letters[1:3], N, replace = T)
  )
  dta$mod <- runif(nrow(dta), -1, 1)
  dta$mod2 <- sample(LETTERS[1:5], N, replace = T)
  dta$y1 <- rbinom(nrow(dta), 1, plogis(runif(5, -1 , 1)[match(dta$state, state.name[1:5])] + runif(5, -1, 1)[match(dta$letter, letters)]))
  dta$y2 <- rbinom(nrow(dta), 1, plogis(runif(5, -1 , 1)[match(dta$state, state.name[1:5])] + runif(5, -1, 1)[match(dta$letter, letters)]))
  dta$mod3 <- dta$mod^2
  dta$y <- ifelse(rbinom(nrow(dta), 1, plogis(dta$mod)) == 1, dta$y1, dta$y2)
  
  #Does it agree with standard LASSO estimates
  est_simple <- FactorHet(formula = y ~ state + letter, design = dta, K = K, 
                          lambda = 0, moderator = ~ mod + mod2 + mod3,
                          control = FactorHet_control(return_data = TRUE))  
  pred_simple <- predict(est_simple)
  X <- est_simple$internal_parameters$data$X

  vcov_cons <- solve(t(X) %*% Diagonal(x = pred_simple * (1-pred_simple)) %*% X)
  
  basis_M <- est_simple$internal_parameters$data$basis_M
  vcov_analytic <- basis_M %*% vcov_cons %*% t(basis_M)
  expect_equivalent(vcov(est_simple), as.matrix(vcov_analytic))
  
  #Does it agree with Hessian of log-posterior with ridge?
  est_simple <- FactorHet(formula = y ~ state + letter, design = dta, K = K, 
                          lambda = 0, moderator = ~ mod + mod2 + mod3,
                          control = FactorHet_control(prior_var_beta = 1/3, 
                                                      log_method = 'standard',
                                                      return_data = TRUE))  
  pred_simple <- predict(est_simple)
  
  vcov_cons <- solve(t(X) %*% Diagonal(x = pred_simple * (1-pred_simple)) %*% X +
    t(basis_M) %*% sparse_diag(c(0,rep(3, nrow(basis_M) - 1))) %*% (basis_M))
  vcov_analytic <- basis_M %*% vcov_cons %*% t(basis_M)
  expect_equivalent(vcov(est_simple), as.matrix(vcov_analytic))
  
})


test_that('K = 3 agrees with manual differentation', {
  
  skip_on_cran()
  skip_on_ci()
  
  N <- 1000
  K <- 3
  dta <- data.frame(
    state = sample(state.name[1:4], N, replace = T),
    letter = sample(letters[1:3], N, replace = T)
  )
  
  dta$group <- rep(sample(1:250, 500, replace = T), each = 2)
  dta$task <- rep(1:500, each = 2)
  dta$prof <- as.vector(sapply(1:500, FUN=function(i){c('l', 'r')[sample(2)]}))
  dta$y <- as.vector(sapply(1:500, FUN=function(i){sample(0:1)}))
  dta$mod <- runif(500, -1, 1)[dta$group]
  dta$mod2 <- sample(LETTERS[1:5], 500, replace = T)[dta$group]
  
  p1 <- plogis(runif(5, -1 , 1)[match(dta$state, state.name[1:5])] + runif(5, -1, 1)[match(dta$letter, letters)])
  p2 <- plogis(runif(5, -1 , 1)[match(dta$state, state.name[1:5])] + runif(5, -1, 1)[match(dta$letter, letters)]) 
  
  dta$y1 <- unlist(lapply(split(p1, rep(1:500, each = 2)), FUN=function(i){sample(0:1, 2, prob = i)}))
  dta$y2 <- unlist(lapply(split(p2, rep(1:500, each = 2)), FUN=function(i){sample(0:1, 2, prob = i)}))

  dta$placeholder <- unsplit(lapply(split(dta$mod, dta$group), FUN=function(i){rep(rbinom(1, 1, plogis(i[1])), length(i))}), dta$group)
  dta$y <- ifelse(dta$placeholder == 1, dta$y1, dta$y2)


  est_simple <- FactorHet(formula = y ~ state + letter, design = dta, K = K, 
                          lambda = 1e-3, moderator = ~ mod + mod2,
                          group = ~ group, task = ~ task, choice_order = ~ prof,                          
                          control = FactorHet_control(
                            prior_var_beta = Inf, single_intercept = TRUE,
                            log_method = 'standard',
                            return_data = TRUE))  
  if (TRUE){
   
    jacobian <- getFromNamespace('jacobian', 'numDeriv')
    hessian <- getFromNamespace('hessian', 'numDeriv')
    
    group_mapping <- sparseMatrix(i = 1:nrow(est_simple$internal_parameters$data$X), 
                                  j = match(est_simple$internal_parameters$group$group_mapping, est_simple$internal_parameters$group$unique_groups), x = 1)
    
    #beta, pi, phi
    check_args <- list(
      X = est_simple$internal_parameters$data$X,
      y = est_simple$internal_parameters$data$y,
      W = est_simple$internal_parameters$data$W,
      weights_W = est_simple$internal_parameters$data$weights_W,
      group_mapping = group_mapping,
      Fmatrix = est_simple$internal_parameters$data$Fmatrix,
      E.prob = NA,
      ridge_phi = 1/est_simple$internal_parameters$control$prior_var_phi,
      ridge_beta = 0,
      adaptive_weight = est_simple$internal_parameters$adaptive_weight,
      gamma = est_simple$parameters$gamma,
      lambda = est_simple$parameters$eff_lambda,
      separate = TRUE
    )
    
    flat_beta <- est_simple$parameters$nullspace_beta
    flat_phi <- est_simple$parameters$phi
    if (K == 1){
      flat_phi <- matrix(nrow = 0, ncol = 0)
    }
    flat_par <- safe_flat_par <- c(
      flat_beta[nrow(flat_beta), 1], 
      as.vector(flat_beta[-nrow(flat_beta),]),
      as.vector(t(flat_phi[-1,,drop=F]))
    )
    
    counter <- 0
    check_args$kern_epsilon <- 1e-4
    ll_numderiv <- function(par, K, other_args, pos = 1){
      counter <<- counter + 1
      if (counter %% 500 == 0){message(counter, appendLF = FALSE)}
      mu <- par[1]
      par <- par[-1]
      beta <- matrix(par[seq_len(K * (ncol(other_args$X) - 1))], ncol = K)
      beta <- rbind(beta, mu)
      if (K != 1){
        phi <- matrix(par[-seq_len(K * (ncol(other_args$X) - 1))], byrow = T,
                      nrow = (K-1))
        phi <- rbind(0, phi)
        pi <- colMeans(softmax_matrix(other_args$W %*% t(phi)))
      }else{
        phi <- NA
        pi <- 1
      }
      return(do.call('evaluate_loglik', 
                     c(list(beta = beta, phi = phi, pi = pi), other_args))[pos])
    }  
    
    direct_optim <- optim(fn = ll_numderiv, par = flat_par, K = K, 
                          control = list(fnscale = -1, reltol = 0, abstol = 0), method = 'BFGS',
                          other_args = check_args)
    
    flat_par <- direct_optim$par
    num_grad <- as.vector(jacobian(func = ll_numderiv, 
                                   x = flat_par, K = K, other_args = check_args))
    num_hessian <- optimHess(fn = ll_numderiv, 
                             par = flat_par, K = K, other_args = check_args)
    num_hessian_numd <- hessian(func = ll_numderiv, 
                                x = flat_par, K = K, other_args = check_args)
    # Check gradient is basically stationary
    expect_lte(max(abs(num_grad)), 1e-4)
    
    recons_par <- direct_optim$par
    recons_mu <- recons_par[1]
    recons_par <- recons_par[-1]
    recons_beta <- matrix(recons_par[seq_len(K * (ncol(check_args$X) - 1))], ncol = K)
    recons_beta <- rbind(recons_beta, recons_mu)
    if (K != 1){
      recons_phi <- matrix(recons_par[-seq_len(K * (ncol(check_args$X) - 1))], byrow = TRUE, nrow = (K-1))
      recons_phi <- rbind(0, recons_phi)
      recons_pi_bar <- colMeans(softmax_matrix(check_args$W %*% t(recons_phi)))
    }else{
      recons_phi <- matrix(0, nrow = 1, ncol = 1)
      recons_pi_bar <- 1
    }
    
    copy_est <- est_simple
    copy_est$parameters$beta <- NA
    copy_est$parameters$pi <- recons_pi_bar
    copy_est$parameters$phi <- recons_phi
    copy_est$parameters$nullspace_beta <- recons_beta
    
    vcov_analytical <- estimate.vcov(copy_est, return_raw_se = TRUE)
    vcov_numerical <-  -solve(num_hessian)
    vcov_numerical_numd <- -solve(num_hessian_numd)
    
    # Check alignment in correlation and absolute terms
    cor_vcov <- cor(as.vector(vcov_analytical$vcov), as.vector(vcov_numerical))
    expect_gte(cor_vcov, 0.999)
    test_lm <- lm(as.vector(vcov_analytical$vcov) ~ as.vector(vcov_numerical))
    expect_equivalent(coef(test_lm), c(0, 1), tol = 1e-4, scale = 1)
    # Get the variance and the percent difference
    var_analytical <- diag(vcov_analytical$vcov)
    var_numerical <- diag(vcov_numerical)
    range_diff_var <- range(100 * abs(var_numerical - var_analytical)/var_numerical)
    expect_true(all(range_diff_var > 0))
    # Within 1 percent.
    expect_lte(range_diff_var[2], 1)

  }
})
