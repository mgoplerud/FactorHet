context('Test moderators')

test_that('Confim working post regression functions (moderator)', {

  dta <- data.frame(
    state = sample(state.name[1:4], 100, replace = T),
    letter = sample(letters[1:3], 100, replace = T)
  )
  
  dta$mod <- runif(nrow(dta), -1, 1)
  dta$mod2 <- sample(letters[1:5], nrow(dta), replace = T)
  dta$y1 <- rbinom(nrow(dta), 1, plogis(runif(5, -1 , 1)[match(dta$state, state.name[1:5])] + runif(5, -1, 1)[match(dta$letter, letters)]))
  dta$y2 <- rbinom(nrow(dta), 1, plogis(runif(5, -1 , 1)[match(dta$state, state.name[1:5])] + runif(5, -1, 1)[match(dta$letter, letters)]))
  
  dta$y <- ifelse(rbinom(nrow(dta), 1, plogis(dta$mod)) == 1, dta$y1, dta$y2)
  
  est_simple <- FactorHet(formula = y ~ state + letter, design = dta, 
     moderator = ~ mod + mod2, 
     initialize = FactorHet_init(short_EM_it = 5),
     control = FactorHet_control(
       prior_var_phi = 4, prior_var_beta = 5,
       iterations = 15),
     K = 2, lambda = 0)
  
  expect_gte(min(diff(logLik(est_simple, 'log_posterior_seq'))), -sqrt(.Machine$double.eps))
  
  est_mod_eff <- tryCatch(moderator_AME(est_simple), error = function(e){NULL})
  est_mod_post <- tryCatch(posterior_FactorHet(est_simple), error = function(e){NULL})
  
  expect_false(is.null(est_mod_eff))  
  expect_false(is.null(est_mod_post)) 
  
  expect_equivalent(est_mod_post$compare$post.predict, 
                    est_simple$posterior$posterior_predictive$cluster_1)
  expect_equivalent(est_mod_post$compare$posterior, 
                    est_simple$posterior$posterior$cluster_1)
  
})

test_that("Gradient of Moderator Stage is Zero at Optimum", {

  skip_on_cran()
  
  N <- 100
  K <- 3
  dta <- data.frame(
    state = sample(state.name[1:4], N, replace = T),
    letter = sample(letters[1:3], N, replace = T)
  )
  dta$mod <- runif(nrow(dta), -1, 1)
  dta$mod2 <- sample(LETTERS[1:5], N, replace = T)
  dta$y1 <- rbinom(nrow(dta), 1, plogis(runif(5, -1 , 1)[match(dta$state, state.name[1:5])] + runif(5, -1, 1)[match(dta$letter, letters)]))
  dta$y2 <- rbinom(nrow(dta), 1, plogis(runif(5, -1 , 1)[match(dta$state, state.name[1:5])] + runif(5, -1, 1)[match(dta$letter, letters)]))
  
  dta$y <- ifelse(rbinom(nrow(dta), 1, plogis(dta$mod)) == 1, dta$y1, dta$y2)
  
  est_simple <- FactorHet(formula = y ~ state + letter, design = dta, K = K, 
    lambda = 1e-2, moderator = ~ mod + mod2, 
    control =  FactorHet_control(
      log_method = 'standard',
      tolerance.parameters = 1e-8, 
      tolerance.logposterior = 0, 
      return_data = TRUE))  

  expect_gte(min(diff(logLik(est_simple, 'log_posterior_seq'))), 
             - sqrt(.Machine$double.eps))
  
  data <- est_simple$internal_parameters$data
  data$norm_weights <- as.vector(data$weights_W/sum(data$weights_W))
  
  expect_equal(
    cpp_obj_phi(par = as.vector(coef(est_simple, 'phi')[-1,,drop=F]),
              K = K, W = data$W, norm_weights = data$norm_weights,
              weights_W = data$weights_W,
              group_E_prob = as.matrix(est_simple$posterior$posterior[,-1]), 
              ridge_penalty = make_TMatrix(K) * 1/est_simple$internal_parameters$control$prior_var_phi,
              gamma = 1, rank_F = data$rank_F,
              power_pi = 1, b_r = data$b_r,
              lambda = est_simple$parameters$eff_lambda),
    objective_moderator(par = as.vector(coef(est_simple, 'phi')[-1,,drop=F]),
      K = K, W = data$W, weights_W = data$weights_W,
      norm_weights = data$norm_weights,
      group_E.prob = as.matrix(est_simple$posterior$posterior[,-1]), 
      ridge_penalty = make_TMatrix(K) * 1/est_simple$internal_parameters$control$prior_var_phi,
      gamma = 1, rank_F = data$rank_F,
      power_pi = 1, b_r = data$b_r,
      lambda = est_simple$parameters$eff_lambda)
  )
  
  grad_estimates <- cpp_gradient_phi(
    par = as.vector(coef(est_simple, 'phi')[-1,,drop=F]),
    K = K, W = data$W, weights_W = data$weights_W, norm_weights = data$norm_weights,
    group_E_prob = as.matrix(est_simple$posterior$posterior[,-1]), 
    ridge_penalty = make_TMatrix(K) * 1/est_simple$internal_parameters$control$prior_var_phi,
    gamma = 1, rank_F = data$rank_F,
    power_pi = 1, b_r = data$b_r,
    lambda = est_simple$parameters$eff_lambda
  )
  #Check that is mostly stationary
  expect_lte(max(abs(grad_estimates)), 1e-2)
  
  update_again <- update_moderator(phi = coef(est_simple, 'phi'), 
   group_E.prob = as.matrix(est_simple$posterior$posterior[,-1]),
   lambda = 1e-2 * nrow(dta), K = K, 
   gamma = 1, Fmatrix = data$Fmatrix, weights_W = data$weights_W,
   maxit_mod = 1000, extra_opt_args = list(reltol = 0, abstol = 0, method = 'BFGS'),
   rank_F = data$rank_F, adaptive_weight = data$adaptive_weight,
   ridge_phi = 1/est_simple$internal_parameters$control$prior_var_phi, 
   beta = est_simple$parameters$nullspace_beta, W = data$W)
  
  again_phi <- update_again$phi
  
  again_grad <- cpp_gradient_phi(
    par = as.vector(again_phi[-1,,drop=F]),
    K = K, W = data$W, norm_weights = data$norm_weights,
    weights_W = data$weights_W,
    group_E_prob = as.matrix(est_simple$posterior$posterior[,-1]), 
    ridge_penalty = make_TMatrix(K) * 1/est_simple$internal_parameters$control$prior_var_phi,
    gamma = 1, rank_F = data$rank_F,
    power_pi = 1, b_r = data$b_r,
    lambda = est_simple$parameters$eff_lambda
  )
  expect_lte(max(abs(again_grad)), 1e-5)
  
  expect_equivalent(est_simple$parameters$pi, 
    colMeans(est_simple$posterior$posterior_predictive[,-1]))
  
})


test_that("Gradient of Moderator Stage is Zero at Optimum (with weights)", {
  
  skip_on_cran()
  
  N <- 100
  K <- 3
  dta <- data.frame(
    state = sample(state.name[1:4], N, replace = T),
    letter = sample(letters[1:3], N, replace = T)
  )
  dta$mod <- runif(nrow(dta), -1, 1)
  dta$mod2 <- sample(LETTERS[1:5], N, replace = T)
  dta$y1 <- rbinom(nrow(dta), 1, plogis(runif(5, -1 , 1)[match(dta$state, state.name[1:5])] + runif(5, -1, 1)[match(dta$letter, letters)]))
  dta$y2 <- rbinom(nrow(dta), 1, plogis(runif(5, -1 , 1)[match(dta$state, state.name[1:5])] + runif(5, -1, 1)[match(dta$letter, letters)]))
  
  dta$y <- ifelse(rbinom(nrow(dta), 1, plogis(dta$mod)) == 1, dta$y1, dta$y2)
  dta$wgt <- abs(rcauchy(nrow(dta)))
  
  est_simple <- FactorHet(formula = y ~ state + letter, design = dta, K = K, 
    lambda = 1e-2, moderator = ~ mod + mod2, 
    weights = ~ wgt,
    control =  FactorHet_control(
      tolerance.parameters = 1e-8, 
      tolerance.logposterior = 0, 
      return_data = TRUE))  
  
  expect_false(all(est_simple$internal_parameters$data$weights_W == 1))
  
  expect_gte(min(diff(logLik(est_simple, 'log_posterior_seq'))), 
             - sqrt(.Machine$double.eps))
  
  data <- est_simple$internal_parameters$data
  data$norm_weights <- as.vector(data$weights_W/sum(data$weights_W))
  
  expect_equal(
    cpp_obj_phi(par = as.vector(coef(est_simple, 'phi')[-1,,drop=F]),
                K = K, W = data$W, norm_weights = data$norm_weights,
                weights_W = data$weights_W,
                group_E_prob = as.matrix(est_simple$posterior$posterior[,-1]), 
                ridge_penalty = make_TMatrix(K) * 1/est_simple$internal_parameters$control$prior_var_phi,
                gamma = 1, rank_F = data$rank_F,
                power_pi = 1, b_r = data$b_r,
                lambda = est_simple$parameters$eff_lambda),
    objective_moderator(par = as.vector(coef(est_simple, 'phi')[-1,,drop=F]),
                        K = K, W = data$W, weights_W = data$weights_W,
                        norm_weights = data$norm_weights,
                        group_E.prob = as.matrix(est_simple$posterior$posterior[,-1]), 
                        ridge_penalty = make_TMatrix(K) * 1/est_simple$internal_parameters$control$prior_var_phi,
                        gamma = 1, rank_F = data$rank_F,
                        power_pi = 1, b_r = data$b_r,
                        lambda = est_simple$parameters$eff_lambda)
  )
  
  grad_estimates <- cpp_gradient_phi(
    par = as.vector(coef(est_simple, 'phi')[-1,,drop=F]),
    K = K, W = data$W, weights_W = data$weights_W, norm_weights = data$norm_weights,
    group_E_prob = as.matrix(est_simple$posterior$posterior[,-1]), 
    ridge_penalty = make_TMatrix(K) * 1/est_simple$internal_parameters$control$prior_var_phi,
    gamma = 1, rank_F = data$rank_F,
    power_pi = 1, b_r = data$b_r,
    lambda = est_simple$parameters$eff_lambda
  )
  #Check that is mostly stationary
  expect_lte(max(abs(grad_estimates)), 1e-2)
  
  update_again <- update_moderator(phi = coef(est_simple, 'phi'), 
                                   group_E.prob = as.matrix(est_simple$posterior$posterior[,-1]),
                                   lambda = 1e-2 * nrow(dta), K = K, 
                                   gamma = 1, Fmatrix = data$Fmatrix, weights_W = data$weights_W,
                                   maxit_mod = 1000, extra_opt_args = list(reltol = 0, abstol = 0, method = 'BFGS'),
                                   rank_F = data$rank_F, adaptive_weight = data$adaptive_weight,
                                   ridge_phi = 1/est_simple$internal_parameters$control$prior_var_phi, 
                                   beta = est_simple$parameters$nullspace_beta, W = data$W)
  
  again_phi <- update_again$phi
  
  again_grad <- cpp_gradient_phi(
    par = as.vector(again_phi[-1,,drop=F]),
    K = K, W = data$W, norm_weights = data$norm_weights,
    weights_W = data$weights_W,
    group_E_prob = as.matrix(est_simple$posterior$posterior[,-1]), 
    ridge_penalty = make_TMatrix(K) * 1/est_simple$internal_parameters$control$prior_var_phi,
    gamma = 1, rank_F = data$rank_F,
    power_pi = 1, b_r = data$b_r,
    lambda = est_simple$parameters$eff_lambda
  )
  expect_lte(max(abs(again_grad)), 1e-5)
  
  norm_weight <- est_simple$internal_parameters$weights$weights_W/sum(est_simple$internal_parameters$weights$weights_W)
  
  expect_equivalent(est_simple$parameters$pi, 
    colSums(Diagonal(x = norm_weight) %*% as.matrix(est_simple$posterior$posterior_predictive[,-1]))
  )
  
})

