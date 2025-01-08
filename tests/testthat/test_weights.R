context("Check Weights")

if (isTRUE(as.logical(Sys.getenv("CI")))){
  # If on CI
  NITER <- 2
  env_test <- "CI"
}else if (!identical(Sys.getenv("NOT_CRAN"), "true")){
  # If on CRAN
  NITER <- 2
  env_test <- "CRAN"
  set.seed(15)
}else{
  # If on local machine
  NITER <- 2000
  env_test <- 'local'
}

test_that('Weights for K = 2, conjoint', {
  
  dta <- data.frame(
    state = sample(state.name[1:4], 1000, replace = T),
    letter = sample(letters[1:3], 1000, replace = T)
  )

  dta$group <- rep(sample(1:250, 500, replace = T), each = 2)
  dta$task <- rep(1:500, each = 2)
  dta$prof <- as.vector(sapply(1:500, FUN=function(i){c('l', 'r')[sample(2)]}))
  dta$y <- as.vector(sapply(1:500, FUN=function(i){sample(0:1)}))
  dta$mod <- runif(500, -1, 1)[dta$group]
  
  weights <- abs(rcauchy(nrow(dta)))
  weights[median(weights) * 500 < weights] <- median(weights) * 500
  
  est_simple <- expect_error(
    FactorHet(formula = y ~ state + letter, design = dta, K = 2, 
                          lambda = 1e-2, moderator = ~ mod, 
                          weights = ~ weights,
                          group = ~ group, task = ~ task, choice_order = ~ prof),
    message = 'weights')
  
  dta$weights <- as.vector(sapply(1:500, FUN=function(i){rep(abs(rcauchy(1)), 2)}))
  
  est_simple <- expect_error(FactorHet(formula = y ~ state + letter, design = dta, K = 2, 
    lambda = 1e-2, moderator = ~ mod, 
    weights = ~ weights,
    group = ~ group, task = ~ task, choice_order = ~ prof),
    message = 'weights')
  
  wgt <- abs(rcauchy(n = 500))
  dta$outcome_weights <- wgt[dta$group]
  
  est_simple <- FactorHet(formula = y ~ state + letter, design = dta, K = 2, 
      lambda = 1e-3, moderator = ~ mod, 
      weights = ~ outcome_weights,
      group = ~ group, task = ~ task, 
      choice_order = ~ prof,
      control = FactorHet_control(
        init_method = 'random_member'),
      init = FactorHet_init(return_all = TRUE, short_EM = FALSE))
  
  expect_false(all(est_simple[[1]]$internal_parameters$weights$weights_W == 1))
  lp_simple <- lapply(est_simple, FUN=function(i){logLik(i, 'log_posterior_seq')})
  expect_true(min(sapply(lp_simple, FUN=function(i){range(diff(i))})) > -sqrt(.Machine$double.eps))

  est_simple_clip <- FactorHet(formula = y ~ state + letter, design = dta, K = 2, 
    lambda = 1e-3, moderator = ~ mod, 
    weights = ~ outcome_weights,
    group = ~ group, task = ~ task, 
    choice_order = ~ prof,
    control = FactorHet_control(tau_method = 'clip', init_method = 'random_member'),
    init = FactorHet_init(return_all = TRUE, nrep =2, short_EM = FALSE))
  
  lp_simple_clip <- lapply(est_simple_clip, FUN=function(i){logLik(i, 'log_posterior_seq')})
  expect_true(min(sapply(lp_simple_clip, FUN=function(i){range(diff(i))})) > -sqrt(.Machine$double.eps))
  
  est_simple_cpp <- suppressMessages(FactorHet(formula = y ~ state + letter, design = dta, K = 2, 
       lambda = 1e-3, moderator = ~ mod, 
       weights = ~ outcome_weights,
       group = ~ group, task = ~ task, 
       choice_order = ~ prof,
       control = FactorHet_control(beta_method = 'cg', init_method = 'random_member'),
       init = FactorHet_init(return_all = TRUE, nrep = 2, short_EM = FALSE)))
  
  lp_simple_cpp <- lapply(est_simple_cpp, FUN=function(i){logLik(i, 'log_posterior_seq')})
  expect_true(min(sapply(lp_simple_cpp, FUN=function(i){range(diff(i))})) > -sqrt(.Machine$double.eps))
  
})

test_that('Weights for K = 2, free intercept', {
  
  dta <- data.frame(
    state = sample(state.name[1:4], 1000, replace = T),
    letter = sample(letters[1:3], 1000, replace = T)
  )
  
  dta$group <- rep(sample(1:250, 500, replace = T), each = 2)
  dta$task <- rep(1:500, each = 2)
  dta$prof <- as.vector(sapply(1:500, FUN=function(i){c('l', 'r')[sample(2)]}))
  dta$y <- as.vector(sapply(1:500, FUN=function(i){sample(0:1)}))
  dta$mod <- runif(500, -1, 1)[dta$group]

  wgt <- abs(rcauchy(n = 500))
  dta$weights <- wgt[dta$group]
  
  est_simple <- FactorHet(formula = y ~ state + letter, design = dta, K = 2, 
      lambda = 1e-3, moderator = ~ mod, 
      weights = ~ weights,
      group = ~ group, task = ~ task, 
      choice_order = ~ prof,
      control = FactorHet_control(single_intercept = FALSE,
        init_method = 'random_member'),
      init = FactorHet_init(return_all = TRUE, short_EM = FALSE))
  
  lp_simple <- lapply(est_simple, FUN=function(i){logLik(i, 'log_posterior_seq')})
  expect_true(min(sapply(lp_simple, FUN=function(i){range(diff(i))})) > -sqrt(.Machine$double.eps))
  
  est_simple_clip <- FactorHet(formula = y ~ state + letter, design = dta, K = 2, 
       lambda = 1e-3, moderator = ~ mod, 
       weights = ~ weights,
       group = ~ group, task = ~ task, 
       choice_order = ~ prof,
       control = FactorHet_control(
         single_intercept = FALSE,
         tau_method = 'clip', init_method = 'random_member'),
       init = FactorHet_init(return_all = TRUE, short_EM = FALSE))
  
  lp_simple_clip <- lapply(est_simple_clip, FUN=function(i){logLik(i, 'log_posterior_seq')})
  expect_true(min(sapply(lp_simple_clip, FUN=function(i){range(diff(i))})) > -sqrt(.Machine$double.eps))

  est_simple_cpp <- FactorHet(formula = y ~ state + letter, 
      design = dta, K = 2, 
      lambda = 1e-3, moderator = ~ mod, 
      weights = ~ weights,
      group = ~ group, task = ~ task, 
      choice_order = ~ prof,
      control = FactorHet_control(
        single_intercept = FALSE, 
        beta_method = 'cpp', init_method = 'random_member'),
      init = FactorHet_init(return_all = TRUE, short_EM = FALSE))
  
  lp_simple_cpp <- lapply(est_simple_cpp, FUN=function(i){logLik(i, 'log_posterior_seq')})
  expect_true(min(sapply(lp_simple_cpp, FUN=function(i){range(diff(i))})) > -sqrt(.Machine$double.eps))

})

test_that('Post-Regression Functions Work with Weights', {
  
  dta <- data.frame(
    state = sample(state.name[1:4], 1000, replace = T),
    letter = sample(letters[1:3], 1000, replace = T)
  )
  
  dta$group <- rep(sample(1:250, 500, replace = T), each = 2)
  dta$task <- rep(1:500, each = 2)
  dta$prof <- as.vector(sapply(1:500, FUN=function(i){c('l', 'r')[sample(2)]}))
  dta$y <- as.vector(sapply(1:500, FUN=function(i){sample(0:1)}))
  dta$mod <- runif(500, -1, 1)[dta$group]
  
  wgt <- abs(rcauchy(n = 500))
  dta$weights <- wgt[dta$group]
  
  est_simple <- FactorHet(formula = y ~ state + letter, design = dta, K = 3, 
    lambda = 1e-3, moderator = ~ mod, 
    weights = ~ weights,
    group = ~ group, task = ~ task, 
    choice_order = ~ prof)

  expect_false(all(est_simple$internal_parameters$weights$weights_W == 1))
  pred_test <- predict(est_simple, newdata = dta)
  pred_AME <- AME(est_simple)
  pred_mfx <- margeff_moderators(est_simple)
})

test_that('Estimation methods for beta work with weights', {
  
  dta <- data.frame(
    state = sample(state.name[1:4], 1000, replace = T),
    letter = sample(letters[1:3], 1000, replace = T)
  )
  
  dta$group <- rep(sample(1:250, 500, replace = T), each = 2)
  dta$task <- rep(1:500, each = 2)
  dta$prof <- as.vector(sapply(1:500, FUN=function(i){c('l', 'r')[sample(2)]}))
  dta$y <- as.vector(sapply(1:500, FUN=function(i){sample(0:1)}))
  dta$mod <- runif(500, -1, 1)[dta$group]
  
  wgt <- abs(rcauchy(n = 500))
  dta$weights <- wgt[dta$group]
  
  est_simple <- FactorHet(formula = y ~ state + letter, design = dta, K = 3, 
    lambda = 1e-3, moderator = ~ mod, 
    weights = ~ weights,
    group = ~ group, task = ~ task, 
    choice_order = ~ prof,
    control = FactorHet_control(return_data = TRUE, init_method = 'random_member'),
    init = FactorHet_init(nrep = 1, short_EM = FALSE))
  
  weights <- est_simple$internal_parameters$data$weights
  y <- est_simple$internal_parameters$data$y
  X <- est_simple$internal_parameters$data$X
  W <- est_simple$internal_parameters$data$W
  group_E.prob <- est_simple$internal_parameters$data$group_E.prob
  beta <- est_simple$parameters$nullspace_beta
  phi <- est_simple$parameters$phi
  K <- ncol(beta)
  group_mapping <- sparseMatrix(i = 1:nrow(X), j = match(est_simple$internal_parameters$group$group_mapping, est_simple$internal_parameters$group$unique_group), x = 1)
  obs.E.prob <- apply(group_E.prob, MARGIN = 2, FUN=function(i){as.vector(group_mapping %*% i)})
  xb <- as.matrix(X %*% beta)
  E.omega <- obs.E.prob/(2 * xb) * tanh(xb/2)
  list_Eridge <- est_simple$internal_parameters$data$E_ridge
  
  cpp_separate <- cpp_beta(K = K, X = X, E_ridge = list_Eridge, y = y, weights = weights, 
    E_omega = E.omega, obs_E_prob = obs.E.prob)
  
  manual_separate <- sapply(1:K, FUN=function(k){
    as.vector(solve(t(X) %*% Diagonal(x = E.omega[,k]) %*% X + list_Eridge[[k]],
          t(X) %*% (weights * obs.E.prob[,k] * (y-1/2))))
  })
  
  expect_equivalent(cpp_separate, manual_separate)
  
  cg_separate <- cg_custom(K = K, X = X, 
    list_ridge = list_Eridge, 
    omega = as.matrix(E.omega), 
    s = weights * (y - 1/2), old_beta = matrix(0, ncol = K, nrow = ncol(X)), 
    weights = as.matrix(obs.E.prob), 
    tol = sqrt(.Machine$double.eps), it_max = 0)

  expect_equivalent(cg_separate$beta, manual_separate, tol = 1e-5)

  blocked_X <- cbind(1, kronecker(diag(K), X[,-ncol(X)]))
  blocked_E <- bdiag(c(0, list_Eridge))
  blocked_E <- blocked_E[-(1 + ncol(X) * 1:K),-(1 + ncol(X) * 1:K)]
  
  cpp_blocked <- cpp_beta_plain(X = blocked_X,
         s = rep(weights * (y - 1/2), K) * as.vector(obs.E.prob),
         K = K, omega = sparse_diag(as.vector(E.omega)),
         ridge = blocked_E
  )
  
  cg_b <- cg_custom(K = 1, X = blocked_X,
      s = rep((y-1/2) * weights, K) * as.vector(obs.E.prob), 
      list_ridge = list(blocked_E), tol = sqrt(.Machine$double.eps),
      it_max = 0, old_beta = matrix(rep(0, ncol(blocked_X))), 
      weights = matrix(0, nrow = 1, ncol = 0),
      omega = matrix(as.vector(E.omega))
  )
  
  blocked_beta <- solve(t(blocked_X) %*% Diagonal(x = as.vector(E.omega)) %*% blocked_X + blocked_E, 
        t(blocked_X) %*% as.vector(rep(y - 1/2, K) * weights * as.vector(obs.E.prob)))
  
  expect_equivalent(cg_b$beta, as.vector(cpp_blocked), tol = 1e-4, scale = 1)
  expect_equivalent(as.vector(blocked_beta), as.vector(cpp_blocked))

  prior_beta <- beta
  prior_beta[,] <- 0
  # Test for separate intercepts
  list_loop <- list()  
  for (mth in c('base', 'cpp', 'cg')){
    list_loop[[mth]] <- update_beta(X=X,y=y, E.omega = E.omega, 
      weights = weights, p_X = ncol(X),
      obs.E.prob = obs.E.prob, cg_it = 0,
      E.ridge = list_Eridge, K= 3, method = mth, 
      global_int = F, prior_beta = prior_beta)
  }
  expect_equivalent(list_loop$base, list_loop$cg, tol = 1e-5)
  expect_equivalent(list_loop$base, list_loop$cpp)
  
  # Test for single intercept
  list_loop <- list()  
  p_X <- ncol(X)
  for (mth in c('base', 'cpp', 'cg')){
    
    list_loop[[mth]] <- update_beta(X=NULL, blocked_X = blocked_X, 
      p_X = p_X, prior_beta = prior_beta,
      weights = weights,
      y=y,E.omega = E.omega, obs.E.prob = obs.E.prob, 
      method = mth,
      E.ridge = list_Eridge, K=K, 
      global_int = TRUE, cg_it = 0)
    
  }
  expect_equivalent(list_loop$base, list_loop$cg, tol = 1e-5)
  expect_equivalent(list_loop$base, list_loop$cpp)
  
  simple_cpp <- simple_logit(y = y, X = X, iterations = 15, weights = weights,
    obs.E.prob = obs.E.prob, beta_method = 'cpp')
  simple_cg <- simple_logit(y = y, X = X, iterations = 15, weights = weights,
    obs.E.prob = obs.E.prob, beta_method = 'cg', beta_cg_it = 0)
  expect_equivalent(simple_cpp, simple_cg, tol = 1e-5)
  
})

test_that('Weights work for gamma = 0', {
  
    dta <- data.frame(
      state = sample(state.name[1:4], 1000, replace = T),
      letter = sample(letters[1:3], 1000, replace = T)
    )
    
    dta$group <- rep(sample(1:250, 500, replace = T), each = 2)
    dta$task <- rep(1:500, each = 2)
    dta$prof <- as.vector(sapply(1:500, FUN=function(i){c('l', 'r')[sample(2)]}))
    dta$y <- as.vector(sapply(1:500, FUN=function(i){sample(0:1)}))
    dta$mod <- runif(500, -1, 1)[dta$group]
    
    wgt <- abs(rcauchy(n = 500))
    dta$weights <- wgt[dta$group]
    
    est_simple <- FactorHet(formula = y ~ state + letter, design = dta, K = 2, 
        lambda = 1e-3, moderator = ~ mod, 
        weights = ~ weights,
        group = ~ group, task = ~ task, 
        choice_order = ~ prof,
        control = FactorHet_control(gamma = 0, 
                                    init_method = 'random_member'),
        init = FactorHet_init(return_all = TRUE, short_EM = FALSE))
    
    expect_false(all(est_simple[[1]]$internal_parameters$weights$weights_W == 1))
    
    lp_simple <- lapply(est_simple, FUN=function(i){logLik(i, 'log_posterior_seq')})
    expect_true(min(sapply(lp_simple, FUN=function(i){range(diff(i))})) > -sqrt(.Machine$double.eps))
    
})

test_that('Weights work for short EM ', {
  
  dta <- data.frame(
    state = sample(state.name[1:4], 1000, replace = T),
    letter = sample(letters[1:3], 1000, replace = T)
  )
  
  dta$group <- rep(sample(1:250, 500, replace = T), each = 2)
  dta$task <- rep(1:500, each = 2)
  dta$prof <- as.vector(sapply(1:500, FUN=function(i){c('l', 'r')[sample(2)]}))
  dta$y <- as.vector(sapply(1:500, FUN=function(i){sample(0:1)}))
  dta$mod <- runif(500, -1, 1)[dta$group]
  
  wgt <- abs(rcauchy(n = 500))
  dta$wgt_interal <- wgt[dta$group]
  
  est_simple <- FactorHet(formula = y ~ state + letter, design = dta, K = 2, 
                          lambda = 1e-3, moderator = ~ mod, 
                          weights = ~ wgt_interal,
                          group = ~ group, task = ~ task, 
                          choice_order = ~ prof)
  
  expect_false(all(est_simple$internal_parameters$weights$weights_W == 1))
  
  lp_simple <- logLik(est_simple, 'log_posterior_seq')
  expect_true(min(diff(lp_simple)) > -sqrt(.Machine$double.eps))
  
})

test_that('Weights work for inital values and MBO', {
  
  dta <- data.frame(
    state = sample(state.name[1:4], 1000, replace = T),
    letter = sample(letters[1:3], 1000, replace = T)
  )
  
  dta$group <- rep(sample(1:250, 500, replace = T), each = 2)
  dta$task <- rep(1:500, each = 2)
  dta$prof <- as.vector(sapply(1:500, FUN=function(i){c('l', 'r')[sample(2)]}))
  dta$y <- as.vector(sapply(1:500, FUN=function(i){sample(0:1)}))
  dta$mod <- runif(500, -1, 1)[dta$group]
  
  wgt <- abs(rcauchy(n = 500))
  dta$wgt_interal <- wgt[dta$group]
  
  fit_MBO <- FactorHet_mbo(formula = y ~ state + letter, 
     design = dta, K = 2, 
     moderator = ~ mod, 
     weights = ~ wgt_interal,
     group = ~ group, task = ~ task, 
     choice_order = ~ prof,
     mbo_control = FactorHet_mbo_control(mbo_design = data.frame(l = c(-3, -2)), 
        iters = 0))

  expect_false(all(fit_MBO$internal_parameters$weights$weights_W == 1))
  
  lp_simple <- logLik(fit_MBO, 'log_posterior_seq')
  expect_true(min(diff(lp_simple)) > -1e-6)
  
})



