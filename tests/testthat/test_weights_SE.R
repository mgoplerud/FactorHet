context("Test weights & SEs")

test_that('Standard errors work with weights (K=1)', {
  
  N <- 200
  dta <- data.frame(
    state = sample(state.name[1:4], N, replace = T),
    letter = sample(letters[1:3], N, replace = T)
  )
  
  dta$group <- rep(sample(1:(N/4), (N/2), replace = T), each = 2)
  dta$task <- rep(1:(N/2), each = 2)
  dta$prof <- as.vector(sapply(1:(N/2), FUN=function(i){c('l', 'r')[sample(2)]}))
  dta$y <- as.vector(sapply(1:(N/2), FUN=function(i){sample(0:1)}))
  dta$mod <- runif((N/2), -1, 1)[dta$group]
  
  wgt <- abs(rcauchy(n = (N/2)))
  dta$weights <- wgt[dta$group]
  
  simple_logit <- FactorHet(formula = y ~ state * letter, design = dta,
                            K = 1, lambda = 0)
  weight_logit <- FactorHet(formula = y ~ state * letter, design = dta,
                            K = 1, lambda = 0, weights = ~ weights,
                            control = FactorHet_control(return_data = TRUE))
  
  X <- weight_logit$internal_parameters$data$X
  y <- weight_logit$internal_parameters$data$y  
  beta <- weight_logit$parameters$nullspace_beta
  weights <- weight_logit$internal_parameters$weights$weights
  
  p <- plogis(as.vector(X %*% beta))
  
  manual_weight <- solve(t(X) %*% Diagonal(x = p * (1-p) * weights) %*% X)
  fit_glm <- suppressWarnings(
    glm(y ~ 0 + as.matrix(X), family = binomial(), weights = weights)  
  )
  expect_equivalent(coef(fit_glm), as.vector(beta), tol = 1e-4)
  
  vcov_glm <- vcov(fit_glm)
  expect_equivalent(vcov_glm, as.matrix(manual_weight), tol = 1e-4)
  
  basis_M <- weight_logit$internal_parameters$data$basis_M
  vcov_analytic <- basis_M %*% manual_weight %*% t(basis_M)
  expect_equivalent(vcov(weight_logit), as.matrix(vcov_analytic))
  
  #Does it agree with Hessian of log-posterior with ridge?
  est_simple <- FactorHet(formula = y ~ state + letter, design = dta, K = 1, 
                          lambda = 0, moderator = ~ mod,
                          weights = ~ weights,
                          control = FactorHet_control(prior_var_beta = 1/3, 
                                                      log_method = 'standard',
                                                      return_data = TRUE))  
  pred_simple <- predict(est_simple)
  X <- est_simple$internal_parameters$data$X
  basis_M <- est_simple$internal_parameters$data$basis_M
  vcov_cons <- solve(t(X) %*% Diagonal(x = est_simple$internal_parameters$weights$weights * pred_simple * (1-pred_simple)) %*% X +
                       t(basis_M) %*% sparse_diag(c(0,rep(3, nrow(basis_M) - 1))) %*% (basis_M))
  vcov_analytic <- basis_M %*% vcov_cons %*% t(basis_M)
  expect_equivalent(vcov(est_simple), as.matrix(vcov_analytic))
})
