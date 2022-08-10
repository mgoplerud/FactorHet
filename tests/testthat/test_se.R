
test_that('K = 1 and lambda = 0 agrees with ridge ', {
  require(Matrix)
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