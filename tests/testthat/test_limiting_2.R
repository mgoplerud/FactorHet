context('Test limiting cases (other)')

test_that('Non-flat prior and lambda = 0, check convergence', {
  
  dta <- data.frame(
    state = sample(state.name[1:4], 200, replace = T),
    letter = sample(letters[1:3], 200, replace = T)
  )
  dta$y <- rbinom(nrow(dta), 1, plogis(runif(5, -1 , 1)[match(dta$state, state.name[1:5])] + runif(5, -1, 1)[match(dta$letter, letters)]))
  dta$mod <- rnorm(nrow(dta))
  for (gamma in c(1,0)){
    for (single_intercept in c(TRUE, FALSE)){
      message(c(gamma, '-', single_intercept))
      
      est_simple <- FactorHet(formula = y ~ state + letter, 
        design = dta, K = 2, lambda = 0, moderator = ~ mod,
        control = FactorHet_control(prior_var_beta = 1, gamma = gamma,
        single_intercept = single_intercept)
      )
      expect_true(min(diff(logLik(est_simple, 'log_posterior_seq'))) > 0)
    }
  }
  
})
