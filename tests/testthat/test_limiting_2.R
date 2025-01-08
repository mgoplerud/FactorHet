context('Test limiting cases (other)')

if (isTRUE(as.logical(Sys.getenv("CI")))){
  # If on CI
  NITER <- 2
  env_test <- "CI"
}else if (!identical(Sys.getenv("NOT_CRAN"), "true")){
  # If on CRAN
  NITER <- 2
  env_test <- "CRAN"
  set.seed(11)
}else{
  # If on local machine
  NITER <- 2000
  env_test <- 'local'
}

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
        initialize = FactorHet_init(nrep = ifelse(NITER > 2, 5, 1), short_EM_it = NITER),
        control = FactorHet_control(prior_var_beta = 1, gamma = gamma,
        single_intercept = single_intercept, iterations = NITER)
      )
      expect_true(min(diff(logLik(est_simple, 'log_posterior_seq'))) > 0)
    }
  }
  
})
