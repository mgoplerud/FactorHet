
context("Model Based Optimization")

if (isTRUE(as.logical(Sys.getenv("CI")))){
  # If on CI
  NITER <- 2
  env_test <- "CI"
}else if (!identical(Sys.getenv("NOT_CRAN"), "true")){
  # If on CRAN
  NITER <- 2
  env_test <- "CRAN"
  set.seed(13)
}else{
  # If on local machine
  NITER <- 2000
  env_test <- 'local'
}

test_that(desc = "Check that MBO runs", {
  
  dta <- data.frame(
    state = sample(state.name[1:4], 100, replace = T),
    letter = sample(letters[1:3], 100, replace = T)
  )
  dta$mod <- runif(nrow(dta), -1, 1)
  dta$y1 <- rbinom(nrow(dta), 1, 
      plogis(runif(5, -1 , 1)[match(dta$state, state.name[1:5])] + runif(5, -1, 1)[match(dta$letter, letters)]))
  dta$y2 <- rbinom(nrow(dta), 1, 
      plogis(runif(5, -1 , 1)[match(dta$state, state.name[1:5])] + runif(5, -1, 1)[match(dta$letter, letters)]))
  
  dta$y <- ifelse(rbinom(nrow(dta), 1, plogis(dta$mod)) == 1, dta$y1, dta$y2)
  dta$wgt <- abs(rcauchy(nrow(dta)))

  m1 <- FactorHet_mbo(
    formula = y ~ state + letter, design = dta, K = 2, 
    moderator = ~ mod, weights = ~ wgt,
    control = FactorHet_control(iterations = NITER, tolerance.logposterior = 0.1),
    mbo_control = FactorHet_mbo_control(iters = 2, mbo_method = 'regr.lm')
  )
  
  m2 <- FactorHet_mbo(
    formula = y ~ state + letter, design = dta, K = 2, 
    moderator = ~ mod, weights = ~ wgt,
    control = FactorHet_control(iterations = NITER, tolerance.logposterior = 0.1),
    mbo_control = FactorHet_mbo_control(iters = 1, mbo_method = 'regr.lm',
                                        mbo_type = 'ridge')
  )
  
  m3 <- FactorHet_mbo(
    formula = y ~ state + letter, design = dta, K = 2, 
    moderator = ~ mod, weights = ~ wgt,
    control = FactorHet_control(iterations = NITER, tolerance.logposterior = 0.1),
    mbo_control = FactorHet_mbo_control(mbo_design = 'random', 
                                        mbo_method = 'regr.lm',
                                        iters = 1)
  )
  
  m4 <- FactorHet_mbo(
    formula = y ~ state + letter, design = dta, K = 1, 
    moderator = ~ mod, 
    control = FactorHet_control(iterations = NITER, tolerance.logposterior = 0.1),
    mbo_control = FactorHet_mbo_control(
      mbo_method = 'regr.lm', mbo_noisy = FALSE,
      mbo_design = data.frame(l = c(-3, -1)), iters = 0)
  )
  
  expect_false(all(m1$internal_parameters$weights$weights_W == 1))
  expect_false(all(m2$internal_parameters$weights$weights_W == 1))
  expect_true(all(m4$internal_parameters$weights$weights_W == 1))
  expect_identical(m2$parameters$eff_lambda, 0)
  
  # Confirm that they ran
  expect_s3_class(m1, 'FactorHet')
  expect_s3_class(m2, 'FactorHet')
  expect_s3_class(m3, 'FactorHet')
  expect_s3_class(m4, 'FactorHet')
  
  vis_mbo <- visualize_MBO(m1)
  expect_s3_class(vis_mbo, 'ggplot')
})
