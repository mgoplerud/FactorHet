
context("Model Based Optimization")


test_that(desc = "Check that MBO runs", {
  
  dta <- data.frame(
    state = sample(state.name[1:4], 100, replace = T),
    letter = sample(letters[1:3], 100, replace = T)
  )
  dta$mod <- runif(nrow(dta), -1, 1)
  dta$y1 <- rbinom(nrow(dta), 1, plogis(runif(5, -1 , 1)[match(dta$state, state.name[1:5])] + runif(5, -1, 1)[match(dta$letter, letters)]))
  dta$y2 <- rbinom(nrow(dta), 1, plogis(runif(5, -1 , 1)[match(dta$state, state.name[1:5])] + runif(5, -1, 1)[match(dta$letter, letters)]))
  
  dta$y <- ifelse(rbinom(nrow(dta), 1, plogis(dta$mod)) == 1, dta$y1, dta$y2)
  dta$wgt <- abs(rcauchy(nrow(dta)))
  

  m1 <- FactorHet_mbo(
    formula = y ~ state + letter, design = dta, K = 2, 
    moderator = ~ mod, weights = ~ wgt,
    mbo_control = FactorHet_mbo_control(iters = 2)
  )
  
  m2 <- FactorHet_mbo(
    formula = y ~ state + letter, design = dta, K = 2, 
    moderator = ~ mod, weights = ~ wgt,
    mbo_control = FactorHet_mbo_control(iters = 2, mbo_type = 'ridge')
  )
  
  m3 <- FactorHet_mbo(
    formula = y ~ state + letter, design = dta, K = 2, 
    moderator = ~ mod, weights = ~ wgt,
    mbo_control = FactorHet_mbo_control(mbo_design = 'random', iters = 2)
  )
  
  m4 <- FactorHet_mbo(
    formula = y ~ state + letter, design = dta, K = 2, 
    moderator = ~ mod, 
    mbo_control = FactorHet_mbo_control(mbo_design = data.frame(l = c(-3, -1)), iters = 0)
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
  
})
