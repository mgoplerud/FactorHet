context('Test limiting cases (lambda = 0, lambda = Inf)')

test_that('Agrees with GLM (lambda = 0)', {
  
  dta <- data.frame(
    state = sample(state.name[1:2], 100, replace = T),
    letter = sample(letters[1:2], 100, replace = T)
  )
  dta$y <- rbinom(nrow(dta), 1, plogis(runif(5, -1 , 1)[match(dta$state, state.name[1:5])] + runif(5, -1, 1)[match(dta$letter, letters)]))
  
  est_simple <- FactorHet(y ~ state + letter + . * ., design = dta, 
    K = 1, lambda = 0, 
    control = FactorHet_control(return_data = TRUE,
      tolerance.logposterior = 0, 
      tolerance.parameters = 0, iterations = 50))
  
  est_glm <- glm(est_simple$internal_parameters$data$y ~ 0 + as.matrix(est_simple$internal_parameters$data$X), family = binomial)
  
  #Check same betas
  expect_equivalent(as.vector(est_simple$parameters$nullspace_beta), coef(est_glm), 
    tolerance = 1e-4, scale = 1)
  #Check same log-lik
  expect_equivalent(logLik(est_simple), as.numeric(logLik(est_glm)), tolerance = 1e-4, scale = 1)
  #Intercept is last column in nullspace  
  expect_equivalent(est_simple$parameters$beta[1], 
    est_simple$parameters$nullspace_beta[nrow(est_simple$parameters$nullspace_beta)],
    tolerance = 1e-4, scale = 1)
})

test_that('Agrees with GLM (lambda = Inf)', {
  
  dta <- data.frame(
    state = sample(state.name[1:2], 100, replace = T),
    letter = sample(letters[1:2], 100, replace = T)
  )
  dta$y <- rbinom(nrow(dta), 1, plogis(runif(5, -1 , 1)[match(dta$state, state.name[1:5])] + runif(5, -1, 1)[match(dta$letter, letters)]))
  
  est_simple <- FactorHet(formula = y ~ state + letter,  design = dta, 
    K = 1, lambda = 10^5, 
    control = FactorHet_control(return_data = TRUE,
      calc_df = FALSE,
      tolerance.logposterior = 0,
      tolerance.parameters = 0, iterations = 50))

  #Intercept is last column in nullspace  
  expect_equivalent(est_simple$parameters$beta[1], 
    est_simple$parameters$nullspace_beta[nrow(est_simple$parameters$nullspace_beta)])
  #Only non-zero parameter is intercept
  expect_equivalent(est_simple$parameters$beta[1], qlogis(mean(dta$y)), tolerance = 1e-6)
  expect_equal(est_simple$parameters$beta[-1], rep(0, nrow(est_simple$parameters$beta) - 1))
  
})

test_that('FactorHet runs with K = 1 without error',{
  
  dta <- data.frame(
    state = sample(state.name[1:4], 50, replace = T),
    letter = sample(letters[1:3], 50, replace = T)
  )
  dta$y <- rbinom(nrow(dta), 1, plogis(runif(5, -1 , 1)[match(dta$state, state.name[1:5])] + runif(5, -1, 1)[match(dta$letter, letters)]))
  
  
  test_repeat_one <- tryCatch(FactorHet(formula = y ~ state * letter, 
                                        design = dta, K = 1, lambda = 0, 
                                        control = FactorHet_control(iterations = 1)), error = function(e){NULL})
  
  expect_false(is.null(test_repeat_one))
})


test_that('FactorHet predicts with excluded interactions',{
  
  dta <- data.frame(
    state = sample(state.name[1:4], 50, replace = T),
    letter = sample(letters[1:3], 50, replace = T),
    stringsAsFactors = F
  )
  dta$y <- rbinom(nrow(dta), 1, plogis(runif(5, -1 , 1)[match(dta$state, state.name[1:5])] + runif(5, -1, 1)[match(dta$letter, letters)]))
  full_dta <- dta
  dta <- subset(dta, !(state == state.name[1] & letter == 'a'))
  
  test_repeat_one <- tryCatch(FactorHet(formula = y ~ state * letter, 
                                        K = 1, lambda = 0, design = dta, 
                                        control = FactorHet_control(tolerance.logposterior = 0,
                                                                    tolerance.parameters = 0, iterations = 1)), 
                              error = function(e){NULL})
  
  coef_simple <- coef(test_repeat_one)
  expect_false(is.null(test_repeat_one))
  
  pred_dta <- suppressWarnings(predict(test_repeat_one, dta))
  pred_full_dta <- suppressWarnings(predict(test_repeat_one, full_dta))
  
  man_state <- coef_simple[match(paste0('state(', full_dta$state, ')'), rownames(coef_simple)),]
  man_state[is.na(man_state)] <- 0
  man_letter <- coef_simple[match(paste0('letter(', full_dta$letter, ')'), rownames(coef_simple)),]
  man_letter[is.na(man_letter)] <- 0
  manual_pred <- coef_simple[1,1] + man_state + man_letter
  manual_inter <- coef_simple[match(paste0(full_dta$state, '-', full_dta$letter), rownames(coef_simple)),]
  manual_inter[is.na(manual_inter)] <- 0
  #Does prediction still match up
  expect_equivalent(pred_full_dta, plogis(manual_pred + manual_inter))
  # Predict runs with one observation
  expect_vector(predict(test_repeat_one, full_dta[sample(1:nrow(full_dta), 1) ,]), numeric())
})

