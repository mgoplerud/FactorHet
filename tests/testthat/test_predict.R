context('Test prediction functions')

test_that('Predict and OOS works for K = 1', {
  
  dta <- data.frame(
    state = sample(state.name[1:4], 1000, replace = T),
    letter = sample(letters[1:3], 1000, replace = T)
  )
  dta$y <- rbinom(nrow(dta), 1, plogis(runif(5, -1 , 1)[match(dta$state, state.name[1:5])] + runif(5, -1, 1)[match(dta$letter, letters)]))

  est_simple <- FactorHet(formula = y ~ state + letter, design = dta, K = 1, 
                          lambda = 1e-3, 
                          initialize = FactorHet_init(nrep =  1),
                          control =  FactorHet_control(iterations = 10, return_data = TRUE))
  coef_simple <- coef(est_simple)
  
  test_data <- data.frame(
    state = state.name[c(1,2,1)],
    letter = letters[1:3],
    mod = runif(3)
  )
  
  manual_pred <- mean(coef_simple[1,]) +
    coef_simple[match(paste0('state(', test_data$state, ')'), rownames(coef_simple)),] +
    coef_simple[match(paste0('letter(', test_data$letter, ')'), rownames(coef_simple)),]
  expect_equivalent(plogis(manual_pred), predict(est_simple, test_data))


})

test_that('Predict and OOS works for K = 2', {
  
  for (single_intercept in c(TRUE, FALSE)){
    dta <- data.frame(
      state = sample(state.name[1:4], 1000, replace = T),
      letter = sample(letters[1:3], 1000, replace = T)
    )
    dta$mod <- runif(nrow(dta), -1, 1)
    dta$mod2 <- sample(LETTERS[1:5], nrow(dta), replace = T)
    dta$y1 <- rbinom(nrow(dta), 1, plogis(runif(5, -1 , 1)[match(dta$state, state.name[1:5])] + runif(5, -1, 1)[match(dta$letter, letters)]))
    dta$y2 <- rbinom(nrow(dta), 1, plogis(runif(5, -1 , 1)[match(dta$state, state.name[1:5])] + runif(5, -1, 1)[match(dta$letter, letters)]))
    
    dta$y <- ifelse(rbinom(nrow(dta), 1, plogis(dta$mod)) == 1, dta$y1, dta$y2)
    
    est_simple <- FactorHet(formula = y ~ state * letter, design = dta, K = 2, 
                            lambda = 1e-3, moderator = ~ mod + mod2, 
                            initialize = FactorHet_init(nrep =  1),
                            control =  FactorHet_control(
                              single_intercept = single_intercept,
                              iterations = 10, return_data = TRUE))
    coef_simple <- coef(est_simple)
    
    test_data <- data.frame(
      state = state.name[c(1,2,1)],
      letter = letters[1:3],
      mod2 = LETTERS[1:3],
      mod = runif(3)
    )
    
    detailed_pred <- predict(est_simple, test_data, return = 'detailed')
    
    manual_pred <- coef_simple[match(paste0('state(', test_data$state, ')'), rownames(coef_simple)),] +
      coef_simple[match(paste0('letter(', test_data$letter, ')'), rownames(coef_simple)),] +
      coef_simple[match(paste0(test_data$state, '-', test_data$letter), rownames(coef_simple)),]
    
    manual_pred <- manual_pred + matrix(coef_simple[1,], byrow = TRUE, nrow = nrow(manual_pred), ncol = ncol(coef(est_simple)))
    
    manual_prob_yes <- apply(manual_pred, MARGIN = 2, plogis)
    
    expect_equivalent(manual_prob_yes, detailed_pred$prediction_by_cluster)
    
    mod_phi <- est_simple$parameters$phi[,'mod']
    mod_phi <- mod_phi[-1] * test_data$mod
    
    mod2_phi <- est_simple$parameters$phi[-1,]
    mod2_phi <- mod2_phi[match(paste0('mod2', test_data$mod2), names(mod2_phi))]
    mod2_phi[is.na(mod2_phi)] <- 0
    
    manual_postpred <- plogis(est_simple$parameters$phi[2,1] + mod_phi + mod2_phi)
    manual_postpred <- cbind(1 - manual_postpred, manual_postpred)
    expect_equivalent(manual_postpred, detailed_pred$posterior_predictive)
    
    expect_equivalent(rowSums(manual_postpred * manual_prob_yes), 
                      predict(est_simple, test_data))
    expect_equivalent(rowSums(manual_postpred * manual_prob_yes), detailed_pred$prediction)

  }

})

test_that('Predict and OOS works for K = 2 and groups', {
  
})
