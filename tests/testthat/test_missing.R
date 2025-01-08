context('Test that model runs on unusual data')

if (isTRUE(as.logical(Sys.getenv("CI")))){
  # If on CI
  NITER <- 2
  env_test <- "CI"
}else if (!identical(Sys.getenv("NOT_CRAN"), "true")){
  # If on CRAN
  NITER <- 2
  env_test <- "CRAN"
  set.seed(20)
}else{
  # If on local machine
  NITER <- 2000
  env_test <- 'local'
}

test_that('FactorHet runs with missing data (ungrouped)', {
  
  dta <- data.frame(
    state = sample(state.name[1:4], 1000, replace = T),
    letter = sample(letters[1:3], 1000, replace = T)
  )
  dta$mod <- runif(nrow(dta), -1, 1)
  dta$y1 <- rbinom(nrow(dta), 1, plogis(runif(5, -1 , 1)[match(dta$state, state.name[1:5])] + runif(5, -1, 1)[match(dta$letter, letters)]))
  dta$y2 <- rbinom(nrow(dta), 1, plogis(runif(5, -1 , 1)[match(dta$state, state.name[1:5])] + runif(5, -1, 1)[match(dta$letter, letters)]))
  
  dta$y <- ifelse(rbinom(nrow(dta), 1, plogis(dta$mod)) == 1, dta$y1, dta$y2)
  
  #Inject Missingness
  dta$y[c(3,56)] <- NA
  dta$mod[38] <- NA
  dta$state[28] <- NA
  dta$letter[1] <- NA
  dta[490,] <- NA
  
  est_simple <- FactorHet(formula = y ~ state + letter, design = dta,
   moderator = ~ mod, K = 2, lambda = 0, initialize = FactorHet_init(nrep = 1),
   control = FactorHet_control(return_data = TRUE, init_method = 'short_EM', 
   iterations = 10))
  #Ensure the estimated data matches what na.omit would provide
  expect_equal(nrow(est_simple$internal_parameters$data$X), nrow(na.omit(dta)))
  expect_equivalent(na.omit(dta)$mod, est_simple$internal_parameters$data$W[,2])
  expect_equivalent(na.omit(dta)$y, est_simple$internal_parameters$data$y)
  
})


test_that('FactorHet runs with white spaces in formula', {
  
  dta <- data.frame(
    `sta te` = sample(state.name[1:4], 1000, replace = T),
    `le-tter` = sample(letters[1:3], 1000, replace = T),
    normal = sample(LETTERS[2:5],  1000, replace = T), check.names = F
  )
  dta$`m od` <- runif(nrow(dta), -1, 1)
  dta$y1 <- rbinom(nrow(dta), 1, plogis(runif(5, -1 , 1)[match(dta[,1], state.name[1:5])] + runif(5, -1, 1)[match(dta[,2], letters)]))
  dta$y2 <- rbinom(nrow(dta), 1, plogis(runif(5, -1 , 1)[match(dta[,1], state.name[1:5])] + runif(5, -1, 1)[match(dta[,2], letters)]))
  
  dta$y <- ifelse(rbinom(nrow(dta), 1, plogis(dta$`m od`)) == 1, dta$y1, dta$y2)
  
  est_simple <- FactorHet(formula = y ~ `sta te` + `le-tter` + normal + . * ., design = dta,
     moderator = ~ `m od`, K = 2, lambda = 0, init = FactorHet_init(nrep = 1),
     control = FactorHet_control(return_data = TRUE, iterations = 10))
  
  n_levels <- lengths(apply(dta[,c("sta te", "le-tter", "normal")], MARGIN = 2, unique))
  check_params <- sum(apply(combn(3,2), MARGIN = 2, FUN=function(i){prod(n_levels[i])})) + sum(n_levels) + 1
  #Does the model have the right number of parameters
  expect_equal(nrow(coef(est_simple)), check_params)

})



test_that('FactorHet drops rare levels', {
  
  dta <- data.frame(
    state = sample(state.name[1:4], 1000, replace = T),
    letter = sample(letters[1:3], 1000, replace = T)
  )
  dta$mod <- runif(nrow(dta), -1, 1)
  dta$y1 <- rbinom(nrow(dta), 1, plogis(runif(5, -1 , 1)[match(dta$state, state.name[1:5])] + runif(5, -1, 1)[match(dta$letter, letters)]))
  dta$y2 <- rbinom(nrow(dta), 1, plogis(runif(5, -1 , 1)[match(dta$state, state.name[1:5])] + runif(5, -1, 1)[match(dta$letter, letters)]))
  
  dta$y <- ifelse(rbinom(nrow(dta), 1, plogis(dta$mod)) == 1, dta$y1, dta$y2)
  
  dta$letter <- ifelse(dta$state == 'Alabama' & dta$letter == 'a', "b", as.character(dta$letter))

  est_simple <- FactorHet(formula = y ~ state * letter, design = dta,
    moderator = ~ mod, K = 2, lambda = 1e-4, initialize = FactorHet_init(nrep = 1),
    control = FactorHet_control(return_data = TRUE, init_method = 'short_EM'))
  
  auto_predict <-  predict(est_simple, data.frame(state = 'Alabama', letter = 'a', mod = 1/2))
  
  disagg_pred <- colSums(coef(est_simple)[rownames(coef(est_simple)) %in% c('(Intercept)', 'state(Alabama)', 'letter(a)'),,drop=F])
  prob_2 <- plogis(sum(coef(est_simple, 'phi')[2,] * c(1, 1/2)))
  
  manual_predict <- sum(plogis(disagg_pred) * c(1-prob_2, prob_2))
  
  expect_equal(auto_predict, manual_predict)
  
})
