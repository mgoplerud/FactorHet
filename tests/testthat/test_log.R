context('test LOG')

test_that('LOG and Standard Agree with No Interactions', {
  
  dta <- data.frame(
    state = sample(state.name[1:4], 1000, replace = T),
    letter = sample(letters[1:3], 1000, replace = T),
    region = sample(state.region, 1000, replace = T)
  )
  dta$group <- rep(sample(1:250, 500, replace = T), each = 2)
  dta$task <- rep(1:500, each = 2)
  dta$prof <- as.vector(sapply(1:500, FUN=function(i){c('l', 'r')[sample(2)]}))
  dta$y <- as.vector(sapply(1:500, FUN=function(i){sample(0:1)}))
  dta$mod <- runif(500, -1, 1)[dta$group]
  
  est_simple_factorial <- FactorHet(formula = y ~ state + letter + region, design = dta, K = 1, 
      lambda = 1e-2, moderator = ~ mod, 
      control = FactorHet_control(iterations = 400, tolerance.parameters = 0,
        tolerance.logposterior = 0, log_method = 'log_ginv',
        return_data = TRUE),
      group = ~ group, task = ~ task, choice = ~ prof)

  est_simple_factorial_std <- FactorHet(formula = y ~ state + letter + region, design = dta, K = 1, 
      lambda = 1e-2, moderator = ~ mod, 
      control = FactorHet_control(iterations = 400, tolerance.parameters = 0, 
        tolerance.logposterior = 0,
        log_method = 'standard', return_data = TRUE),
      group = ~ group, task = ~ task, choice = ~ prof)
  # Check the weights are the same
  expect_equal(est_simple_factorial_std$internal_parameters$adaptive_weight,
  est_simple_factorial$internal_parameters$adaptive_weight)
  # Check the estimates are the same
  expect_equal(coef(est_simple_factorial), coef(est_simple_factorial_std))
  # Check the SIZE of the nullspace is the same, i.e. no extra coefficients
  expect_equal(length(est_simple_factorial$parameters$nullspace_beta),
               length(est_simple_factorial_std$parameters$nullspace_beta))
  expect_equal(length(est_simple_factorial$internal_parameters$data$Fmatrix), 3)
})

test_that('Test Projection Method', {
  
  # Projection methods all agree with no B&R weights and no interactions
  message('ADD TEST HERE for B&R Weights')

})

test_that('B&R weights align with analytical', {
  
  
  dta <- data.frame(
    state = sample(state.name[1:4], 1000, replace = T),
    letter = sample(letters[1:3], 1000, replace = T)
  )
  dta$group <- rep(sample(1:250, 500, replace = T), each = 2)
  dta$task <- rep(1:500, each = 2)
  dta$prof <- as.vector(sapply(1:500, FUN=function(i){c('l', 'r')[sample(2)]}))
  dta$y <- as.vector(sapply(1:500, FUN=function(i){sample(0:1)}))
  dta$mod <- runif(500, -1, 1)[dta$group]
  dta$fake_task <- paste0(dta$task, dta$prof)
  est_simple_factorial <- FactorHet(formula = y ~ state + letter, design = dta, K = 2, 
    lambda = 1e-2, moderator = ~ mod, 
    control = FactorHet_control(iterations = 1, return_data = TRUE),
    initialize = FactorHet_init(nrep = 1),
    group = ~ group, task = ~ fake_task)
  
  # Check the weights for a factorial design
  # conform with Bondell and Reich (2009)
  levels_dta_state <- table(dta$state)
  br_state <- 1/(length(levels_dta_state) + 1) * outer(levels_dta_state, levels_dta_state, FUN=function(x,y){sqrt(x + y)}) *
    1/sqrt(nrow(dta))
  br_state <- br_state[upper.tri(br_state)]
  expect_equal(
    est_simple_factorial$internal_parameters$adaptive_weight[[1]]$state,
    br_state
  )
  levels_dta_letter <- table(dta$letter)
  br_letter <- 1/(length(levels_dta_letter) + 1) * outer(levels_dta_letter, levels_dta_letter, FUN=function(x,y){sqrt(x + y)}) *
    1/sqrt(nrow(dta))
  br_letter <- br_letter[upper.tri(br_letter)]
  expect_equal(
    est_simple_factorial$internal_parameters$adaptive_weight[[1]]$letter,
    br_letter
  )
  
  # Check that even in a model with interactions,
  # the weights conform for the non-interacted factors
  dta$fake_one <- sample(LETTERS[1:5], 1000, replace = T)
  dta$fake_two <- sample(LETTERS[1:5], 1000, replace = T)
  est_fact_2 <- FactorHet(formula = y ~ state + letter + fake_one * fake_two, design = dta, K = 2, 
      lambda = 1e-2, moderator = ~ mod, 
      init = FactorHet_init(nrep = 1),
      control = FactorHet_control(iterations = 1, return_data = TRUE),
      group = ~ group, task = ~ fake_task)
  
  expect_equal(est_fact_2$internal_parameters$adaptive_weight[[1]]$state, br_state)
  expect_equal(est_fact_2$internal_parameters$adaptive_weight[[1]]$letter, br_letter)
  # Check they DO NOT align with the interacted one
  levels_dta_fake <- table(dta$fake_one)
  br_fake <- 1/(length(levels_dta_fake) + 1) * outer(levels_dta_fake, levels_dta_fake, FUN=function(x,y){sqrt(x + y)}) *
    1/sqrt(nrow(dta))
  br_fake <- br_fake[upper.tri(br_fake)]
  expect_true(
    !isTRUE(all.equal(est_fact_2$internal_parameters$adaptive_weight[[1]]$fake_one,
    br_fake))
  )
  expect_equal(length(est_fact_2$internal_parameters$adaptive_weight[[1]]$log), 2 * choose(5,2))
  
  est_fact_2$internal_parameters$data$X
  expect_equal(length(est_fact_2$internal_parameters$Fmatrix$log), 2 * choose(5,2))
  
  est_fact_2_std <- FactorHet(formula = y ~ state + letter + fake_one * fake_two, design = dta, K = 2, 
    lambda = 1e-2, moderator = ~ mod, 
    init = FactorHet_init(nrep = 1),
    control = FactorHet_control(iterations = 1, log_method = 'standard', return_data = TRUE),
    group = ~ group, task = ~ fake_task)
  expect_equal(est_fact_2_std$internal_parameters$adaptive_weight[[1]]$letter, br_letter)
  expect_equal(est_fact_2_std$internal_parameters$adaptive_weight[[1]]$state, br_state)
  # Confirm STD does not CREATE latent overlapping groups
  expect_equal(length(est_fact_2_std$internal_parameters$Fmatrix), 4)
  
})

test_that('Check LOG recovers correct differences', {
  
  dta <- data.frame(
    state = sample(state.name[1:7], 1000, replace = T),
    letter = sample(letters[1:5], 1000, replace = T)
  )
  dta$group <- rep(sample(1:250, 500, replace = T), each = 2)
  dta$task <- rep(1:500, each = 2)
  dta$prof <- as.vector(sapply(1:500, FUN=function(i){c('l', 'r')[sample(2)]}))
  dta$y <- as.vector(sapply(1:500, FUN=function(i){sample(0:1)}))
  dta$mod <- runif(500, -1, 1)[dta$group]
  dta$fake_task <- paste0(dta$task, dta$prof)
  
  dta$fake_one <- sample(LETTERS[1:6], 1000, replace = T)
  dta$fake_two <- sample(LETTERS[1:3], 1000, replace = T)
  est_fact_2 <- FactorHet(formula = y ~ state + letter + fake_one * fake_two, design = dta, K = 2, 
                          lambda = 1e-2, moderator = ~ mod, 
                          control = FactorHet_control(return_data = TRUE),
                          group = ~ group, task = ~ fake_task)

  D <- do.call('rbind', lapply(est_fact_2$internal_parameters$penalty$D, FUN=function(i){do.call('rbind',i)}))
  orig_beta <- est_fact_2$parameters$beta
  diff_orig <- D %*% orig_beta
  expand_beta <- est_fact_2$internal_parameters$data$basis_M %*% est_fact_2$parameters$nullspace_beta
  expand_diff <- expand_beta[-seq_len(nrow(orig_beta)),]
  expand_log <- expand_diff[-seq_len(nrow(diff_orig)),]
  expand_diff <- expand_diff[1:nrow(diff_orig),]
  
  lengths_factors <- lengths(est_fact_2$internal_parameters$penalty$D)
  lengths_nrow_i <- sapply(est_fact_2$internal_parameters$penalty$D, FUN=function(i){mean(sapply(i, nrow))})
  lengths_nonint <- sum(lengths_factors[1:2])
  
  # Confirm the first p_orig elements of the "expanded" are the original
  expect_equal(as.vector(expand_beta[1:nrow(orig_beta),]), as.vector(orig_beta))
  # Confirm the NON interacted elements exactly recover
  expect_equal(diff_orig[seq_len(lengths_nonint),], expand_diff[seq_len(lengths_nonint),])
  
  manual_positions <- unlist(mapply(lengths_factors[-1:-2], lengths_nrow_i[-1:-2], FUN=function(i,j){
    unlist(lapply(1:i, FUN=function(k){seq_len(j)}))
  }))
  # Adjust for the NON interacted terms
  manual_positions <- which(manual_positions == 1) + lengths_nonint
  # Combine together the "main" diff PLUS the LOG
  reconstitue_diff <- expand_diff + as.matrix(Matrix::t(Matrix::sparseMatrix(i = 1:length(manual_positions), j = manual_positions,
   x = 1, dims = c(length(manual_positions), nrow(expand_diff)))) %*% expand_log)
  # Confirm that D beta = D "main" + D LOG
  expect_equal(diff_orig, reconstitue_diff)
})

