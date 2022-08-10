context('Test that Average Marginal (Interaction) Effects can be calculated')

test_that('Test that AME runs', {
  
  dta <- data.frame(
    state = sample(state.name[1:4], 50, replace = T),
    letter = sample(letters[1:3], 50, replace = T),
    region = sample(state.region, 50, replace = T)
  )
  dta$mod <- runif(nrow(dta), -1, 1)
  dta$y1 <- rbinom(nrow(dta), 1, plogis(runif(5, -1 , 1)[match(dta$state, state.name[1:5])] + runif(5, -1, 1)[match(dta$letter, letters)]))
  dta$y2 <- rbinom(nrow(dta), 1, plogis(runif(5, -1 , 1)[match(dta$state, state.name[1:5])] + runif(5, -1, 1)[match(dta$letter, letters)]))
  
  dta$y <- ifelse(rbinom(nrow(dta), 1, plogis(dta$mod)) == 1, dta$y1, dta$y2)
  
  est_simple <- FactorHet(formula = y ~ state + letter + region, design = dta, 
    init = FactorHet_init(nrep = 2),
    K = 2, lambda = 1e-3, 
    moderator = ~ mod, 
    control =  FactorHet_control(iterations = 10))
  
  est_AME <- marginal_AME(est_simple)
  est_ACE <- marginal_ACE(est_simple, baseline = list(letter = 'a', state = 'Arizona'))
  est_AMIE <- marginal_AMIE(est_simple, baseline = list(letter = 'a', state = 'Arizona'))
  est_AMIE <- marginal_AMIE(est_simple)
  
  expect_s3_class(est_AME, 'FactorHet_vis')
  expect_s3_class(est_ACE, 'FactorHet_vis')
  expect_s3_class(est_AMIE, 'FactorHet_vis')
})

test_that('Test that AME works with restriction', {
  
  NLETTER <- sample(c(3,4,5,2,7,8), 1)
  dta <- data.frame(
    state = sample(state.name[1:4], 1000, replace = T),
    letter = sample(letters[1:NLETTER], 1000, replace = T),
    extra = sample(state.abb[1:5], 1000, replace = TRUE)
  )
  induce <- which(dta$state == 'Alabama' & dta$letter == 'a')
  dta$state[induce] <- sample(state.name[2:4], length(induce), replace = T)
  
  dta$group <- rep(sample(1:250, 500, replace = T), each = 2)
  dta$task <- rep(1:500, each = 2)
  dta$prof <- as.vector(sapply(1:500, FUN=function(i){c('l', 'r')[sample(2)]}))
  dta$y <- as.vector(sapply(1:500, FUN=function(i){sample(0:1)}))
  dta$mod <- runif(500, -1, 1)[dta$group]
  
  est_simple <- FactorHet(formula = y ~ state + letter + extra, design = dta, K = 2, 
    lambda = 1e-3, moderator = ~ mod, 
    group = ~ group, task = ~ task, 
    choice_order = ~ prof,
    initialize = FactorHet_init(nrep = 1),
    control =  FactorHet_control(return_data = TRUE, iterations = 15))
  
  expect_null(print(est_simple))
  expect_false(is.null(est_simple$internal_parameters$rare$rare_fmt_col))
  
  simple_AME <- suppressMessages(
    suppressWarnings(marginal_AME(est_simple, verbose = FALSE, ignore_restrictions = TRUE))
  )
  #Do AME manually
  manual_AME <- data.frame()
  for (v in state.name[2:4]){
    
    copy_dta <- dta
    copy_dta$state <- as.character(copy_dta$state)
    baseline_dta <- mod_dta <- copy_dta
    # For the RIGHT profile, set the "treatment" to v
    # And then to baseline and compare
    # Marginalizes over all other factors
    mod_dta$state <- ifelse(mod_dta$prof == 'r', v, mod_dta$state)
    baseline_dta$state <- ifelse(baseline_dta$prof == 'r', 'Alabama', baseline_dta$state)
    AME_right <- colMeans(predict(est_simple, newdata = mod_dta, by_cluster = TRUE) -
      predict(est_simple, newdata = baseline_dta, by_cluster = TRUE))
    # Do the same for LEFT
    baseline_dta <- mod_dta <- copy_dta
    mod_dta$state <- ifelse(mod_dta$prof == 'l', v, mod_dta$state)
    baseline_dta$state <- ifelse(baseline_dta$prof == 'l', 'Alabama', baseline_dta$state)
    # Note that this is Pr(Y_i = 1 | T) - Pr(Y_i = | BASE)
    # So if we +1/-1  - Pr(Y_i = 0 | T) + Pr(Y_i = 0 | Base)
    # So if we take negative, Pr(Y_i 0 | T) - Pr(Y_i = 0 | Base)
    AME_left <- colMeans(predict(est_simple, newdata = mod_dta, by_cluster = TRUE) -
      predict(est_simple, newdata = baseline_dta, by_cluster = TRUE))
    AME <- (AME_right - AME_left)/2
    manual_AME <- rbind(manual_AME, 
      data.frame(state = v, marginal_effect = AME, cluster = c(1,2)))
  }
  
  implemented_AME <- subset(simple_AME$data, !baseline & factor == 'state')
  order_AME <- apply(implemented_AME[, c('level', 'cluster')], MARGIN = 1, paste, collapse =' ')

  manual_AME$order <- apply(manual_AME[, c('state', 'cluster')], MARGIN = 1, paste, collapse =' ')
  
  expect_equivalent(
    implemented_AME$marginal_effect,
    manual_AME$marginal_effect[match(order_AME, manual_AME$order)]
  )

  #Do AME manually *with* restricted randomization
  simple_AME <- marginal_AME(est_simple)
  manual_AME <- data.frame()
  for (v in state.name[2:4]){
    
    baseline_dta_right <- dta
    baseline_dta_right$invalid <- unsplit(lapply(split(baseline_dta_right, paste(baseline_dta_right$group, baseline_dta_right$task)),
     FUN=function(i){
       rep(any(i$prof == 'r' & i$letter == 'a'), 2)
     }), paste(baseline_dta_right$group, baseline_dta_right$task))
    baseline_dta_right <- baseline_dta_right[!baseline_dta_right$invalid, ]
    # dplyr_filter <- baseline_dta_right %>% dplyr::group_by(group, task) %>% 
    #   dplyr::mutate(invalid = any(prof == 'r' & letter == 'a')) %>%
    #   dplyr::filter(!invalid)
    # expect_equivalent(baseline_dta_right, dplyr_filter)
    baseline_dta_right$state <- ifelse(baseline_dta_right$prof == 'r', 'Alabama', as.character(baseline_dta_right$state))
    mod_dta_right <- baseline_dta_right
    mod_dta_right$state <- ifelse(mod_dta_right$prof == 'r', v, as.character(mod_dta_right$state))
    mod_dta_right$state <- factor(mod_dta_right$state, levels = state.name[1:4])
    baseline_dta_right$state <- factor(baseline_dta_right$state, levels = state.name[1:4])
    
    AME_right <- colMeans(predict(est_simple, newdata = mod_dta_right, by_cluster = TRUE) -
                            predict(est_simple, newdata = baseline_dta_right, by_cluster = TRUE))
    
    baseline_dta_left <- dta
    baseline_dta_left$invalid <- unsplit(lapply(split(baseline_dta_left, 
      paste(baseline_dta_left$group, baseline_dta_left$task)),
       FUN=function(i){
         rep(any(i$prof == 'l' & i$letter == 'a'), 2)
       }), paste(baseline_dta_left$group, baseline_dta_left$task))
    baseline_dta_left <- baseline_dta_left[!baseline_dta_left$invalid, ]
    
    baseline_dta_left$state <- ifelse(baseline_dta_left$prof == 'l', 'Alabama', as.character(baseline_dta_left$state))
    mod_dta_left <- baseline_dta_left
    mod_dta_left$state <- ifelse(mod_dta_left$prof == 'l', v, as.character(mod_dta_left$state))
    mod_dta_left$state <- factor(mod_dta_left$state, levels = state.name[1:4])
    baseline_dta_left$state <- factor(baseline_dta_left$state, levels = state.name[1:4])
    
    AME_left <- colMeans(predict(est_simple, newdata = mod_dta_left, by_cluster = TRUE) -
                           predict(est_simple, newdata = baseline_dta_left, by_cluster = TRUE))
    AME <- ( AME_right - AME_left )/2
    
    manual_AME <- rbind(manual_AME, data.frame(state = v, marginal_effect = AME, cluster = c(1,2)))
  }
  
  implemented_AME <- subset(simple_AME$data, !baseline & factor == 'state')
  order_AME <- apply(implemented_AME[, c('level', 'cluster')], MARGIN = 1, paste, collapse =' ')
  manual_AME$order <- apply(manual_AME[, c('state', 'cluster')], MARGIN = 1, paste, collapse =' ')
  
  expect_equivalent(
    implemented_AME$marginal_effect,
    manual_AME$marginal_effect[match(order_AME, manual_AME$order)]
  )
})

