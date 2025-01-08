context('Test grouped')

if (isTRUE(as.logical(Sys.getenv("CI")))){
  # If on CI
  NITER <- 2
  env_test <- "CI"
}else if (!identical(Sys.getenv("NOT_CRAN"), "true")){
  # If on CRAN
  NITER <- 2
  env_test <- "CRAN"
  set.seed(2)
}else{
  # If on local machine
  NITER <- 2000
  env_test <- 'local'
}

test_that('Runs with groups and tasks; post-estimation functions work', {

  N <- 1000
  dta <- data.frame(
    state = sample(state.name[1:4], N, replace = T),
    letter = sample(letters[1:3], N, replace = T)
  )
  dta$group <- rep(sample(1:(N/4), N/2, replace = T), each = 2)
  dta$task <- rep(1:(N/2), each = 2)
  dta$prof <- as.vector(sapply(1:(N/2), FUN=function(i){c('l', 'r')[sample(2)]}))
  dta$y <- as.vector(sapply(1:(N/2), FUN=function(i){sample(0:1)}))
  dta$mod <- runif((N/2), -1, 1)[dta$group]
  
  est_simple <- FactorHet(formula = y ~ state + letter, design = dta, K = 2, 
    lambda = 1e-2, moderator = ~ mod, 
    group = ~ group, task = ~ task, choice_order = ~ prof)
  #Check monotonic increase
  expect_gte(min(diff(logLik(est_simple, 'log_posterior_seq'))), -sqrt(.Machine$double.eps))
  #Check that prediction works  
  test_pred <- tryCatch(predict(object = est_simple), error = function(e){NULL})
  expect_false(is.null(test_pred))
  #Check that conjoint style plot works
  test_plot <- tryCatch(cjoint_plot(est_simple), error = function(e){NULL})
  expect_false(is.null(test_plot))
  #Check that getting marginal treatment effects works
  test_mfe <- tryCatch(AME(est_simple, verbose = FALSE), error = function(e){NULL})
  expect_false(is.null(test_mfe))
  #Same for ACE, AMIE
  test_mfe <- tryCatch(ACE(est_simple, verbose = FALSE, baseline = list(state = 'Arizona', letter = 'c')), error = function(e){NULL})
  expect_false(is.null(test_mfe))
  test_mfe <- tryCatch(AMIE(est_simple, verbose = FALSE, baseline = list(state = 'Arizona', letter = 'c')), error = function(e){NULL})
  expect_false(is.null(test_mfe))
})

test_that('Fails if incorrectly specified groups, task, choice', {
  
  N <- 200
  dta <- data.frame(
    state = sample(state.name[1:4], N, replace = T),
    letter = sample(letters[1:3], N, replace = T)
  )
  dta$group <- rep(sample(1:(N/4), N/2, replace = T), each = 2)
  dta$task <- rep(1:(N/2), each = 2)
  dta$prof <- as.vector(sapply(1:(N/2), FUN=function(i){c('l', 'r')[sample(2)]}))
  dta$y <- as.vector(sapply(1:(N/2), FUN=function(i){sample(0:1)}))
  dta$mod <- runif((N/2), -1, 1)[dta$group]
  
  expect_error(
    FactorHet(formula = y ~ state + letter, design = dta, K = 2, 
            lambda = 1e-2, moderator = ~ mod, 
            group = ~ group, task = ~ task), regex = 'More than one'
  )
  expect_error(
    FactorHet(formula = y ~ state + letter, design = dta, K = 2, 
              lambda = 1e-2, moderator = ~ mod, 
              group = ~ group), regex = 'More than one'
  )
  # Run conjoint as factorial by treating each as a task
  dta$fake_task <- paste0(dta$task, dta$prof)
  expect_s3_class(
    FactorHet(formula = y ~ state + letter, design = dta, K = 2, 
              lambda = 1e-2, moderator = ~ mod, 
              control = FactorHet_control(iterations = 1),
              initialize = FactorHet_init(nrep = 1, short_EM_it = 1),
              group = ~ group, task = ~ fake_task), 'FactorHet'
  )
  
})
