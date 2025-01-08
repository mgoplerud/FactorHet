context('Check for initialization and packaging data')

if (isTRUE(as.logical(Sys.getenv("CI")))){
  # If on CI
  NITER <- 2
  env_test <- "CI"
}else if (!identical(Sys.getenv("NOT_CRAN"), "true")){
  # If on CRAN
  NITER <- 2
  env_test <- "CRAN"
  set.seed(16)
}else{
  # If on local machine
  NITER <- 2000
  env_test <- 'local'
}

test_that('Check that packaged / non-packaged agree', {
  
  dta <- data.frame(
    state = sample(state.name[1:4], 1000, replace = T),
    letter = sample(letters[1:3], 1000, replace = T)
  )
  
  dta$group <- rep(sample(1:250, 500, replace = T), each = 2)
  dta$task <- rep(1:500, each = 2)
  dta$prof <- as.vector(sapply(1:500, FUN=function(i){c('l', 'r')[sample(2)]}))
  dta$y <- as.vector(sapply(1:500, FUN=function(i){sample(0:1)}))
  dta$mod <- runif(500, -1, 1)[dta$group]
  
  wgt <- abs(rcauchy(n = 500))
  dta$wgt_internal <- wgt[dta$group]
  
  package_data <- prepare_regression_data(formula = y ~ state + letter, design = dta,
    moderator = ~ mod, group = ~ group, weights = ~ wgt_internal,
    task = ~ task, choice_order = ~ prof, 
    )
  
  package_init <- prepare_deterministic_init(data = package_data, 
        K = 2, method = 'mclust')
  
  est_1 <- FactorHet(formula = y ~ state + letter, design = package_data,
            moderator = ~ mod, group = ~ group, 
            lambda = 1e-3, K = 2, weights = ~ wgt_internal,
            control = FactorHet_control(init_method = 'mclust'),
            task = ~ task, choice_order = ~ prof)
  est_2 <- FactorHet(formula = y ~ state + letter, design = dta,
            moderator = ~ mod, group = ~ group, 
            lambda = 1e-3, K = 2, weights = ~ wgt_internal,
            control = FactorHet_control(init_method = 'mclust'),
            task = ~ task, choice_order = ~ prof)
  
  est_3 <- FactorHet(formula = y ~ state + letter, design = dta,
     moderator = ~ mod, group = ~ group, 
     lambda = 1e-3, K = 2, weights = ~ wgt_internal,
     control = FactorHet_control(init_method = package_init),
     task = ~ task, choice_order = ~ prof)
  
  est_1$internal_parameters$timing <- NULL
  est_1$internal_parameters$control$init_method <- NULL
  est_2$internal_parameters$timing <- NULL
  est_2$internal_parameters$control$init_method <- NULL
  est_3$internal_parameters$timing <- NULL
  est_3$internal_parameters$control$init_method <- NULL
  
  est_1$internal_parameters$refit$make_X_refit$add_col <- NULL
  est_2$internal_parameters$refit$make_X_refit$add_col <- NULL
  est_3$internal_parameters$refit$make_X_refit$add_col <- NULL

  est_1$internal_parameters$control$lambda_scale <- NULL
  est_2$internal_parameters$control$lambda_scale <- NULL
  est_3$internal_parameters$control$lambda_scale <- NULL
  
  # Equivalent except for timing
  expect_true(all.equal(est_1, est_2))
  expect_true(all.equal(est_1, est_3))
})

test_that('Check that all methods for weighted initialization run', {
  
  dta <- data.frame(
    state = sample(state.name[1:4], 1000, replace = T),
    letter = sample(letters[1:3], 1000, replace = T)
  )
  
  dta$group <- rep(sample(1:250, 500, replace = T), each = 2)
  dta$task <- rep(1:500, each = 2)
  dta$prof <- as.vector(sapply(1:500, FUN=function(i){c('l', 'r')[sample(2)]}))
  dta$y <- as.vector(sapply(1:500, FUN=function(i){sample(0:1)}))
  dta$mod <- runif(500, -1, 1)[dta$group]
  
  wgt <- abs(rcauchy(n = 500))
  dta$wgt_internal <- wgt[dta$group]
  
  package_data <- prepare_regression_data(formula = y ~ state + letter, design = dta,
                moderator = ~ mod, group = ~ group, weights = ~ wgt_internal,
                task = ~ task, choice_order = ~ prof
  )
  
  for (m in c('mm_mclust_prob', 'mm_mclust', 'mclust',
              'mm_spectral', 'mm_spectral_prob', 'spectral')){
    expect_type(suppressWarnings(prepare_deterministic_init(package_data, K = 3, method = m)), 'list')
  }
  
})
