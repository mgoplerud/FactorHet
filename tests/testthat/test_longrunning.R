if (isTRUE(as.logical(Sys.getenv("CI")))){
  # If on CI
  NITER <- 2
  env_test <- "CI"
}else if (!identical(Sys.getenv("NOT_CRAN"), "true")){
  # If on CRAN
  NITER <- 2
  env_test <- "CRAN"
  set.seed(4)
}else{
  # If on local machine
  NITER <- 2000
  env_test <- 'local'
}

test_that('Agrees with FlexMix (K = 2)', {
  
  skip_on_cran()
  skip_on_ci()
  
  dta <- data.frame(
    state = sample(state.name[1:4], 100, replace = T),
    letter = sample(letters[1:3], 100, replace = T)
  )
  dta$y <- rbinom(nrow(dta), 1, plogis(runif(5, -1 , 1)[match(dta$state, state.name[1:5])] + runif(5, -1, 1)[match(dta$letter, letters)]))
  
  for (single_intercept in c(TRUE, FALSE)){
    message(single_intercept)
    
    est_simple <- FactorHet(formula = y ~ state + letter, design = dta, K = 2, lambda = 0, 
                            control = FactorHet_control(return_data = TRUE, 
                                                        single_intercept = single_intercept)
    )
    
    expect_gt(min(diff(logLik(est_simple, 'log_posterior_seq'))), -sqrt(.Machine$double.eps))
    if (TRUE){
      
      fmt_data <- data.frame(y = est_simple$internal_parameters$data$y, fmt_X = as.matrix(est_simple$internal_parameters$data$X))
      
      flexmix <- getFromNamespace('flexmix', 'flexmix')
      FLXMRglmfix <- getFromNamespace('FLXMRglmfix', 'flexmix')
      fm_param <- getFromNamespace('parameters', 'modeltools')
      fm_posterior <- getFromNamespace('posterior', 'modeltools')
      if (single_intercept){
        #Initialize flexmix from FactorHet (otherwise FH does better and hard to get 'agreement')
        est_fm <- suppressWarnings(flexmix(cbind(y, 1-y) ~ 0 + fmt_X.1 + fmt_X.2 + fmt_X.3 + fmt_X.4 + fmt_X.5, k = 2, 
                                                    data = fmt_data, cluster = as.matrix(est_simple$posterior$posterior[,-1]),
                                                    model = list(FLXMRglmfix(fixed = ~ fmt_X.6, family = 'binomial'))))
      }else{
        #Initialize flexmix from FactorHet (otherwise FH does better and hard to get 'agreement')
        est_fm <- suppressWarnings(flexmix(cbind(y, 1-y) ~ fmt_X.1 + fmt_X.2 + fmt_X.3 + fmt_X.4 + fmt_X.5, k = 2, 
                                                    data = fmt_data,  cluster = as.matrix(est_simple$posterior$posterior[,-1]),
                                                    model = list(FLXMRglmfix(family = 'binomial'))))
        
      }
      
      ll_fm <- est_fm@logLik
      
      expect_equivalent(as.numeric(ll_fm), logLik(est_simple), tolerance = 0.01)
      
      param_fm <- fm_param(est_fm)
      param_fm <- rbind(param_fm[-1,], param_fm[1,])
      
      cor_all <- cor(sort(as.vector(param_fm)), sort(as.vector(est_simple$parameters$nullspace_beta)))
      expect_gt(cor_all, 0.90)
      
      cor_posterior <- abs(cor(fm_posterior(est_fm)[,1], est_simple$posterior$posterior$group_1))
      expect_gt(cor_posterior, 0.90)
      
      cor_11 <- cor(param_fm[,1], est_simple$parameters$nullspace_beta[,1])
      cor_12 <- cor(param_fm[,1], est_simple$parameters$nullspace_beta[,2])
      
      if (cor_11 < cor_12){
        param_fm <- param_fm[,2:1]
      }
      
      cor_1 <- cor(param_fm[,1], est_simple$parameters$nullspace_beta[,1])
      cor_2 <- cor(param_fm[,2], est_simple$parameters$nullspace_beta[,2])
      
      expect_gt(cor_1, 0.85)
      expect_gt(cor_2, 0.85)
    }
    
  }
  
})

test_that('Test Double Restrictions', {
  
  age <- sample(c('18-29', '30-39', '40-49'), 10^4, replace = T)
  edu <- sample(c('HS', 'Col', "PostGrad"), 10^4, replace = T)
  job <- sample(c('banker', 'janitor', 'doctor', 'gardener'), 10^4, replace = T)
  
  df <- data.frame(age, edu, job, stringsAsFactors = F)
  df <- subset(df,
    ! (job %in% 'doctor' & edu %in% c('HS', 'Col') )
  )
  df <- subset(df, 
    ! (job %in% 'banker' & age %in% '18-29' )
  )
  with(df, table(age, edu))
  with(df, table(age, job))
  with(df, table(edu, job))
  df$outcome <- rbinom(nrow(df), 1, 0.5)
  
  fit_FH <- FactorHet(
    formula = outcome ~ age + edu + job, design = df,
    K = 1, lambda = 1e-3
  )
  
  
  fit_default <- AME(fit_FH, baseline = list('job' = 'banker'))
  fit_manual_res <- AME(fit_FH, baseline = list('job' = 'banker'),
    design = subset(df, (age != '18-29') & !(edu %in% c('Col', 'HS')))
  )
  
  expect_equal(fit_default$data, fit_manual_res$data)
  expect_equal(fit_default$plot$data, fit_manual_res$plot$data)
  
  fit_default <- AME(fit_FH, baseline = list('edu' = 'Col'))
  fit_manual_res <- AME(fit_FH, baseline = list('edu' = 'Col'),
    design = subset(df, job != 'doctor')
  )
  
  expect_equal(fit_default$data, fit_manual_res$data)
  expect_equal(fit_default$plot$data, fit_manual_res$plot$data)
  
})

