context("Misc. Tests")

test_that('Test Double Resrictions', {
  
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
  df$outcome <- rbinom(nrow(df), 1, 0.5)
  
  fit_FH <- FactorHet(
    formula = outcome ~ age + edu + job, design = df,
    K = 1, lambda = 1e-3,
  )
  
  
  fit_default <- marginal_AME(fit_FH, baseline = list('job' = 'banker'))
  fit_manual_res <- marginal_AME(fit_FH, baseline = list('job' = 'banker'),
    design = subset(df, (age != '18-29') & !(edu %in% c('Col', 'HS')))
  )
  
  expect_equal(fit_default$data, fit_manual_res$data)
  expect_equal(fit_default$plot$data, fit_manual_res$plot$data)
  
  fit_default <- marginal_AME(fit_FH, baseline = list('edu' = 'Col'))
  fit_manual_res <- marginal_AME(fit_FH, baseline = list('edu' = 'Col'),
    design = subset(df, job != 'doctor')
  )
  
  expect_equal(fit_default$data, fit_manual_res$data)
  expect_equal(fit_default$plot$data, fit_manual_res$plot$data)
  
})

