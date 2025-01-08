context("Test ordered categories")

if (isTRUE(as.logical(Sys.getenv("CI")))){
  # If on CI
  NITER <- 2
  env_test <- "CI"
}else if (!identical(Sys.getenv("NOT_CRAN"), "true")){
  # If on CRAN
  NITER <- 2
  env_test <- "CRAN"
  set.seed(14)
}else{
  # If on local machine
  NITER <- 2000
  env_test <- 'local'
}

test_that("Ordered fusion runs and limits F accordingly", {
  
  dta <- data.frame(
    state = sample(state.name[1:4], 100, replace = T),
    letter = sample(letters[1:25], 100, replace = T),
    stringsAsFactors = FALSE
  )
  dta$y <- rbinom(nrow(dta), 1, plogis(runif(5, -1 , 1)[match(dta$state, state.name[1:4])] + sort(runif(26, -1, 1))[match(dta$letter, letters)]))
  dta$letter <- factor(dta$letter, levels = letters[1:25], ordered = TRUE)
  
  est_simple <- FactorHet(y ~ state + letter + . * ., design = dta, 
                          K = 2, lambda = 1e-3, 
                          initialize = FactorHet_init(nrep = 1),
                          control = FactorHet_control(return_data = TRUE, 
                                                      iterations = 50))
  expect_equivalent(est_simple$internal_parameters$ordered_factors, c(FALSE, TRUE))
  expect_equal(names(est_simple$internal_parameters$penalty[["F"]]$letter),
               sapply(1:24, FUN=function(i){paste0(letters[i], '-', letters[i+1])}))


  dta$letter <- factor(dta$letter, levels = levels(dta$letter), ordered = FALSE)
  
  est_simple <- FactorHet(y ~ state + letter + . * ., design = dta, 
                          K = 2, lambda = 1e-3, 
                          initialize = FactorHet_init(nrep = 1),
                          control = FactorHet_control(return_data = TRUE, 
                                                      iterations = 50))
  expect_equivalent(est_simple$internal_parameters$ordered_factors, c(FALSE, FALSE))
  expect_equal(length(est_simple$internal_parameters$penalty[["F"]]$letter),
               ncol(combn(25, 2)))
  
})
