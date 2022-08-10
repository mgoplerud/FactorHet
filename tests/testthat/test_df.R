context("Degrees of Freedom")


test_that(desc = "Check Rcpp and Base R agree for trace method", {
  
  base_R_kc <- function(X, E.ridge, beta){
    df_per_component <- rep(NA, ncol(beta))
    
    for (k in 1:ncol(beta)){
      p_k <- plogis(as.vector(X %*% beta[,k]))
      H_k <- t(X) %*% Diagonal(x = p_k * (1-p_k)) %*% X
      
      trace_df <- sum(diag(solve(H_k + E.ridge[[k]], H_k)))
      df_per_component[k] <- trace_df
    }
    return(df_per_component)
  }
  
  N <- 1000
  p <- 25
  K <- 2
  
  X <- rsparsematrix(nrow = N, ncol = p, density = 0.8)
  beta <- matrix(rnorm(p * K), ncol = K)
  list_ridge <- lapply(1:K, FUN=function(i){drop0(rWishart(n = 1, df = ncol(X) + 5, Sigma = diag(ncol(X)))[,,1])})
  
  est_base <- base_R_kc(X = X, beta = beta, E.ridge = list_ridge)
  est_cpp <- sapply(1:K, FUN=function(i){trace_df_cpp(X = X, beta = beta[,i], ridge = list_ridge[[i]])})
  
  expect_equal(est_base, est_cpp)
})

test_that('Check that with ONLY intercept, model runs and gives right DF', {
  
  dta <- data.frame(
    state = sample(state.name[1:2], 50, replace = T),
    letter = sample(letters[1:2], 50, replace = T)
  )
  dta$y <- rbinom(nrow(dta), 1, plogis(runif(5, -1 , 1)[match(dta$state, state.name[1:5])] + runif(5, -1, 1)[match(dta$letter, letters)]))

  
  est_simple <- FactorHet(formula = y ~ state + letter, moderator = ~ 1, 
    K = 1, lambda = 10^5, 
    design = dta, control = FactorHet_control(return_data = TRUE,
      iterations = 5))
  
  expect_equal(est_simple$information_criterion$df, rep(1))
  
  est_simple <- FactorHet(formula = y ~ state + letter, design = dta, 
    K = 2, lambda = 10, 
    initialize = FactorHet_init(nrep = 1),
    control = FactorHet_control(iterations = 5, single_intercept = TRUE))
  
  expect_gte(min(diff(logLik(est_simple, 'log_posterior_seq'))), -1e-7)
  
  expect_equal(est_simple$information_criterion$df_beta, rep(1,1))
  expect_equal(est_simple$information_criterion$df, rep(2,1))

  est_simple <- FactorHet(formula = y ~ state + letter, design = dta, 
    K = 2, lambda = 10, 
    initialize = FactorHet_init(nrep = 1),
    control = FactorHet_control(iterations = 5, single_intercept = FALSE))
  expect_gte(min(diff(logLik(est_simple, 'log_posterior_seq'))), -1e-7)
  expect_equal(est_simple$information_criterion$df_beta, rep(2,1))
  expect_equal(est_simple$information_criterion$df, rep(3,1))
  
})
