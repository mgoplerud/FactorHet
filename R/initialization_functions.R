
#' Prepare Data
#' 
#' Package data to use in FactorHet. Called internally by `FactorHet_mbo`;
#' rarely called directly by the user.
#' 
#' @keywords internal
prepare_regression_data <- function(formula, design, moderator = NULL, 
  group = NULL, task = NULL, 
  choice_order = NULL, weights = NULL,
  single_intercept = NULL,
  forced_randomize = NULL,
  rare_threshold = NULL, weight_dlist = FALSE,
  rare_verbose = NULL){
  
  default_options <- FactorHet_control()
  if (is.null(forced_randomize)){
    forced_randomize <- default_options$forced_randomize
  }
  if (is.null(rare_threshold)){
    rare_threshold <- default_options$rare_threshold
  }
  if (is.null(rare_verbose)){
    rare_verbose <- default_options$rare_verbose
  }
  
  # Prepare the data
  prep_data <- prepare_formula(fmla_main = formula, fmla_moderator = moderator, weights = weights,
                               group = group, task = task, choice_order = choice_order, design = design)
  
  do_interactions <- prep_data$do_interactions
  design <- prep_data$design
  y <- design[[prep_data$outcome]]
  
  weights <- prep_data$weights
  
  if (!is.null(prep_data$name_weights)){
    if (prep_data$name_weights %in% names(design)){
      stop('Name of "weights" duplicated in design factors or outcome. Rename "weights"!')
    }
    
    design[[prep_data$name_weights]] <- weights
  }
  
  group <- data.frame(design[prep_data$name_group])
  task <- data.frame(design[prep_data$name_task])
  choice_order <- data.frame(design[prep_data$name_choice_order])
  
  
  if (ncol(group) == 0){group <- NULL}else{group <- as.character(as.vector(group[,1]))}
  if (ncol(task) == 0){task <- NULL}else{task <- as.character(as.vector(task[,1]))}
  if (ncol(choice_order) == 0){choice_order <- unique_choice <- NULL}else{
    choice_order <- as.character(as.vector(choice_order[,1]))
    unique_choice <- sort(unique(choice_order))
    if (any(is.na(choice_order))){stop('No missingness permitted for "choice_order"')}
    if (length(unique_choice) != 2){stop('choice_order must have exactly two unique values indexing the "left" and "right" profiles.')}
    
  }
  
  conjoint_names <- prep_data[c('name_group', 'name_task', 'name_choice_order')]
  factor_names <- prep_data$factor_names
  formula_recons <- prep_data$formula_recons
  formula_mod <- prep_data$formula_mod
  
  rm(prep_data)
  
  if (is.null(group)){
    group <- 1:nrow(design)
    null_group <- TRUE
  }else{
    null_group <- FALSE
  }
  unique_group <- unique(group)
  if (any(is.na(group))){stop('No Missingness Allowed for Group')}
  n_G <- length(unique(group))
  
  #Build the moderator variables
  if (is.null(moderator)){
    moderator <- ~ 1
  }
  
  
  concom_W <- create_moderator(design = design, moderator = moderator, group = group, unique_group = unique_group)
  W <- concom_W$W
  args_W <- concom_W$args_W
  ncol_W <- concom_W$ncol_W
  
  #From the design, get the level of each factor
  factor_levels <- lapply(dplyr::select(.data = as.data.frame(design), !! factor_names), FUN=function(i){
    if (any(class(i) == 'factor')){
      return(levels(i))
    }else{
      return(sort(unique(i)))
    }
  })
  
  ordered_factors <- sapply(dplyr::select(.data = as.data.frame(design), !! factor_names), FUN=function(i){
    is.ordered(i)
  })
  
  #Create the sparsity penalty
  #and the ANOVA sum to zero constraints
  penalty_for_regression <- create_penalty(factor_levels, weight_dlist = weight_dlist,
                                           ordered_factors = ordered_factors,
                                           make_interactions = do_interactions)
  X <- create_data(design, penalty_for_regression)
  
  #Remove rare levels (i.e. below rare_threshold # of observations)
  rare_output <- remove_rare_levels(X = X, penalty_for_regression = penalty_for_regression, 
                             rare_threshold = rare_threshold, verbose = rare_verbose)
  
  X <- rare_output$X
  penalty_for_regression <- rare_output$penalty_for_regression
  rare_col <- rare_output$rare
  rare_fmt_col <- rare_output$fmt_rare
  rm(rare_output)
  
  if (length(factor_names) > 1){
    add_restrictions <- c()
    combo_factors <- combn(factor_names, 2)
    for (i in 1:ncol(combo_factors)){
      combo_i <- combo_factors[,i]
      tab_cols <- table(design[[combo_i[1]]], design[[combo_i[2]]])
      if (any(tab_cols == 0)){
        
        ind_restrict <- which(tab_cols == 0, arr.ind = T)
        ind_restrict <- cbind(rownames(tab_cols)[ind_restrict[,'row']],
              colnames(tab_cols)[ind_restrict[,'col']])
        ind_restrict <- paste(combo_i[1], '(', ind_restrict[,1], ')-', combo_i[2], '(', ind_restrict[,2] ,')', sep = '')
        add_restrictions <- c(add_restrictions, ind_restrict)
      }
    }
    # Are there any additional restrictions?
    new_restrictions <- setdiff(add_restrictions, rare_fmt_col)
    if (length(new_restrictions) > 0){
      message('Adding ', length(new_restrictions), ' additional randomization restrictions because no observations at the intersection of these levels.')
      message(paste0(new_restrictions, collapse='\n'))
      rare_fmt_col <- c(rare_fmt_col, new_restrictions)
    }
  }

  Fmatrix <- penalty_for_regression[["F"]]
  Dlist <- penalty_for_regression$D
  constraint <- penalty_for_regression$constraint
  basis_M <- penalty_for_regression$constraint_basis
  
  term_position <- penalty_for_regression$term_position
  
  coef_names <- penalty_for_regression$coef
  
  
  #Clip any small values below 1e-15 to be zero
  clip_tiny <- floor(-log(.Machine$double.eps)/log(10))
  if (clip_tiny < 15){
    clip_tiny <- 15
  }

  basis_M <- drop0(zapsmall(basis_M, clip_tiny))
  #If choice order provided, DO FORCED CHOICE
  if (!is.null(choice_order)){
    
    if (null_group){stop('Group must be provided for forced choice')}
    forced_output <- make_forced_choice(y = y, X = X, group = group, task = task, 
      choice_order = choice_order, unique_choice = unique_choice,
      weights = weights,
      estimate_intercept = TRUE, #require intercept
      randomize = forced_randomize)
    X <- forced_output$X
    y <- forced_output$y
    weights <- forced_output$weights
    group <- forced_output$group
    task <- forced_output$task
    
    rm(forced_output)
    gc()
    use_forced_choice <- TRUE
  }else{
    use_forced_choice <- FALSE
  }
  
  if (is.null(single_intercept)){#If not provided (default),
    #then do a varying intercept if factorial and
    #single intercept if conjoint (forced choice)
    single_intercept <- use_forced_choice
  }
  group_mapping <- sparseMatrix(i = 1:nrow(X), j = match(group, unique_group), x = 1)
  

  
  
  data_output <- mget(ls())
  class(data_output) <- 'FactorHet_data'
  return(data_output)
}

prepare_deterministic_init <- function(data, K, method, iterations = NULL){
  
  if (!inherits(data, 'FactorHet_data')){
    stop('prepare_deterministic_init needs an object from "prepare_regression_data".')
  }
  if (K == 1){
    out <- matrix(1, nrow = length(data$unique_group))
  }else if (method == 'spectral'){
    out <- spectral_init(data$W, K)
  }else if (method == 'mclust'){
    out <- mclust_init(data$W, K)
  }else if (grepl(method, pattern='^mm')){
    
    main_X <- data$X
    levels_X <- c('(Intercept)', unlist(mapply(names(data$factor_levels), data$factor_levels, SIMPLIFY = FALSE, FUN=function(i,j){paste0(i, '(', j, ')')})))
    if (length(setdiff(levels_X, colnames(main_X))) != 0){
      stop('Missing levels in X?')
    }
    main_X <- main_X[ , colnames(main_X) %in% levels_X]

    if (grepl(method, pattern='^mm_spectral')){
      mm_method <- 'spectral'
    }else if (grepl(method, pattern='^mm_mclust')){
      mm_method <- 'mclust'
    }else{stop('mm_ must have "spectral" or "mclust" after')}
    
    out <- murphy_murphy_initialize(y = as.numeric(data$y), X = main_X, W = data$W, K = K,
      group_mapping = data$group_mapping, method = mm_method,
      weights = data$weights, iterations = iterations,
      probabilistic = grepl(method, pattern='prob'))

  }else{stop('Invalid deterministic init method')}
  
  out <- list('group_E.prob' = data.frame(
    group = data$unique_group, out))
  
  return(out)
}

#' @importFrom stats sd kmeans
spectral_init <- function(W, K){
  if (ncol(W) == 1 & all(W[,1] == 1)){
    stop('Spectral initialization with no moderators is not permitted.')
  }
  if (!requireNamespace('FNN', quietly = TRUE) & !requireNamespace('RSpectra', quietly = TRUE)){
    stop('Spectral initialization requires "FNN" and "RSpectra" installed.')
  }

  scale_W <- apply(W, MARGIN = 2, FUN=function(i){
    if (sd(i) == 0){
      return(rep(0, length(i)))
    }else{
      return(scale(i))
    }
  })
  cov_W <- var(scale_W)
  diag(cov_W)[which(diag(cov_W) == 0)] <- 1
  cov_W <- eigen(cov_W)
  scale_W <- scale_W %*% (cov_W$vectors) %*% diag(x = 1/sqrt(cov_W$values))
  
  knn_W <- FNN::get.knn(scale_W, k=10)
  knn_W <- merge(reshape2::melt(knn_W$nn.index, value.name = 'index'),
                     reshape2::melt(knn_W$nn.dist, value.name = 'dist'), 
                     by = c('Var1', 'Var2'))
  graph_W <- sparseMatrix(i = knn_W$Var1, j = knn_W$index, x = 1, dims = c(nrow(scale_W), nrow(W)))
  graph_W <- graph_W + t(graph_W)
  graph_W@x[graph_W@x != 2] <- 0
  # graph_W@x[graph_W@x == 1] <- 2
  graph_W <- drop0(graph_W)
  graph_W <- graph_W/2
  graph_W <- as(graph_W, 'dgTMatrix')
  graph_W <- with(attributes(graph_W), cbind(i, j) + 1)
  knn_W$id <- paste(knn_W$Var1, knn_W$index, sep = '-')
  rotate_K <- knn_W
  rotate_K$id <- paste(knn_W$index, knn_W$Var1, sep = '-')
  knn_W <- rbind(knn_W, rotate_K)
  weight <- exp(-1/2 * knn_W$dist[match(paste(graph_W[,1], graph_W[,2], sep = '-'), knn_W$id)]^2)
  graph_W <- sparseMatrix(i = graph_W[,1], j = graph_W[,2], x = weight,
                          dims = c(nrow(scale_W), nrow(W)))
  diag(graph_W)[] <- 1
  #Un-normalized graph Laplacian
  graph_L <- Diagonal(x = rowSums(graph_W)) - graph_W
  method <- 'ShiMalik'    
  if (method == 'ShiMalik'){
    L_norm <- Diagonal(x= 1/rowSums(graph_W)) %*% graph_L 
    eigen_L <- RSpectra::eigs(A = L_norm, k = K, which = 'SM')$vectors
  }
  #Separate "real" and imaginary component.
  #Shouldn't be much imaginary but...
  eigen_L <- cbind(Re(eigen_L), Im(eigen_L))
  eigen_only <- eigen_L
  
  eigen_L <- cbind(scale_W, eigen_L)
  
  rotate_L <- apply(eigen_L, MARGIN = 2, FUN=function(i){
    if (sd(i) == 0){
      return(rep(0, length(i)))
    }else{
      return(scale(i))
    }
  })
  rotate_L <- rotate_L[, which(apply(rotate_L, MARGIN = 2, FUN=function(i){all(abs(i) < 1e-6)}) != TRUE)]
  
  cov_W <- var(rotate_L)
  diag(cov_W)[which(diag(cov_W) == 0)] <- 1
  cov_W <- eigen(cov_W)
  cov_W$values[cov_W$values < 0] <- 0
  fmt_eigen <- ifelse(cov_W$values <= 0, 0, 1/sqrt(cov_W$values))
  rotate_L <- rotate_L %*% (cov_W$vectors) %*% diag(x = fmt_eigen)
  rotate_L <- rotate_L[, which(apply(rotate_L, MARGIN = 2, FUN=function(i){all(abs(i) < 1e-6)}) != TRUE)]
  
  eigen_kmeans <- kmeans(eigen_only, centers = K, 
   iter.max = 50, 
   nstart = 5000)
  
  group_E.prob <- as.matrix(sparseMatrix(i = 1:nrow(W), j = eigen_kmeans$cluster, x = 1))

  return(group_E.prob)
}

mclust_init <- function(W, K){
  if (requireNamespace('mclust', quietly = TRUE)){
    hcVVV <- mclust::hcVVV
    classif <- as.vector(mclust::hclass(mclust::hc(data = W, use = "SVD"), G = K))
    group_E.prob <- sparseMatrix(i = 1:nrow(W), j = classif, x = 1)
    group_E.prob <- as.matrix(group_E.prob)
  }else{
    stop('"mclust" initialization requires "mclust" installed.')
  }
  return(group_E.prob)
}

# Initialization Based on Murphy/Murphy (2020); see the appendix of Goplerud et
# al. for discussion.
#' @importFrom stats median
murphy_murphy_initialize <- function(y, X, W, K, group_mapping, weights,
                                     method = 'mclust', iterations = NULL,
                                     probabilistic = FALSE){
  if (is.null(iterations)){
    if (probabilistic){iterations <- 100}else{iterations <- 50}
  }
  if (K == 1){stop('No initalization needed or permitted if K = 1')}
  # Step 1: Cluster based on moderators 
  if (method == 'mclust'){
    if (requireNamespace('mclust', quietly = TRUE)){
      if (ncol(W) == 1 & all(W[,1] == 1)){
        classif <- sample(1:K, nrow(W), replace = T)
      }else{
        hcVVV <- mclust::hcVVV
        classif <- as.vector(mclust::hclass(mclust::hc(data = W, use = "SVD"), G = K))
      }
    }else{
      stop('"mclust" initialization requires "mclust" installed.')
    }
  }else if (method == 'spectral'){
    classif <- spectral_init(W, K)
    classif <- rowSums(classif %*% Diagonal(x = 1:K))
  }else{stop('Invalid method to Murphy/Murphy')}
  init_classif <- classif
  
  if (probabilistic){
    init_classif <- classif <- as.matrix(sparseMatrix(i = 1:nrow(W), j = classif, x = 1))
  }
  # Step 2a: Loop over logistic regression and reclassify
  counter <- 0
  converge <- FALSE
  while (counter < iterations & converge == FALSE){
    
    if (probabilistic){
      group_E.prob <- classif
    }else{
      group_E.prob <- as.matrix(sparseMatrix(i = 1:nrow(W), j = classif, x = 1))
    }
    
    obs.E.prob <- apply(group_E.prob, MARGIN = 2, FUN=function(i){as.vector(group_mapping %*% i)})
    
    if (probabilistic){
      
      beta <- sapply(list_from_cols(obs.E.prob), FUN=function(k){
        simple_logit(y = y, X = X, obs.E.prob = k, iterations = 15,
                     weights = weights,
                     beta_method = 'cpp', beta_cg_it = 10, prec_ridge = 1/4)
      })
      
    }else{
      beta <- sapply(list_from_cols(obs.E.prob), FUN=function(k){
        simple_logit(y = y[which(k == 1)], X = X[which(k == 1),], iterations = 25, 
                     weights = weights[which(k == 1)],
                     beta_method = 'cpp', beta_cg_it = 10, prec_ridge = 1/4)
      })
    }
    rownames(beta) <- NULL
    pmat <- plogis(as.matrix(X %*% beta))
    brier <- apply(pmat, MARGIN = 2, FUN=function(i){ifelse(y == 1, log(i), log(1-i))})
    sum_brier <- t(group_mapping) %*% brier
    
    if (probabilistic){
      new_classif <- softmax_matrix(as.matrix(sum_brier))
      
      change <- max(abs(classif - new_classif))
      median_change <- median(abs(classif - new_classif))
      if (change < 1e-2){
        converge <- TRUE; break
      }
    }else{
      
      new_classif <- apply(sum_brier, MARGIN = 1, which.max)
      change <- sum(classif != new_classif)
      if (change == 0){converge <- TRUE; break}
      
    }
    counter <- counter + 1
    classif <- new_classif
  }
  
  if (converge == FALSE){
    message('Murphy/Murphy initalization did not fully converge.')
    if (probabilistic){
      message(paste0('The largest change in any cluster probability was ', round(change, 2)))
      message(paste0('The median change in cluster probability was ', round(median_change, 3)))
    }else{
      message(paste0(change, ' observations out of ', length(classif), ' changed cluster membership at the final iteration.'))
    }
  }
  
  if (probabilistic){
    group_E.prob <- as.matrix(new_classif)
    
  }else{
    group_E.prob <- sparseMatrix(i = 1:nrow(W), j = classif, x = 1)
    group_E.prob <- as.matrix(group_E.prob)
    
  }
  return(group_E.prob)
}