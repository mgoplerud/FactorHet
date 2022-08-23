# Internal function for fitting analysis
# Only called via FactorHet; see that documentation
#' @import tictoc rlang Matrix
#' @importFrom stats rexp var runif rnorm
#' @importFrom dplyr select
#' @importFrom utils combn
#' @importFrom methods as
EM_analysis <- function(formula, design, K, lambda, weights = NULL,
    moderator = NULL, group = NULL, task = NULL, choice_order = NULL,
    control = FactorHet_control()){
  
  tic.clear()
  tic.clearlog()
  tic('Initial Preparation')
  
  group <- enquo(group)
  task <- enquo(task)
  choice_order <- enquo(choice_order)
  formula_weight <- paste(as.character(weights), collapse = ' ')
  
  if (!inherits(control, 'FactorHet_control')){
    stop('control must be made using FactorHet_control()')
  }
  if (lambda < 0){stop('lambda must be greater than or equal to zero.')}
  # Load arguments from control
  iterations <- control$iterations
  maxit_pi <- control$maxit_pi
  optim_phi_ctrl <- control$optim_phi_controls
  
  gamma <- control$gamma
  repeat_beta <- control$repeat_beta
  tau_truncate <- control$tau_truncate
  
  debug <- control$debug

  tolerance.parameters <- control$tolerance.parameters
  tolerance.logposterior <- control$tolerance.logposterior
  
  log_method <- control$log_method
  if (lambda == 0 & log_method != 'standard'){
    log_method <- 'standard'
    message('Ridge enforces "standard" project method.')
  }
  weight_dlist <- control$weight_dlist
  
  beta_method <- control$beta_method
  beta_cg_it <- control$beta_cg_it
  
  tau_method <- control$tau_method
  tau_stabilization <- control$tau_stabilization
  tau_equality_threshold <- control$tau_equality_threshold
  
  single_intercept <- control$single_intercept
  
  do_SQUAREM <- control$do_SQUAREM

  quiet_tictoc <- control$quiet_tictoc

  adaptive_weight <- control$adaptive_weight
  #If adaptive weight is *NOT* none or basic,
  #get Bondell & Reich standardization weights
  if (identical(adaptive_weight, 'none') | identical(adaptive_weight, 'basic')){
    do_BR <- FALSE
  }else{
    do_BR <- TRUE
  }
  if (control$override_BR){
    do_BR <- FALSE
  }

  if (debug){
    check_all_ll <- TRUE
  }else{
    check_all_ll <- FALSE
  }
  
  #If the "design" is a data.frame, do the usual parsing.
  #If comes from prepare_regression_data (e.g. for repeated reuse)
  #and is of class FactorHet_data -> then simply parse to environment
  if (!inherits(design, 'FactorHet_data')){
    # Parse the formula, omit missing data
    prep_data <- prepare_formula(fmla_main = formula, fmla_moderator = moderator, weights = weights,
          group = group, task = task, choice_order = choice_order, design = design)
    
    design <- prep_data$design
    do_interactions <- prep_data$do_interactions
    y <- design[[prep_data$outcome]]
    y <- as.numeric(y)
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
    if (is.null(moderator) | K == 1){
      if (K == 1){message('Provided moderators are ignored since K = 1.')}
      moderator <- ~ 1
    }
    
    tic('Prep Moderator')
    concom_W <- create_moderator(design = design, moderator = moderator, group = group, unique_group = unique_group)
    W <- concom_W$W
    args_W <- concom_W$args_W
    ncol_W <- concom_W$ncol_W
    
    ridge_phi <- 1/control$prior_var_phi
    ridge_beta <- 1/control$prior_var_beta
    
    if (grepl(log_method, pattern='^log') & lambda == 0 & ridge_beta == 0){
      message('Cannot Use LOG if lambda == 0 and prior_var_beta == 0; converting to "standard')
      log_method <- 'standard'
    }
    
    if (lambda != 0 & ridge_beta != 0){
      warning('ridge_beta ignored when lambda > 0; only used to stabilize lambda = 0 case.')
    }
    rm(concom_W)
    toc(quiet = quiet_tictoc, log = TRUE)
    
    if (lambda == 0){
      message('Doing Mixture of Regressions (Lambda = 0)')
    }
    
    if (K == 1){
      message('Doing Single Fusion (No Mixture) (K = 1)')
    }
    
    tic('Prepare Design')
    factor_names <- enquo(factor_names)
    
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
    
    #Remove rare combinations [for interactions?] (i.e. below clip_threshold # of observations)
    rare_output <- remove_rare_levels(X = X, penalty_for_regression = penalty_for_regression, 
                                rare_threshold = control$rare_threshold, verbose = control$rare_verbose)
    
    X <- rare_output$X
    penalty_for_regression <- rare_output$penalty_for_regression
    rare_col <- rare_output$rare
    rare_fmt_col <- rare_output$fmt_rare
    rm(rare_output)
    gc()    
    
    if (length(penalty_for_regression$J_names) > 1){
      add_restrictions <- c()
      combo_factors <- combn(penalty_for_regression$J_names, 2)
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
    
    toc(quiet = quiet_tictoc, log = TRUE)
    
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
    tic('Prepare Forced Choice')
    basis_M <- drop0(zapsmall(basis_M, clip_tiny))
    #If choice order provided, DO FORCED CHOICE
    
    if (!is.null(choice_order)){
      
      if (null_group){stop('Group must be provided for forced choice')}
      forced_output <- make_forced_choice(y = y, X = X, group = group, task = task, 
            choice_order = choice_order, unique_choice = unique_choice,
            weights = weights,
            estimate_intercept = TRUE, #require intercept
            randomize = control$forced_randomize)
      X <- forced_output$X
      y <- forced_output$y
      y <- as.numeric(y)
      weights <- forced_output$weights
      if (length(weights) != nrow(X)){stop('weight length mis-aligned.')}
      group <- forced_output$group
      task <- forced_output$task
      
      rm(forced_output)
      use_forced_choice <- TRUE
    }else{
      use_forced_choice <- FALSE
    }
    
    if (is.null(single_intercept)){
      #If not provided (default),
      #then do a varying intercept if factorial and
      #single intercept if conjoint (forced choice)
      single_intercept <- use_forced_choice
    }
    
    group_mapping <- sparseMatrix(i = 1:nrow(X), j = match(group, unique_group), x = 1)
    
    toc(quiet = quiet_tictoc, log = TRUE)
  }else{
    ridge_phi <- 1/control$prior_var_phi
    ridge_beta <- 1/control$prior_var_beta
    
    list2env(design, base::environment())
    gc()
    
    ridge_phi <- 1/control$prior_var_phi
    ridge_beta <- 1/control$prior_var_beta
    
    if (grepl(log_method, pattern='^log') & lambda == 0 & ridge_beta == 0){
      message('Cannot Use LOG if lambda == 0 and prior_var_beta == 0; converting to "standard')
      log_method <- 'standard'
    }
    
    if (lambda != 0 & ridge_beta != 0){
      warning('ridge_beta ignored when lambda > 0; only used to stabilize lambda = 0 case.')
    }
    
  }
  
  # Prepare Weights
  if (any(is.na(weights))){stop('weights may not contain any missing values.')}
  if (any(weights <= 0)){
    stop('"weights" must be non-negative.')
  }
  
  weights_W <- Diagonal(x = 1/colSums(group_mapping)) %*% t(group_mapping) %*% weights
  weights_sq_W <- Diagonal(x = 1/colSums(group_mapping)) %*% t(group_mapping) %*% weights^2
  if (any(abs(as.vector(weights_sq_W - weights_W^2)) > sqrt(.Machine$double.eps))){
    stop('weights are not constant within group (e.g. person) across tasks.')
  }
  std_weights <- 1/sum(weights_W) * nrow(W)
  weights_W <- weights_W * std_weights
  weights <- weights * std_weights
  
  weights_W <- as.vector(weights_W)
  weights <- as.vector(weights)
  
  # print(head(weights_W))
  if (abs(sum(weights_W) - nrow(W)) > sqrt(.Machine$double.eps)){
    warning('Weird Standardization Issue')
  }
  names(weights) <- NULL

  original_X <- X
  
  if (lambda == 0){
    sd_X <- apply(original_X, MARGIN = 2, FUN=function(i){sqrt(sum(i^2))})
    sd_X[sd_X == 0] <- 1
    sd_X[] <- 1
  }
  # Scale lambda based on design size
  orig_lambda <- lambda
  lambda <- control$lambda_scale(lambda, nrow(design))
  
  if (log_method == 'standard'){
    
    #Project into NULLSPACE for regression
    
    X <- X %*% basis_M  
    p_X <- ncol(X)
    
    if (!all(X[,ncol(X)] == 1)){
      stop('Weird alignemnt issue with INTERCEPT not being last column after project. Examine X * basis_M')
    }
    
    tic('Proj F')
    
    Fmatrix_orig <- Fmatrix
    
    if (!isTRUE(control$skip_check_rank)){
      # Faster than inbuilt.
      # if (rankMatrix(X, method = 'qr') != ncol(X)){
      #   warning('Projected X not full Rank')
      # }
      rank_X <- rank_via_null(X, outer = TRUE)
      if (rank_X != ncol(X)){
        warning('Projected X not full rank')
      }
    }
    
    Fmatrix <- lapply(Fmatrix, FUN=function(F.j){
      lapply(F.j, FUN=function(l){
        t(basis_M) %*% l %*% basis_M
      })
    })
    Fmatrix <- lapply(Fmatrix, FUN=function(F.j){
      lapply(F.j, FUN=function(l){
        l <- drop0(zapsmall(l, clip_tiny))
        #l@x <- round(l@x, 10)
        return(l)
      })
    })
    
    #Get various ranks from \sum_\ell F_\ell and F_\ell
    #to use later in the algorithm.
    rank_F <- Reduce("+", lapply(Fmatrix, FUN=function(k){Reduce("+", k)}))
    rank_F <- rank_via_null(rank_F)
    
    if (rank_F != (ncol(X) - 1)){
      print(c(ncol(X), rank_F))
      warning('Rank F is unusual?')
      stop()
    }
    
    raw_F <- lapply(Fmatrix, FUN=function(F.j){
      raw_F.j <- lapply(F.j, FUN=function(ell){
        with(attributes(as(ell, 'dgTMatrix')), cbind(i, j, x))
      })
      size_raw.j <- sapply(raw_F.j, nrow)
      raw_F.j <- do.call('rbind', raw_F.j)
      return(list(raw = raw_F.j, size = size_raw.j))
    })
    raw_F <- list(
      raw = do.call('rbind', lapply(raw_F, FUN=function(i){i$raw})),
      size = do.call('c', lapply(raw_F, FUN=function(i){i$size}))
    )
    raw_F$dim <- dim(Fmatrix[[1]][[1]])

    toc(quiet = quiet_tictoc, log = TRUE)
  }else if (grepl(log_method, pattern='^log_')){
    
    # #Project into NULLSPACE for regression

    D <- do.call('rbind', lapply(Dlist, FUN=function(i){do.call('rbind', i)}))
    
    copy_main_D <- unlist(lapply(Dlist, FUN=function(i){sapply(i, nrow)}))
    need_to_copy <- which(copy_main_D != 1)
    copy_main_D <- c(cumsum(c(1, copy_main_D)))
    names(copy_main_D) <- NULL
    copy_main_D <- copy_main_D[-length(copy_main_D)]
    all_main_D <- copy_main_D
    copy_main_D <- copy_main_D[need_to_copy]
    
    if (length(copy_main_D) > 0){
      sparse_main_D <- sparseMatrix(i = copy_main_D, j = 1:length(copy_main_D), x = -1,
                                    dims = c(nrow(D), length(copy_main_D)))
      aug_restrictions <- rbind(
        cbind(D, -Diagonal(n = nrow(D)), sparse_main_D),
        cbind(constraint,
              drop0(sparseMatrix(i = 1, j = 1, x = 0,
                                 dims = c(nrow(constraint), length(copy_main_D) + nrow(D))))
        )
      )
    }else{
    
      aug_restrictions <- rbind(
        cbind(D, -Diagonal(n = nrow(D))),
        cbind(constraint,
              drop0(sparseMatrix(i = 1, j = 1, x = 0,
                                 dims = c(nrow(constraint), length(copy_main_D) + nrow(D))))
        )
      )
    }
    
    basis_aug <- calculate_nullspace_basis(aug_restrictions)
    basis_aug <- drop0(zapsmall(basis_aug, clip_tiny))

    p_orig <- ncol(X)
    
    if (log_method == 'log_ginv'){#Duplicate columns *after* projection
      
      M_ginv <- rbind(Diagonal(n = ncol(D)), D)
      M_ginv <- solve(crossprod(M_ginv), t(M_ginv))
      M_ginv <- cbind(M_ginv, M_ginv[,p_orig + copy_main_D])
      
      if (!do_BR){
        proj_X <- drop0(zapsmall(M_ginv %*% basis_aug, 15))
        X <- drop0(zapsmall(X %*% proj_X, clip_tiny))
      }else{
        X <- drop0(zapsmall(X %*% M_ginv, clip_tiny))
      }
      
      
    }else if (log_method == 'log_0'){#Z = [X, 0]
      
      X <- cbind(X, sparseMatrix(i=1,j=1,x=0, dims = c(nrow(X), length(copy_main_D) + nrow(D)))) 
      if (!do_BR){
        X <- drop0(zapsmall(X %*% basis_aug, clip_tiny))
      }
    }else if (log_method == 'log_random'){
      
      M <- rbind(Diagonal(n = ncol(D)), D)
      M_ginv <- solve(crossprod(M), t(M))
      part_X1 <- drop0(zapsmall(X %*% M_ginv, clip_tiny))
      proj_M <- Diagonal(n = nrow(M)) - M %*% M_ginv
      part_X2 <- rsparsematrix(nrow = nrow(X), ncol = ncol(proj_M), density = 0.1) %*% t(proj_M)
      part_X2[,1] <- 0
      part_X2 <- drop0(part_X2)
      X <- part_X1 + part_X2
      X <- cbind(X, X[ , p_orig + copy_main_D])
      
      if (!do_BR){
        X <- drop0(zapsmall(X %*% basis_aug))
      }
      
      rm(part_X1, part_X2)
      gc()
    }else{
      stop()
    }
    
    if (do_BR){
      v <- unlist(lapply(Dlist, FUN=function(i){sapply(i,nrow)}))
      names(v) <- NULL
      weight_duplicated <- apply(X[, c(p_orig + copy_main_D)], 
        MARGIN = 2, FUN=function(i){(sum(i^2))})
      #Get the *GROUPED* columns and calculate frobeniuns norm
      weight_group <- mapply(v, all_main_D, FUN=function(i,j){
        (sum( X[, p_orig + j:(i + j -1)]^2 ))
        # sqrt(sum(diag(crossprod( X[, i:(i + j -1)]))))
      })
      
      weight_group <- sqrt(weight_group)/sqrt(nrow(X))
      weight_duplicated <- sqrt(weight_duplicated)/sqrt(nrow(X))
      #Confirm that columns with dimension 1 are extracted correctly
      # stopifnot(length(setdiff(which(v == 1), which(weight_group[copy_main_D] == weight_duplicated))) == 0)
      
      X <- X %*% basis_aug
      X <- drop0(zapsmall(X, clip_tiny))
    }
    
    basis_M <- basis_aug
    p_X <- ncol(X)
    
    if (!all(X[,ncol(X)] == 1)){
      stop('Weird alignemnt issue with INTERCEPT not being last column after project. Examine X * basis_M')
    }
    
    tic('Proj F')
    
    Fmatrix_orig <- Fmatrix
    
    Fmatrix <- create_standard_group(Dlist, weight_dlist = weight_dlist)
    #Get the standard group penalties
    if (length(copy_main_D) != 0){
      extra_bdiag <- sparseMatrix(i=1,j=1,x=0,dims= rep(length(copy_main_D), 2))
    }else{
      extra_bdiag <- matrix(nrow = 0, ncol = 0)
    }
    Fmatrix <- lapply(Fmatrix, FUN=function(F.j){
      lapply(F.j, FUN=function(l){
        bdiag(sparseMatrix(i=1,j=1,x=0, dims = c(p_orig, p_orig)),
              l,
              extra_bdiag)
      })
    })
    
    if (length(copy_main_D) != 0){
      
      # *ADD* the group restrictions based on the NEW main effect only columns
      Fmatrix <- c(Fmatrix, list('log' = lapply(1:length(copy_main_D),
        FUN=function(i){
          bdiag(sparseMatrix(i=1,j=1,x=0,dims=rep(nrow(D) + p_orig,2)),
                sparseMatrix(i=i,j=i,x=1, dims = rep(length(copy_main_D),2)))
        })))
    }

    Fmatrix <- lapply(Fmatrix, FUN=function(F.j){
      lapply(F.j, FUN=function(l){
        t(basis_aug) %*% l %*% basis_aug
      })
    })
    
    Fmatrix <- lapply(Fmatrix, FUN=function(F.j){
      lapply(F.j, FUN=function(l){
        l <- drop0(zapsmall(l, clip_tiny))
        #l@x <- round(l@x, 10)
        return(l)
      })
    })
  
    #Get various ranks from \sum_\ell F_\ell and F_\ell
    #to use later in the algorithm.
    rank_F <- Reduce("+", lapply(Fmatrix, FUN=function(k){Reduce("+", k)}))
    rank_F <- rank_via_null(rank_F)
    
    if (rank_F != (ncol(X) - 1)){
      warning('Rank F is unusual?')
    }
    
    raw_F <- lapply(Fmatrix, FUN=function(F.j){
      raw_F.j <- lapply(F.j, FUN=function(ell){
        with(attributes(as(ell, 'dgTMatrix')), cbind(i, j, x))
      })
      size_raw.j <- sapply(raw_F.j, nrow)
      raw_F.j <- do.call('rbind', raw_F.j)
      return(list(raw = raw_F.j, size = size_raw.j))
    })
    raw_F <- list(
      raw = do.call('rbind', lapply(raw_F, FUN=function(i){i$raw})),
      size = do.call('c', lapply(raw_F, FUN=function(i){i$size}))
    )
    raw_F$dim <- dim(Fmatrix[[1]][[1]])
    
    toc(quiet = quiet_tictoc, log = TRUE)
  }else if (log_method == 'recast_0'){

    D <- do.call('rbind', lapply(Dlist, FUN=function(i){do.call('rbind', i)}))
    
    aug_restrictions <- rbind(
      cbind(D, -Diagonal(n = nrow(D))),
      cbind(constraint, 
            drop0(sparseMatrix(i = 1, j = 1, x = 0, 
                               dims = c(nrow(constraint), nrow(D))))
      )
    )
    
    aug_X <- cbind(X, sparseMatrix(i = 1, j =1, x = 0, dims = c(nrow(X), ncol = nrow(D))))
    basis_aug <- calculate_nullspace_basis(aug_restrictions)
    basis_aug <- drop0(zapsmall(basis_aug, clip_tiny))
    
    p_orig <- ncol(X)
      
    X <- aug_X %*% basis_aug
    X <- drop0(zapsmall(X, clip_tiny))
    
    basis_M <- basis_aug
    p_X <- ncol(X)
    
    if (!all(X[,ncol(X)] == 1)){
      stop('Weird alignemnt issue with INTERCEPT not being last column after project. Examine X * basis_M')
    }
    
    tic('Proj F')
    
    Fmatrix_orig <- Fmatrix
    
    if (!isTRUE(control$skip_check_rank)){
      # Faster than inbuilt.
      # if (rankMatrix(X, method = 'qr') != ncol(X)){
      #   warning('Projected X not flull Rank')
      # }
      rank_X <- rank_via_null(X, outer = TRUE)
      if (rank_X != ncol(X)){
        warning('Projected X not full rank')
      }
    }
    
    Fmatrix <- create_standard_group(Dlist, weight_dlist = weight_dlist)

    Fmatrix <- lapply(Fmatrix, FUN=function(F.j){
      lapply(F.j, FUN=function(l){
        bdiag(drop0(sparseMatrix(i=1,j=1,x=0, dims = c(p_orig, p_orig))), l)
      })
    })
    Fmatrix <- lapply(Fmatrix, FUN=function(F.j){
      lapply(F.j, FUN=function(l){
        t(basis_aug) %*% l %*% basis_aug
      })
    })
    
    Fmatrix <- lapply(Fmatrix, FUN=function(F.j){
      lapply(F.j, FUN=function(l){
        l <- drop0(zapsmall(l, clip_tiny))
        #l@x <- round(l@x, 10)
        return(l)
      })
    })
    
    #Get various ranks from \sum_\ell F_\ell and F_\ell
    #to use later in the algorithm.
    rank_F <- Reduce("+", lapply(Fmatrix, FUN=function(k){Reduce("+", k)}))
    rank_F <- rank_via_null(rank_F)
    
    if (rank_F != (ncol(X) - 1)){
      warning('Rank F is unusual?')
    }
    
    raw_F <- lapply(Fmatrix, FUN=function(F.j){
      raw_F.j <- lapply(F.j, FUN=function(ell){
        with(attributes(as(ell, 'dgTMatrix')), cbind(i, j, x))
      })
      size_raw.j <- sapply(raw_F.j, nrow)
      raw_F.j <- do.call('rbind', raw_F.j)
      return(list(raw = raw_F.j, size = size_raw.j))
    })
    raw_F <- list(
      raw = do.call('rbind', lapply(raw_F, FUN=function(i){i$raw})),
      size = do.call('c', lapply(raw_F, FUN=function(i){i$size}))
    )
    
    toc(quiet = quiet_tictoc, log = TRUE)
    
  }
  
  if (lambda == 0){
    if (single_intercept){
      rank_F <- (ncol(X) - 1)/2
    }else{
      rank_F <- ncol(X)/2
    }
  }
  #Number of levels for each factor
  n_levels <- sapply(factor_levels, length)
  tic('Build Adaptive')

  if (do_BR){
    
    if (log_method == 'standard'){
      br_weights <- standardization_weights(D_list = Dlist, X = original_X)
    }else if (grepl(log_method, pattern='^log')){
      br_weights <- split(c(weight_group, weight_duplicated), rep(names(Fmatrix), lengths(Fmatrix)))
      br_weights <- br_weights[names(Fmatrix)]
      
      if (log_method == 'log_0'){
        message('B&R Weights Ignored if log_0')
        br_weights <- lapply(br_weights, FUN=function(i){i[] <- 1; return(i)})  
      }
      if (any(unlist(br_weights) == 0)){
        br_weights <- lapply(br_weights, FUN=function(i){
          i[which(i == 0)] <- 1
          return(i)
        })
      }
    }else{stop('Invalid "log" method.')}
    br_weights <- rep(list(br_weights), K)
  }else{
    
    weight_func <- function(L){
      ifelse(L == 1, 1, sqrt(L) * (L + 1))
    }
    
    if ('log' %in% names(Fmatrix) & !('log' %in% names(n_levels))){
      n_levels <- c(n_levels, 'log' = 1)
    }
    
    br_weights <- mapply(Fmatrix, n_levels, SIMPLIFY = FALSE, FUN=function(F.j, L.j){
      weight_Lj <- weight_func(L.j)
      sapply(F.j, FUN=function(ell){
        adapt_wj <- weight_Lj
        return(1/adapt_wj)
      })
    })
    
    br_weights <- rep(list(br_weights), K)
  }

  if (!inherits(adaptive_weight, 'character')){
    
    if (inherits(adaptive_weight, c('matrix')) | inherits(adaptive_weight, 'dgeMatrix')){
      
      if (!(ncol(adaptive_weight) %in% c(1, K))){
        stop('If adaptive_weight is a matrix, must have either 1 or K columns.')
      }
      names(coef_names) <- NULL
      r_aw <- rownames(adaptive_weight)
      names(r_aw) <- NULL
      if (!identical(r_aw, coef_names)){
        print(coef_names)
        stop('Row names on adaptive_weight must match coef_names')
      }

      clip_weight <- 10^5 #1/sqrt(.Machine$double.eps)
      

      adaptive_weight <- lapply(1:K, FUN=function(k){
        if (ncol(adaptive_weight) == K){
          b <- adaptive_weight[,k]
        }else{
          b <- adaptive_weight[,1]
        }
        a_br <- br_weights[[k]]
        if (grepl(log_method, pattern='^log')){
          
          b <- c(F_EM_update(b, Fmat = Fmatrix_orig, lambda = 1, 
                      pi.gamma = 1, exo_weights = NULL),
          list('log' = cbind(1/abs(as.vector(D[copy_main_D,] %*% b)), NA)))
          if (length(copy_main_D) == 0){
            b <- b[-length(b)]
          }
        }else{
          b <- F_EM_update(b, Fmat = Fmatrix_orig, lambda = 1, 
                           pi.gamma = 1, exo_weights = NULL)
        }
        b <- mapply(b, a_br, SIMPLIFY = FALSE, FUN=function(i, w){
          i <- i[,1] * w
          i[i > clip_weight] <- clip_weight
          return(i)
        })
        
      })
      
      if (length(adaptive_weight) == 1){
        adaptive_weight <- rep(adaptive_weight, K)
      }
    }else{
      stop('Must provide coefficients from FactorHet object or matrix with names.')
    }

    if ('log' %in% names(Fmatrix) & !('log' %in% names(n_levels))){
      n_levels <- c(n_levels, 'log' = 1)
    }
    
    adaptive_weight <- mapply(Fmatrix, n_levels, SIMPLIFY = FALSE, FUN=function(F.j, L.j){
      weight_Lj <- weight_func(L.j)
      sapply(F.j, FUN=function(ell){
        adapt_wj <- weight_Lj
        return(1/adapt_wj)
      })
    })
    adaptive_weight <- rep(list(adaptive_weight), K)
  }else if (adaptive_weight == 'none'){
    adaptive_weight <- mapply(Fmatrix, SIMPLIFY = FALSE, FUN=function(F.j){
      sapply(F.j, FUN=function(ell){
        return(1)
      })
    })
    adaptive_weight <- rep(list(adaptive_weight), K)
  }else if (adaptive_weight == 'B&R'){
    adaptive_weight <- br_weights
  }else{stop('')}

  toc(quiet = quiet_tictoc, log = TRUE)
  ##Group Level z_{i,k}

  if (length(group) != nrow(X)){stop('Group is wrong length')}


  if (tau_method == 'nullspace'){
    binding_null_basis <- NULL
    # if (log_method == 'standard'){
    #   binding_lookup <- build_binding_lookup(Fmatrix = Fmatrix, factor_levels = factor_levels, term_position = term_position, coef_names = coef_names)
    # }
    existing_restrictions <- lapply(1:K, FUN=function(i){
      i <- lapply(1:length(Fmatrix), FUN=function(j){c()})
      names(i) <- names(Fmatrix)
      return(i)
    })
    tracking_restrictions <- matrix(NA, nrow = iterations, ncol = 2)
  }else if (tau_method == 'clip'){
    binding_null_basis <- NULL
    tracking_restrictions <- NULL
  }else{stop('tau_method must be clip or nullspace.')}
  
  #Initialize Parameters  
  tic('Make Init')
  update_mod <- list()
  phi <- matrix(0, nrow = K, ncol = ncol(W))
  
  if (K == 1){
    group_E.prob <- matrix(1, nrow = n_G)
    obs.E.prob <- matrix(1, nrow = nrow(X))
    beta <- simple_logit(y = y, X = X, iterations = 10, weights = weights,
      beta_method = 'cg', beta_cg_it = 10, prec_ridge = 1)
    beta <- matrix(as.vector(beta))  
  }else if (inherits(control$init_method, "matrix")){
    stop('init_method must be list or character.')
  }else if (inherits(control$init_method, "list")){
    valid_init <- all(c('pi', 'beta', 'group_E.prob', 'phi') %in% names(control$init_method))
    if (!valid_init){
      
      if (!identical('group_E.prob', names(control$init_method))){
        stop('A list-based initialization must contain either (i) a single data.frame of probabilities for each group/unit with the column names "group" and then "cluster_[0-9]}" or "[0-9]", (ii) a full specification of initialization (uncommon)')
      }
      if (length(setdiff(unique_group, control$init_method$group_E.prob$group)) > 0){
        stop('Groups missing from init_method that are in design.')
      }
      if (length(setdiff(control$init_method$group_E.prob$group, unique_group)) > 0){
        stop('Groups found in init_method that are not in design.')
      }
      group_E.prob <- control$init_method$group_E.prob[match(unique_group, control$init_method$group_E.prob$group),]
      if (inherits(group_E.prob$group, "factor")){
        group_E.prob$group <- as.character(group_E.prob$group)
      }
      stopifnot(identical(group_E.prob$group, unique_group))
      group_E.prob <- dplyr::select(group_E.prob, -group)
      group_E.prob <- as.matrix(group_E.prob)
      rownames(group_E.prob) <- NULL
      
      if (ncol(group_E.prob) != K){stop('init_method as list must contain 1 + K columns (group + probabilities)')}
      
      obs.E.prob <- apply(group_E.prob, MARGIN = 2, FUN=function(i){as.vector(group_mapping %*% i)})
      
      if (all(obs.E.prob %in% c(0,1))){
        beta <- sapply(list_from_cols(obs.E.prob), FUN=function(k){
          simple_logit(y = y[which(k == 1)], weights = weights[which(k == 1)],
            X = X[which(k == 1),,drop=F], iterations = 10, 
            beta_method = 'cg', beta_cg_it = 10, prec_ridge = 1)
        })
      }else{
        beta <- sapply(list_from_cols(obs.E.prob), FUN=function(k){
          simple_logit(y = y, X = X, obs.E.prob = k, weights = weights,
                       iterations = 10, beta_method = 'cg', beta_cg_it = 10, prec_ridge = 1)
        })
      }
      rownames(beta) <- NULL
      
    }else{#FULL initialization; usually short EM internal usage.
      pi <- control$init_method$pi
      phi <- control$init_method$phi
      beta <- control$init_method$beta
      group_E.prob <- control$init_method$group_E.prob
    }
    
  }else if (control$init_method %in% c('mclust')){
    
    if (requireNamespace('mclust', quietly = TRUE)){
      hcVVV <- mclust::hcVVV
      mclust_assign <- as.vector(mclust::hclass(mclust::hc(W, use = "SVD"), G = K))
      group_E.prob <- as.matrix(sparseMatrix(i = 1:nrow(W), j = mclust_assign, x = 1))
      obs.E.prob <- apply(group_E.prob, MARGIN = 2, FUN=function(i){as.vector(group_mapping %*% i)})
    }else{
      stop("'mclust' must be installed.")
    }
    
    
    beta <- sapply(list_from_cols(obs.E.prob), FUN=function(k){
      simple_logit(y = y[which(k == 1)], X = X[which(k == 1),], iterations = 10, 
          weights = weights[which(k == 1)],
          beta_method = 'cg', beta_cg_it = 10, prec_ridge = 1)
      #coef(bayesglm(y ~ 0 + as.matrix(X), weights = k[match(group, unique_group)], family = binomial))
    })
    rownames(beta) <- NULL
    
  }else if (control$init_method == 'kmeans'){
    #Do K-means and hard assignment   
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
    
    if (ncol(scale_W) == 1 & all(scale_W[,1] == 0)){
      warning('With no moderators, k-means replaced by random assignment; it is best to respecify to ensure short_EM is done as expected')
      group_E.prob <- sample(1:K, nrow(W), replace = T)
      group_E.prob <- sparseMatrix(i = 1:nrow(W), j = group_E.prob, x = 1)
    }else{
      init_k <- suppressWarnings(kmeans(x = scale_W, centers = K, iter.max = 25, nstart = 1500)$cluster)
      group_E.prob <- as.matrix(sparseMatrix(i = 1:nrow(W), j = init_k, x = 1))
    }
    obs.E.prob <- apply(group_E.prob, MARGIN = 2, FUN=function(i){as.vector(group_mapping %*% i)})
    
    beta <- sapply(list_from_cols(obs.E.prob), FUN=function(k){
      simple_logit(y = y[which(k == 1)], X = X[which(k == 1),], weights = weights[which(k == 1)],
                   iterations = 10, beta_method = 'cg', beta_cg_it = 10, prec_ridge = 1)
      #coef(bayesglm(y ~ 0 + as.matrix(X), weights = k[match(group, unique_group)], family = binomial))
    })
    rownames(beta) <- NULL
  }else if (control$init_method == 'random_pi'){
    pi <- rexp(K)
    pi <- pi/sum(pi)
    
    simK <- sample(1:K, n_G, prob = pi, replace = T)
    
    if (length(unique(simK)) != K){
      warning('Unbalanced Initialization for Random Pi: Doing with balanced 1/K')
      pi <- rep(1/K, K)
      simK <- sample(1:K, n_G, prob = pi, replace = T)
    }
    group_E.prob <- as.matrix(sparseMatrix(i = 1:n_G, j = simK, x = 1, dims = c(n_G, K)))
    obs.E.prob <- apply(group_E.prob, MARGIN = 2, FUN=function(i){as.vector(group_mapping %*% i)})
    
    beta <- sapply(list_from_cols(obs.E.prob), FUN=function(k){
      simple_logit(y = y[which(k == 1)], X = X[which(k == 1),], iterations = 10, 
                   weights = weights[which(k == 1)],
                   beta_method = 'cg', beta_cg_it = 10, prec_ridge = 1)
      #coef(bayesglm(y ~ 0 + as.matrix(X), weights = k[match(group, unique_group)], family = binomial))
    })
    rownames(beta) <- NULL
  }else if (control$init_method == 'random_member'){
    simK <- sample(1:K, n_G, replace = T)
    if (length(unique(simK)) != K){
      warning('Degenerat random allocation; initializing each probabilistically')
      group_E.prob <- t(sapply(1:n_G, FUN=function(i){i <- rexp(K); return(i/sum(i))}))
      obs.E.prob <- apply(group_E.prob, MARGIN = 2, FUN=function(i){as.vector(group_mapping %*% i)})
      beta <- sapply(list_from_cols(obs.E.prob), FUN=function(k){
        simple_logit(y = y, X = X, 
                     weights = weights,
                     obs.E.prob = k, iterations = 10, beta_method = 'cg', beta_cg_it = 10, prec_ridge = 1)
        #coef(bayesglm(y ~ 0 + as.matrix(X), weights = k[match(group, unique_group)], family = binomial))
      })
      
    }else{
      group_E.prob <- as.matrix(sparseMatrix(i = 1:n_G, j = simK, x = 1, dims = c(n_G, K)))
      obs.E.prob <- apply(group_E.prob, MARGIN = 2, FUN=function(i){as.vector(group_mapping %*% i)})
      beta <- sapply(list_from_cols(obs.E.prob), FUN=function(k){
        simple_logit(y = y[which(k == 1)], X = X[which(k == 1),], iterations = 10, 
                  weights = weights[which(k == 1)],
                  beta_method = 'cg', beta_cg_it = 10, prec_ridge = 1)
      })
    }
    rownames(beta) <- NULL
  }else if (control$init_method == 'random_beta'){
    #Random initialization
    beta <- matrix(rnorm(ncol(X) * K), ncol = K)
    #Generate cluster membership probabilities randomly.
    group_E.prob <- matrix(pi, ncol = K, nrow = n_G, byrow=T)
  }else{stop('Invalid initialization method.')}

  # If no provided phi, use initializations to get estimate.
  do_init_phi <- tryCatch(is.null(control$init_method$phi), error = function(e){TRUE})
  if (do_init_phi){
    if (K == 1){
      pi <- 1
    }else{
      #If initialization not provided, do update to calibrate phi
      update_phi <- update_phi(phi = phi, W = W,
           ridge_phi = ridge_phi, weights_W = weights_W,
           group_E.prob = group_E.prob, K = K,
           maxit_phi = 10)
      pi <- update_phi$pi_bar
      phi <- update_phi$phi
    }
  }
  
  toc(quiet = quiet_tictoc, log = T)
  
  #Set up placeholders for tracking progress of algorithm.
  old.beta <- beta
  old.beta[,] <- -Inf
  running_ll <- rep(-Inf, 3)
  store_ll <- matrix(nrow = iterations, ncol = 3)
  store_beta <- array(NA, dim = c(iterations, nrow(beta), K))
  
  toc(log = TRUE, quiet = TRUE)

  progress <- floor(iterations / 20)
  progress <- max(c(1, progress))

  if (do_SQUAREM){
    SQUAREM_counter <- 0
    SQUAREM_list <- list()
    SQUAREM_backtrack <- control$backtrack_SQUAREM
    SQUAREM_track <- matrix(NA, nrow = iterations, ncol = 2)
  }else{
    SQUAREM_track <- NULL
  }
  
  for (it in 1:iterations){
    
    if (it %% progress == 0){message('.', appendLF = FALSE)}
    
    for (m in 1:repeat_beta){#Repeat multiple times.
      tic('Main E Step')
      ##############################
      ############E-Step############
      ##############################
      xb <- as.matrix(X %*% beta) 
      
      #Get log[L(y_i | \beta_r) * pi_r] = [\sum_{t} ll(y_{it} | \beta_r)] + log(pi_r)
      loglik.k <- plogis(xb * (2 * y - 1), log.p = TRUE) 
      
      group_E.prob <- calculate_posterior_zi(
        loglik.k = loglik.k, group_mapping = group_mapping, K = K, ncol_W = ncol_W, 
        pi = pi, phi = phi, W = W
      ) 
      
      #Reconstitute E[z_{i}] for each observation (i,t) for weighting
      #regression
      obs.E.prob <- apply(group_E.prob, MARGIN = 2, FUN=function(i){as.vector(group_mapping %*% i)})
      
      
      #Get E[1/tau^2_k | -] for weight
      
      if (lambda > 0){
        tic('Calculate Tau')
        
        E.tau <- mapply(list_from_cols(beta), pi^gamma, adaptive_weight, SIMPLIFY = FALSE, FUN=function(b, pi.g, aw){
          F_EM_update(beta = b, Fmat = Fmatrix, lambda = lambda, pi.gamma = pi.g, exo_weights = aw)
        })
        
        #Clip E[1/tau^2_k] at big # (1e10) to prevent numerical instability
        #When tau_method = 'nullspace', only relevant for initial stabilization steps
        #as tau_equality_threshold governs when two levels are fused.
        
        E.tau <- clip_Etau(E.tau, threshold = tau_truncate)
        toc(quiet = quiet_tictoc, log = TRUE)
      }

      #Get E[\omega_{ik} | -] for polya-gamma
      E.omega <- obs.E.prob/(2 * xb) * tanh(xb/2)
      
      if (any(abs(xb) < 1e-10)){
        #Deal with possible numerical instability for x_i^T\beta \approx 0
        #Note that lim_{x \to 0} 1/(2 * x) tanh(x/2) = 1/4
        E.omega[which(abs(xb) < 1e-10)] <- obs.E.prob[which(abs(xb) < 1e-10)] * 1/4
      }
      
      E.omega <- Diagonal(x = weights) %*% E.omega
      
      #Get E[\sum_{l} F_l / tau^2_{l,r}] = \sum_l F_l * E[1/tau^2_{l,r}]
      tic('Prepare Ridge')
      if (lambda == 0){#if lambda = 0, add slight ridge stabilization
        #Create a ridge prior on the UNTRANSFORMED space and then project
        #into the nullspace
        if (it == 1){
          # -1/2 \ln(2 * pi * sigma^2_b) - 1/(2 sigma^2_b) b_j^2
          #
          ridge_beta <- c(0, rep(ridge_beta, nrow(basis_M) - 1))

          ridge_beta <- t(basis_M) %*% sparse_diag(ridge_beta) %*% basis_M
          ridge_beta <- drop0(zapsmall(ridge_beta, clip_tiny))
          ridge_beta <- as(ridge_beta, 'dgCMatrix')
          raw.E.ridge <- lapply(1:K, FUN=function(i){ridge_beta})
          
        }
        
        E.ridge <- mapply(raw.E.ridge, pi, SIMPLIFY = FALSE, FUN=function(r, p){
          r * p^gamma
        })

      }else{
        E.ridge <- list()
        for (k in 1:K){#For each cluster
            #Sum over factors j \in 1....J
            E.ridge[[k]] <- make_ridge(E.tauk = E.tau[[k]], raw_F = raw_F, method = 'raw')
            
            # ##Checksum to make sure "raw method works"
            # ridge_2 <- Reduce('+', mapply(E.tau[[k]], Fmatrix, SIMPLIFY = FALSE, FUN=function(etau, F.j){
            #   #Sum over the L (l) F matrices for each factor.
            #   ridge.j <- Reduce('+', mapply(etau[,1], F.j, SIMPLIFY = FALSE, FUN=function(etau.l, l){
            #     etau.l * l
            #   }))
            # }))
            # stopifnot(all.equal(ridge_2, E.ridge[[k]]))
        }
      }
      toc(quiet = quiet_tictoc, log = TRUE)
      
      
      toc(quiet = quiet_tictoc, log = TRUE)
      
      if (check_all_ll){
        tic('Check LL') #Compare LL vs. prior iteration. Should increase!
        new_ll <- evaluate_loglik(beta = beta, pi = pi, X = X, y = y, group_mapping = group_mapping,
                                  gamma = gamma, E.prob = group_E.prob, adaptive_weight = adaptive_weight,
                                  phi = phi, W = W, ridge_phi = ridge_phi, ridge_beta = ridge_beta,
                                  weights_W = weights_W,
                                  Fmatrix = Fmatrix, lambda = lambda, rank_F = rank_F, separate = TRUE)
        running_ll <- compare_ll(new_ll, running_ll, stage = 'E', debug_ll = debug)
        toc(quiet = quiet_tictoc, log = TRUE)
      }
      #####M-Step#####
      ###Update Beta
      #Temporary to freeze E[1/tau^2_{l,r}] to check performance
      #E.ridge <- lapply(E.ridge, FUN=function(i){Diagonal(x = rep(1/25, ncol(i)))})
      tic('Update Beta')
      old_beta <- beta
      if (it == 1){#Create large blocked X for each later iteration
        if (single_intercept){
          blocked_X <- cbind(1, kronecker(diag(K), X[,-ncol(X)]))
        }else{
          blocked_X <- kronecker(diag(K), X)
        }
      }
      
      if (lambda == 0 | tau_method == 'clip' | it < tau_stabilization){
        if (debug){message('.', appendLF = F)}
        if (single_intercept){
          
          beta <- update_beta(X=NULL, blocked_X = blocked_X, p_X = p_X, prior_beta = beta,
                              weights = weights,
                              y=y,E.omega = E.omega, obs.E.prob = obs.E.prob, method = beta_method,
                              E.ridge = E.ridge, K=K, global_int = single_intercept, cg_it = beta_cg_it)
        }else{
          
          beta <- update_beta(X=X,y=y,E.omega = E.omega, obs.E.prob = obs.E.prob,
              p_X = p_X, method = beta_method, weights = weights,
              E.ridge = E.ridge, K=K, global_int = single_intercept, prior_beta = beta)
        }
        
      }else if (tau_method == 'nullspace'){
        
        # # Get L2 distance
        # binding_restrictions <- apply(beta, MARGIN = 2, 
        #   FUN=function(b){F_EM_update(b, Fmatrix, lambda = 1, pi.gamma = 1, exo_weights = 1)})
        # binding_restrictions <- lapply(binding_restrictions, FUN=function(i){
        #   lapply(i, FUN=function(j){which(j[,1] > 1/tau_equality_threshold)})
        # })
        
        # if (log_method != 'standard'){
        binding_restrictions <- lapply(E.tau, FUN=function(etau.k){
          lapply(etau.k, FUN=function(k){which(k[,1] >= tau_truncate)})
        })
        # }else{
        #   binding_restrictions <- apply(basis_M %*% beta, MARGIN = 2, FUN=function(b){
        #     br <- calc_distance_binding(b, binding_lookup)
        #     # br <- lapply(prepare_fusion(factor_levels = factor_levels, term_position = term_position, coef_names = coef_names, beta = b, simplify = FALSE), FUN=function(i){
        #     #   sapply(i, FUN=function(j){max(unlist(j))})
        #     # })
        #     br <- lapply(br, FUN=function(j){which(j <= tau_equality_threshold)})
        #     return(br)
        #   })
        # }

        if (length(unlist(binding_restrictions)) == 0){#If nothing binds, do standard "clip".
          if (debug){message('|', appendLF = F)}
          if (single_intercept){
            tic('cg_est')
            beta <- update_beta(X=NULL, blocked_X = blocked_X, p_X = p_X, method = beta_method, cg_it = beta_cg_it,
                                y=y,E.omega = E.omega, obs.E.prob = obs.E.prob, prior_beta = beta,
                                weights = weights,
                                E.ridge = E.ridge, K=K, global_int = single_intercept)
            toc(quiet = quiet_tictoc, log = T)
          }else{
            beta <- update_beta(X=X,y=y,E.omega = E.omega, obs.E.prob = obs.E.prob, p_X = p_X, method = beta_method,
                                prior_beta = beta, cg_it = beta_cg_it,
                                weights = weights,
                                E.ridge = E.ridge, K=K, global_int = single_intercept)
          }
          
        }else{
          
          if (do_SQUAREM){
            #Verify that no restrictions have been *REMOVED* if so, reset nullspace.
            removed_restrictions <- mapply(existing_restrictions, binding_restrictions, SIMPLIFY = F, FUN=function(e,b){
              mapply(e,b, FUN=function(e_j, b_j){setdiff(e_j,b_j)}, SIMPLIFY = F)
            })
            if (length(unlist(removed_restrictions) > 0)){
              existing_restrictions <- binding_restrictions
              reset_nullspace <- TRUE
            }else{
              reset_nullspace <- FALSE
            }
          }else{reset_nullspace <- FALSE}
          
          all_restrictions <- mapply(existing_restrictions, binding_restrictions, SIMPLIFY = F, FUN=function(e,b){
            mapply(e,b, FUN=function(e_j, b_j){union(e_j,b_j)}, SIMPLIFY = F)
          })
          
          new_restrictions <- mapply(all_restrictions, existing_restrictions, SIMPLIFY = F, FUN=function(e,b){
            mapply(e,b, FUN=function(e_j, b_j){setdiff(e_j,b_j)}, SIMPLIFY = F)
          })
          
          if (debug){
            if (reset_nullspace){
              message('r', appendLF = F)
              message(length(unlist(all_restrictions)), appendLF = F)
            }else{
              message(length(unlist(all_restrictions)), appendLF = F)

            }
            message('-', appendLF = F)
          }
          
          # Checksum for basis
          # b2 <- drop0(new_null_basis %*% joint_basis[-1:-ncol(binding_null_basis),])
          #
          # if (!all.equal(b1, b2)){
          #   b1 <- b1
          #   b2 <- b2
          #   stop('Null Spaces do not match ')
          # }
          
          if (single_intercept){
            if (is.null(binding_null_basis) | reset_nullspace){
              binding_null_basis <- lapply(all_restrictions, FUN=function(r_k){
                binding_k <- do.call('rbind', mapply(r_k, Fmatrix, SIMPLIFY = FALSE, FUN=function(i, F_j){
                  do.call('rbind', F_j[i])
                }))
                binding_k <- binding_k[,-p_X]
                if (is.null(binding_k)){
                  return(sparse_diag(rep(1, p_X - 1)))
                }else{
                  drop0(zapsmall(calculate_nullspace_basis(binding_k), clip_tiny))
                }
              })
              binding_null_basis <- bdiag(c(1, bdiag(binding_null_basis)))
            }else if (length(unlist(new_restrictions)) > 0){
              new_null_basis <- lapply(new_restrictions, FUN=function(r_k){
                binding_k <- do.call('rbind', mapply(r_k, Fmatrix, SIMPLIFY = FALSE, FUN=function(i, F_j){
                  do.call('rbind', F_j[i])
                }))
                binding_k <- binding_k[,-ncol(X)]
                if (is.null(binding_k)){
                  sparse_diag(rep(1, ncol(X) - 1))
                }else{
                  drop0(zapsmall(calculate_nullspace_basis(binding_k), clip_tiny))
                }
              })
              new_null_basis <- bdiag(c(1, new_null_basis))
              
              joint_basis <- calculate_nullspace_basis(as(crossprod(cbind(binding_null_basis, - new_null_basis)), 'dgCMatrix'))
              b1 <- drop0(binding_null_basis %*% joint_basis[1:ncol(binding_null_basis),])
              binding_null_basis <- drop0(zapsmall(b1, clip_tiny))
              rm(b1)
            }
          }else{
            if (is.null(binding_null_basis) | reset_nullspace){
              binding_null_basis <- lapply(all_restrictions, FUN=function(r_k){
                binding_k <- do.call('rbind', mapply(r_k, Fmatrix, SIMPLIFY = FALSE, FUN=function(i, F_j){
                  do.call('rbind', F_j[i])
                }))
                if (is.null(binding_k)){
                  sparse_diag(rep(1, ncol(X)))
                }else{
                  drop0(zapsmall(calculate_nullspace_basis(binding_k), clip_tiny))
                }
              })
              binding_null_basis <- bdiag(binding_null_basis)
            }else if (length(unlist(new_restrictions)) > 0){
              new_null_basis <- lapply(new_restrictions, FUN=function(r_k){
                binding_k <- do.call('rbind', mapply(r_k, Fmatrix, SIMPLIFY = FALSE, FUN=function(i, F_j){
                  do.call('rbind', F_j[i])
                }))
                if (is.null(binding_k)){
                  sparse_diag(rep(1, ncol(X)))
                }else{
                  drop0(zapsmall(calculate_nullspace_basis(binding_k), clip_tiny))
                }
              })
              new_null_basis <- bdiag(new_null_basis)
              
              joint_basis <- calculate_nullspace_basis(as(crossprod(cbind(binding_null_basis, - new_null_basis)), 'dgCMatrix'))
              b1 <- drop0(binding_null_basis %*% joint_basis[1:ncol(binding_null_basis),])
              binding_null_basis <- drop0(zapsmall(b1, clip_tiny))
              rm(b1)
            }
          }
          
          existing_restrictions <- all_restrictions
          tracking_restrictions[it,] <- c(length(unlist(new_restrictions)), ncol(binding_null_basis))
        
          if (single_intercept){
            blocked_E <- bdiag(c(0, E.ridge))
            blocked_E <- blocked_E[-(1 + ncol(X) * 1:K),-(1 + ncol(X) * 1:K)]
          }else{
            blocked_E <- bdiag(E.ridge)
          }
            
          if (beta_method == 'cg'){

            tic('cg_init')
            if (single_intercept){
              cg_flat_beta <- c(beta[p_X,1], as.vector(beta[-p_X,]))
              cg_flat_null <- matrix(as.vector(solve(crossprod(binding_null_basis), t(binding_null_basis) %*% cg_flat_beta)))
            }else{
              cg_flat_beta <- as.vector(beta)
              cg_flat_null <- matrix(as.vector(solve(crossprod(binding_null_basis), t(binding_null_basis) %*% cg_flat_beta)))
            }
            toc(quiet = quiet_tictoc, log = T)
            
            tic('cg_est')
            
            cg_b <- cg_custom(K = 1,X = blocked_X %*% binding_null_basis,
                  s = rep(weights * (y - 1/2), K) * as.vector(obs.E.prob),
                  omega = matrix(as.vector(E.omega)),
                  list_ridge = list(t(binding_null_basis) %*% blocked_E %*% binding_null_basis), 
                  tol = sqrt(.Machine$double.eps), weights = matrix(0, nrow = 1, ncol = 0),
                  it_max = beta_cg_it, old_beta = cg_flat_null)

            blocked_beta_null <- cg_b$beta
              
            toc(quiet = quiet_tictoc, log = T)

          }else if (beta_method == 'cpp'){
            
            # if (!single_intercept){stop()}
            
            blocked_beta_null <- cpp_beta_plain(X = blocked_X %*% binding_null_basis,
              s = rep(weights * (y - 1/2), K) * as.vector(obs.E.prob),
              K = K, omega = sparse_diag(as.vector(E.omega)),
              ridge = t(binding_null_basis) %*% blocked_E %*% binding_null_basis
            )
            
          }else{stop('invalid beta method: cpp or cg')}
          
          if (single_intercept){
            blocked_beta <- as.vector(binding_null_basis %*% blocked_beta_null)
            mu <- blocked_beta[1]
            
            beta <- matrix(blocked_beta[-1], ncol = K)
            beta <- rbind(beta, mu)
          }else{
            blocked_beta <- as.vector(binding_null_basis %*% blocked_beta_null)
            beta <- matrix(blocked_beta, ncol = K)
          }

        }
      }
      
      
    }
    
    store_beta[it,,] <- beta
    toc(quiet = quiet_tictoc, log = TRUE)
    
    if (check_all_ll){
      tic('Check LL') #Check progress of LL
      new_ll <- evaluate_loglik(beta = beta, pi = pi, X = X, y = y, group_mapping = group_mapping,
                                gamma = gamma, E.prob = group_E.prob, phi = phi,adaptive_weight = adaptive_weight,
                                W = W, ridge_phi = ridge_phi, ridge_beta = ridge_beta, weights_W = weights_W,
                                Fmatrix = Fmatrix, lambda = lambda, rank_F = rank_F, separate = TRUE)
      running_ll <- compare_ll(new_ll, running_ll, stage = 'beta', debug_ll = debug)    
      toc(quiet = quiet_tictoc, log = TRUE)
    }
    if (any(is.na(beta))){
      stop('NaN found in beta update.')
    }
    change.beta <- apply(abs(beta - old.beta), MARGIN = 2, max)

    ###Update Pi
    old.phi <- phi
    tic('Update Pi')
    if (K == 1){
      pi <- 1
    }else{
      #If doing gamma > 0, then must do another
      #E Step for E[z_ir | \beta, pi] to maintain EM algorithm (AECM)
      xb <- as.matrix(X %*% beta) 
      
      #Get log[L(y_i | \beta_r) * pi_r] = [\sum_{t} ll(y_{it} | \beta_r)] + log(pi_r)
      loglik.k <- plogis(xb * (2 * y - 1), log.p = TRUE) 
      
      group_E.prob <- calculate_posterior_zi(
        loglik.k = loglik.k, group_mapping = group_mapping, K = K, ncol_W = ncol_W, 
        pi = pi, phi = phi, W = W
      ) 

      if (K == 1){
        #Do nothing
      }else if (gamma == 0){
      #If gamma = 0, parametric independence so prior on \phi depends only on 
      # moderators
        update_phi <- update_phi(phi = phi, W = W, 
                          ridge_phi = ridge_phi, weights_W = weights_W,
                          group_E.prob = group_E.prob, K = K,
                          maxit_phi = maxit_pi)
        pi <- update_phi$pi_bar
        phi <- update_phi$phi
        
      }else{
      #If gamma = 1, parametric dependence so prior on \phi
      # depends on beta and pi so more complicated to update
        update_mod <- update_moderator(K = K, phi = phi, beta = beta,
         W = W, weights_W = weights_W,
         maxit_mod = maxit_pi, adaptive_weight = adaptive_weight,
         rank_F = rank_F, Fmatrix = Fmatrix, group_E.prob = group_E.prob,
         ridge_beta = ridge_beta, single_intercept = single_intercept,
         ridge_phi = ridge_phi, gamma = gamma, lambda = lambda,
         extra_opt_args = optim_phi_ctrl)

        phi <- update_mod$phi
        pi <- update_mod$pi_bar
        
      }
      colnames(phi) <- colnames(W)
      rownames(phi) <- 1:K
    }
    if (ncol_W > 0 & K > 1){
      change.phi <- max(abs(phi - old.phi))
    }else{change.phi <- 0}
    toc(quiet = quiet_tictoc, log = TRUE)

    tic('Check LL') #Check progress of loglik
    new_ll <- evaluate_loglik(beta = beta, pi = pi, X = X, y = y, group_mapping = group_mapping,
      E.prob =  group_E.prob, phi = phi, W = W, weights_W = weights_W,
      ridge_phi = ridge_phi, ridge_beta = ridge_beta,adaptive_weight = adaptive_weight,
      gamma = gamma, Fmatrix = Fmatrix, lambda = lambda, 
      rank_F = rank_F, separate = TRUE)
    toc(quiet = quiet_tictoc, log = TRUE)
    
    if (do_SQUAREM){
      
      temp_ll <- compare_ll(new_ll, running_ll, stage = 'pi', debug_ll = debug)
      tic('SQUAREM')
      SQUAREM_counter <- SQUAREM_counter + 1
      SQUAREM_list[[SQUAREM_counter]] <- list(beta = beta, phi = phi)
      
      
      if (it %% 3 == 0){
        SQUAREM_counter <- 0
        SQUAREM_list <- SQUAREM_list
        
        #Prepare SQUAREM update
        proposed_SQUAREM <- prepare_SQUAREM(SQUAREM_list, scale_alpha = 2)
        
        if (proposed_SQUAREM$alpha != -1){
          
          if (K == 1){
            proposed_pi <- 1 #mean(proposed_postpred)
          }else{
            proposed.loglik.k <- plogis(as.matrix(X %*% proposed_SQUAREM$update$beta) * 
                (2 * y - 1), log.p = TRUE) 
            
            proposed_pi <- colSums((Diagonal(x = as.vector(weights_W)/sum(weights_W)) %*% 
                softmax_matrix(W %*% t(proposed_SQUAREM$update$phi))))
          }
          
          square_EM_ll <- evaluate_loglik(beta = proposed_SQUAREM$update$beta, 
            pi = proposed_pi, X = X, y = y, group_mapping = group_mapping, 
            adaptive_weight = adaptive_weight, weights_W = weights_W,
            E.prob =  NA, phi = proposed_SQUAREM$update$phi, W = W, 
            ridge_phi = ridge_phi, ridge_beta = ridge_beta,
            gamma = gamma, Fmatrix = Fmatrix, lambda = lambda, 
            rank_F = rank_F, separate = TRUE)
          
          if (square_EM_ll[1] > new_ll[1]){
            beta <- proposed_SQUAREM$update$beta
            pi <- proposed_pi
            phi <- proposed_SQUAREM$update$phi
            new_ll <- square_EM_ll
            backtrack_counter <- 0
            backtrack_alpha <- proposed_SQUAREM$alpha
          }else{
            backtrack_counter <- 1
            backtrack_alpha <- proposed_SQUAREM$alpha
            while(backtrack_counter <= SQUAREM_backtrack){
              
              backtrack_alpha <- (backtrack_alpha - 1)/2
              
              SQUAREM_update <- lapply(proposed_SQUAREM$terms, FUN=function(i){
                i$init - 2 * backtrack_alpha * i$r + backtrack_alpha^2 * i$v
              })
              
              if (K == 1){
                proposed_pi <- 1 #mean(proposed_postpred)
              }else{
                proposed.loglik.k <- plogis(as.matrix(X %*% SQUAREM_update[[1]]) * 
                                              (2 * y - 1), log.p = TRUE) 
                
                proposed_pi <- colSums((Diagonal(x = as.vector(weights_W)/sum(weights_W)) %*% 
                                          softmax_matrix(W %*% t(SQUAREM_update[[2]]))))
              }
              
              square_EM_ll <- evaluate_loglik(beta = SQUAREM_update[[1]], 
                pi = proposed_pi, 
                X = X, y = y, group_mapping = group_mapping,
                adaptive_weight = adaptive_weight,
                E.prob =  NA, weights_W = weights_W,
                phi = SQUAREM_update[[2]], W = W, 
                ridge_phi = ridge_phi, ridge_beta = ridge_beta,
                gamma = gamma, Fmatrix = Fmatrix, lambda = lambda, rank_F = rank_F, separate = TRUE)
              
              if (square_EM_ll[1] > new_ll[1]){
                break
              }
              
              backtrack_counter <- backtrack_counter + 1
            }
            if (backtrack_counter <= SQUAREM_backtrack){
              beta <- SQUAREM_update[[1]]
              pi <- proposed_pi
              phi <- SQUAREM_update[[2]]
              new_ll <- square_EM_ll
              
            }
          }
        }else{backtrack_counter <- -1; backtrack_alpha <- proposed_SQUAREM$alpha}
        SQUAREM_list <- list()
        SQUAREM_track[it,] <- c(backtrack_alpha, backtrack_counter)
      }
      toc(quiet = quiet_tictoc, log = TRUE)
      
      change_logposterior <- new_ll[1] - running_ll[1]
      running_ll <- compare_ll(new_ll, running_ll, stage = 'SQUAREM', debug_ll = debug)    
    }else{
      change_logposterior <- new_ll[1] - running_ll[1]
      running_ll <- compare_ll(new_ll, running_ll, stage = 'pi', debug_ll = debug)    
    }
    
    
    store_ll[it,] <- running_ll
    old.beta <- beta
    
    max.all <- max(c(change.beta, change.phi))
    # if (change_logposterior < 0){change.logpos <- Inf}
    # print(change_logposterior)
    if (max.all < tolerance.parameters | (change_logposterior > 0 & (change_logposterior < tolerance.logposterior))){
      if (it > control$tau_stabilization){
        message('Converged')
        break
      }
    }
  }
  if (it == iterations){
    message('\n', appendLF = F)
    message('Ended without Convergence')
  }
  

  recons.beta <- basis_M %*% beta
  full_recons <- recons.beta
  
  # recons.beta <<- recons.beta
  # coef_names <<- coef_names
  
  recons.beta <- recons.beta[1:length(coef_names),,drop=F]
  rownames(recons.beta) <- coef_names
  
  tic('Fuse')

  fusion_analysis <- mapply(list_from_cols(recons.beta), 1:ncol(recons.beta), SIMPLIFY = FALSE, FUN=function(b, cl){
    out <- prepare_fusion(factor_levels, term_position, coef_names, beta = b)
    out$cluster <- cl
    return(out)
  })
  fusion_analysis <- do.call('rbind', fusion_analysis)
  attributes(fusion_analysis)$factor_levels <- factor_levels
  toc(quiet = quiet_tictoc, log = TRUE)

  tic('Calculate Information Criterion')
  
  if (lambda == 0){
    binding_restrictions <- NULL
  }else{
    binding_restrictions <- lapply(E.tau, FUN=function(etau.k){
      lapply(etau.k, FUN=function(k){which(k[,1] >= tau_truncate)})
    })
  }
  
  
  if (control$calc_df){
    
    if (lambda == 0){
      
      binding_null_basis <- Diagonal(n = ncol(blocked_X))
      
      df_fp <- NA
      
      # rEw <<- raw.E.ridge
      
      E.ridge <- mapply(raw.E.ridge, pi, SIMPLIFY = FALSE, FUN=function(r, p){
        r * p^gamma
      })
      
      
      if (single_intercept){
        blocked_E <- bdiag(c(0, E.ridge))
        blocked_E <- blocked_E[-(1 + p_X * 1:K),-(1 + p_X * 1:K)]
      }else{
        blocked_E <- bdiag(E.ridge)
      }
      
    }else{
      if (tau_method == 'clip' | is.null(binding_null_basis)){
        
        # binding_lookup <- build_binding_lookup(Fmatrix = Fmatrix, factor_levels = factor_levels, term_position = term_position, coef_names = coef_names)
        existing_restrictions <- lapply(1:K, FUN=function(i){
          i <- lapply(1:length(Fmatrix), FUN=function(j){c()})
          names(i) <- names(Fmatrix)
          return(i)
        })

        if (single_intercept){
          blocked_X <- cbind(1, kronecker(diag(K), X[,-ncol(X)]))
        }else{
          blocked_X <- kronecker(diag(K), X)
        }
        
        # binding_restrictions <- apply(basis_M %*% beta, MARGIN = 2, FUN=function(b){
        #   br <- calc_distance_binding(b, binding_lookup)
        #   # br <- lapply(prepare_fusion(factor_levels = factor_levels, term_position = term_position, coef_names = coef_names, beta = b, simplify = FALSE), FUN=function(i){
        #   #   sapply(i, FUN=function(j){max(unlist(j))})
        #   # })
        #   br <- lapply(br, FUN=function(j){which(j <= tau_equality_threshold)})
        #   return(br)
        # })
        
        
        if (single_intercept){
          binding_null_basis <- lapply(binding_restrictions, FUN=function(r_k){
            binding_k <- do.call('rbind', mapply(r_k, Fmatrix, SIMPLIFY = FALSE, FUN=function(i, F_j){
              do.call('rbind', F_j[i])
            }))
            binding_k <- binding_k[,-p_X]
            if (is.null(binding_k)){
              return(sparse_diag(rep(1, p_X - 1)))
            }else{
              drop0(zapsmall(calculate_nullspace_basis(binding_k), clip_tiny))
            }
          })
          binding_null_basis <- bdiag(c(1, bdiag(binding_null_basis)))
          
          blocked_E <- bdiag(c(0, E.ridge))
          blocked_E <- blocked_E[-(1 + p_X * 1:K),-(1 + p_X * 1:K)]
          
        }else{
          binding_null_basis <- lapply(binding_restrictions, FUN=function(r_k){
            binding_k <- do.call('rbind', mapply(r_k, Fmatrix, SIMPLIFY = FALSE, FUN=function(i, F_j){
              do.call('rbind', F_j[i])
            }))
            if (is.null(binding_k)){
              return(sparse_diag(rep(1, p_X)))
            }else{
              drop0(zapsmall(calculate_nullspace_basis(binding_k), clip_tiny))
            }
          })
          binding_null_basis <- bdiag(binding_null_basis)
          
          blocked_E <- bdiag(E.ridge)
        }
        
      }
      

    }

    if (single_intercept){
      flat_beta <- c(beta[p_X, 1], as.vector(beta[-p_X,]))
    }else{
      flat_beta <- as.vector(beta)
    }
    
    if (control$df_method[1] == 'all'){
      control$df_method <- c('free_param', 'IRLS', 'EM') 
    }
    
    total_df <- c()
    df_beta <- c()
    df_name <- c()
    if (lambda == 0){
      control$df_method <- setdiff(control$df_method, 'free_param')
    }
    
    if ('free_param' %in% control$df_method){
      df_fp <- calculate_df_freeparam(blocked_X = blocked_X,
        binding_null_basis = binding_null_basis)
      df_beta <- c(df_beta, df_fp)
      df_name <- c(df_name, 'free_param')
    }
    
    if ('IRLS' %in% control$df_method){
      df_kc <- calculate_df_kc(X = blocked_X, binding_null_basis = binding_null_basis,
                               ridge = blocked_E, beta = matrix(flat_beta))
      df_beta <- c(df_beta, df_kc)
      df_name <- c(df_name, 'IRLS')
    }
    
    if ('EM' %in% control$df_method){
      df_EM <- calculate_df_EM(X = blocked_X, binding_null_basis = binding_null_basis,
                               ridge = blocked_E, omega_weight = as.vector(E.omega))
      df_beta <- c(df_beta, df_EM)
      df_name <- c(df_name, 'EM')
      
    }
    
    total_df <- df_beta + (K - 1) * ncol(W)

    final_ll <- running_ll[2]
    
    est_IC <- data.frame(method = df_name,
                         df = total_df,
                         df_beta = df_beta, stringsAsFactors = F)
    
    est_IC$BIC <- -2 * final_ll + log(nrow(X)) * est_IC$df
    est_IC$BIC_group <- -2 * final_ll + log(length(unique_group)) * est_IC$df
    est_IC$AIC <- -2 * final_ll + 2 * est_IC$df
    
    est_IC$GCV <- -2 * final_ll / (nrow(X) * (1 - est_IC$df/nrow(X))^2)
    est_IC$GCV_group <- -2 * final_ll / (length(unique_group) * (1 - est_IC$df/length(unique_group))^2)
    est_IC$ll <- final_ll
    est_IC$N <- nrow(X)
    est_IC$N_group <- length(unique_group)
    est_IC$K <- K
  }else{
    est_IC <- NULL
  }
  
  toc(quiet = quiet_tictoc, log = TRUE)
  tic('Final Collection')
  
  if (ncol_W != 0){
    group_postpred_prob <- calculate_posterior_zi(loglik.k = loglik.k, group_mapping = group_mapping, W = W, 
         K = K, ncol_W = ncol_W, pi = pi, phi = phi, return = 'postpred')
    group_postpred_prob <- data.frame(group = unique_group, as.matrix(group_postpred_prob), stringsAsFactors = F)
    names(group_postpred_prob)[-1] <- paste0('cluster_', 1:K)
  }else{
    group_postpred_prob <- NA
  }
  
  group_output <- data.frame(group = unique_group, as.matrix(group_E.prob), stringsAsFactors = F)
             
  names(group_output)[-1] <- paste0('cluster_', 1:K)
  
  if (is.na(control$return_data)){
    data_list <- NULL
  }else if (control$return_data){
    data_list <- list(design = design, X = X, y = y,  W = W, rank_F = rank_F,
      Fmatrix = Fmatrix, adaptive_weight = adaptive_weight,
      group_E.prob = group_E.prob, basis_M = basis_M,
      weights = weights, weights_W = weights_W,
      group = group, b_r = update_mod$b_r)
    data_list$E_ridge <- E.ridge
  }else{
    data_list <- list(design = design)
  }
  
  internal_parameters <- list(
    ordered_factors = ordered_factors,
    W = list(args = args_W), 
    rare = list(rare_fmt_col = rare_fmt_col, rare_col = rare_col),
    unique_choice = unique_choice,
    interactions = do_interactions,
    weights = list(weights = weights, weights_W = weights_W),
    single_intercept = single_intercept,
    group = list(null_group = null_group), 
    misc = list(clip_tiny = clip_tiny, nullspace = tracking_restrictions))
  internal_parameters$control <- control
  internal_parameters$SQUAREM <- list(SQUAREM_track = SQUAREM_track)
  internal_parameters$data <- data_list
  internal_parameters$adaptive_weight <- adaptive_weight

  internal_parameters$convergence <- list(beta = change.beta, phi = change.phi, 
      log.posterior = change_logposterior)
  internal_parameters$trajectory <- list(beta = store_beta, ll = store_ll)
  
  internal_parameters$fusion <- fusion_analysis
  internal_parameters$factor_levels <- factor_levels
  internal_parameters$formula <- list(het = formula_recons, mod = formula_mod,
                                      weights = formula_weight,
                                      other_parameters = conjoint_names)
  internal_parameters$group <- list(unique_groups = unique_group, 
    group_mapping = group)
  internal_parameters$penalty <- penalty_for_regression
  internal_parameters$Fmatrix <- Fmatrix
  internal_parameters$use_forced_choice <- use_forced_choice
  internal_parameters$basis_M <- basis_M
  
  parameters <- list(beta = recons.beta, pi = pi, 
      eff_lambda = lambda, orig_lambda = orig_lambda,
      phi = phi, nullspace_beta = beta, 
      gamma = gamma, full_recons = full_recons)
  
  posterior <- list('posterior' = group_output, 
                     'posterior_predictive' = group_postpred_prob)
  
  output <- list(parameters = parameters, K = K,
                 posterior = posterior, 
                 information_criterion = est_IC, 
                 internal_parameters = internal_parameters)
  class(output) <- 'FactorHet'
  toc(quiet = quiet_tictoc, log = TRUE)
  
  tic('Calculate SE')
  if (control$calc_se){
    if (!(control$return_data %in% TRUE)){
      orig_data <- output$internal_parameters$data
      output$internal_parameters$data <- list(design = design, X = X, y = y,  W = W, rank_F = rank_F,
                                              Fmatrix = Fmatrix, adaptive_weight = adaptive_weight,
                                              group_E.prob = group_E.prob, basis_M = basis_M,
                                              group = group, b_r = update_mod$b_r, binding_restrictions = binding_restrictions)
      output$vcov <- estimate.vcov(output)
      output$internal_parameters$data <- orig_data
    }else{
      output$vcov <- estimate.vcov(output)
    }
  }
  toc(quiet = quiet_tictoc, log = TRUE)
  
  tic.clear()
  timing_list <- unlist(tic.log())
  timing_list <- do.call('rbind', strsplit(timing_list, split=': | sec elapsed', perl = TRUE))
  timing_list <- do.call('rbind', lapply(split(timing_list[,2], timing_list[,1]), FUN=function(i){data.frame(t(c(length(i), as.vector(summary(as.numeric(i))))))}))
  colnames(timing_list) <- c('number_of_times', 'min', 'first_q', 'median', 'mean', 'third_q', 'max')
  timing_list <- cbind('stage' = rownames(timing_list), timing_list)
  timing_list$total_time <- with(timing_list, number_of_times * mean)
  rownames(timing_list) <- NULL
  
  tic.clearlog()
  
  output$internal_parameters$timing <- timing_list
  
  return(output)
}

# List from matrix columns
list_from_cols <- function(M){
  lapply(1:ncol(M), FUN=function(i){M[,i]})
}

F_prior_weight <- function(beta, Fmat, rank_F, lambda, gamma, pi, adaptive_weight, kern_epsilon = 0){
  #(lambda pi_k^gamma)^{rank(F)} exp(- \lambda pi^gamma sqrt(b^T F b))
  #rank(F) * [log(lambda) + gamma log(pi_k)] 
  #-lambda * pi^gamma sqrt(b^T F b)
  #kernel_weight <- apply(beta, MARGIN = 2, F_prior_kernel, Fmat = Fmat, log = TRUE)
  kernel_weight <- mapply(list_from_cols(beta), adaptive_weight, FUN=function(b, aw){
    F_prior_kernel(beta = b, Fmat = Fmat, aw = aw, kern_epsilon = kern_epsilon, log = TRUE)
  })
  logprior_kernel <- sum(lambda * pi^gamma * kernel_weight)
  
  logprior_normcons <- rank_F * log(lambda) + rank_F * gamma * log(pi)
  logprior_normcons <- sum(logprior_normcons)
  return(logprior_kernel + logprior_normcons)
  
}

F_prior_kernel <- function(beta, Fmat, aw, kern_epsilon = 0, disagg = FALSE, log = TRUE){
  ncol_beta <- ncol(beta)
  if (is.null(ncol_beta)){ncol_beta <- 1}
  if (ncol_beta > 1){
    stop('beta must be vector or n x 1 matrix')
  }
  weight <- mapply(Fmat, aw, FUN=function(F.j, exo_j){
    weight.j <- sapply(F.j, FUN=function(l){
      as.vector(t(beta) %*% l %*% beta) + kern_epsilon
    })
    #Correct the possibility of very small negative weights
    weight.j[weight.j < 0 & weight.j > -1e-10] <- 0
    if (any(weight.j < 0)){
      stop('Nontrivial Negative weight j in F_prior_kernel. Estimation terminated.')
    }
    return(sum(exo_j * sqrt(weight.j)))
  })
  if (disagg){return(weight)}
  log.prior <- - sum(weight)
  if (log){
    return(log.prior)
  }else{
    return(exp(log.prior))
  }
}

# Internal function to do EM update
# Takes: beta, Fmatrix, lambda, pi.gamma and
# weights
F_EM_update <- function(beta, Fmat, lambda, pi.gamma, exo_weights){
  
  if (is.null(exo_weights)){
    exo_weights <- as.list(rep(1, length(Fmat)))
  }

  EM_weight <- mapply(Fmat, exo_weights, SIMPLIFY = FALSE, FUN=function(F.j, exo.j){
    weight.j <- sapply(F.j, FUN=function(l){
      as.vector(t(beta) %*% l %*% beta)
    })
    weight.j[weight.j < 0 & weight.j > -1e-10] <- 0
    weight.j <- sqrt(weight.j)
    if (any(is.na(weight.j))){stop('Negative weight j in E-Step')}
    effective_lambda <- lambda * pi.gamma * exo.j
    EM.tau <- cbind(effective_lambda / weight.j, weight.j / (effective_lambda) + 1/(effective_lambda)^2)
    colnames(EM.tau) <- c('Einv.tau', 'Etau')
    return(EM.tau)
  })
  return(EM_weight)
}

# Get the log posterior (un-normalized) of the (i) observed likelihood and (ii) the
# single augmented log posterior to monitor convergence
# 
evaluate_loglik <- function(beta, pi, phi, weights_W,
    X, y, group_mapping, W, adaptive_weight,
    gamma, Fmatrix, lambda, ridge_phi, ridge_beta,
    kern_epsilon = 0,
    calc_single = FALSE, rank_F = NULL, E.prob = NULL, separate = FALSE){
  
  if (length(pi) != ncol(beta)){
    stop('pi must be a K-length vector')
  }
  #If no rank(F) is provided, add. This is necessary for \gamma > 0 to correctly
  #evaluate the log-prior on \beta that depends on \pi.
  if (is.null(rank_F)){
    rank_F <- Reduce("+", lapply(Fmatrix, FUN=function(k){Reduce("+", k)}))
    rank_F <- rank_via_null(rank_F)
  }
  ncol_W <- ncol(W)
  #Get the loglikelihood for y_i given cluster assignment.
  xb <- as.matrix(X %*% beta) 
  K <- ncol(beta)
  
  loglik.k <- plogis(xb * (2 * y - 1), log.p = TRUE) 
  
  group_loglik.k <- calculate_posterior_zi(
    loglik.k = loglik.k, group_mapping = group_mapping, K = K, ncol_W = ncol_W, 
    pi = pi, phi = phi, return = 'loglik',W = W
  ) 

  #Evaluate the log-likelihood (observed)
  loglik.obs <- apply(group_loglik.k, MARGIN =1, FUN=function(i){
    #log(sum(exp(i)))
    #log(sum_k L_k)
    #log(L_max * sum_k L_k/L_max)
    #ll_max + log(sum_k L_k/L_max)#
    #ll_max + log(sum_k exp(ll_k - ll_max))
    max.i <- max(i)
    return(max.i + log(sum(exp(i - max.i))))
  })
  
  loglik.obs <- sum(weights_W * loglik.obs)

  if (lambda == 0){
    logprior.beta <- -1/2 * apply(beta, MARGIN = 2, 
      FUN=function(i){as.numeric(t(i) %*% ridge_beta %*% i)})
    logprior.beta <- sum(pi^gamma * logprior.beta) +
      sum(rank_F * gamma * log(pi)) 
  }else{
    #Get the log prior for each cluster p(\beta_r) = c(F) (\lambda pi_k^\gamma)^m \exp(- \sum_{l} \lambda \pi_r^\gamma\sqrt{\beta_r^T F_l \beta_r})
    logprior.beta <- F_prior_weight(beta = beta, Fmat = Fmatrix, lambda = lambda, 
                                    adaptive_weight = adaptive_weight,
                                    kern_epsilon = kern_epsilon,
                                    pi = pi, gamma = gamma, rank_F = rank_F)
  }
  if (K != 1 & ridge_phi > 0 & ncol_W != 0){
    ridge_penalty <- ridge_phi * make_TMatrix(K)
    zeroed_phi <- phi[-1,,drop=F]
    logprior.phi <- -1/2 * sum(diag(t(zeroed_phi) %*% ridge_penalty %*% zeroed_phi))
  }else{
    logprior.phi <- 0
  }
  
  logposterior.obs <- loglik.obs + logprior.beta + logprior.phi
  
  if (separate){
    return(c(logposterior.obs, loglik.obs, logprior.phi + logprior.beta))
  }else{
    return(logposterior.obs)
    
  }
}

clip_Etau <- function(E.tau, threshold){
  lapply(E.tau, FUN=function(E.tauk){
    lapply(E.tauk, FUN=function(l){
      l[,1] <- ifelse(l[,1] > threshold, threshold, l[,1])
      return(l)
    })
  })
}

compare_ll <- function(new_ll, running_ll, stage, debug_ll, tol = 1e-8){
  if (running_ll[1] > new_ll[1] + tol){
    if (debug_ll){
      print('New')
      print(new_ll)
      print('Running')
      print(running_ll)
      print(new_ll - running_ll)
      stop(paste0(stage, ': Decrease in observed loglik'))
    }
  }
  return(new_ll)
}

make_ridge <- function(E.tauk, Fmatrix = NULL, raw_F = NULL, method = 'unlist'){
  if (method == 'old'){
    ridge.k <- Reduce('+', mapply(E.tauk, Fmatrix, SIMPLIFY = FALSE, FUN=function(etau, F.j){
      #Sum over the L (l) F matrices for each factor.
      ridge.j <- Reduce('+', mapply(etau[,1], F.j, SIMPLIFY = FALSE, FUN=function(etau.l, l){
        etau.l * l
      }))
    }))
  }else if (method == 'unlist'){
    unlist_tau <- unlist(lapply(E.tauk, FUN=function(i){i[,1]}))
    unlist_F <- unlist(Fmatrix)
    
    L <- length(unlist_tau)
    
    ridge.k <- sparseMatrix(i = 1, j = 1, x = 0, dims = dim(unlist_F[[1]]))
    for (l in 1:L){
      ridge.k <- ridge.k + unlist_tau[l] * unlist_F[[l]]
    }
    ridge.k <- drop0(ridge.k)
  }else if (method == 'raw'){
    unlist_tau <- unlist(lapply(E.tauk, FUN=function(i){i[,1]}))  
    raw_size <- raw_F$size
    if (length(raw_size) != length(unlist_tau)){stop('Raw Method Misaligned')}  
    raw_matrix <- raw_F$raw
    raw_matrix[,3] <- raw_matrix[,3] * rep(unlist_tau, times = raw_size)
    ridge.k <- sparseMatrix(i = raw_matrix[,1] + 1, j = raw_matrix[,2] + 1, x = raw_matrix[,3], dims = raw_F$dim)
  }else{stop()}
  return(ridge.k)
}

calculate_posterior_zi <- function(loglik.k, group_mapping, K, 
   W, ncol_W, pi, phi, return = 'prob'){
  #Get the summed \sum_t ll(y_{it} | \beta_r) for each group
  if (K == 1){
    if (return == 'prob'){
      return(Matrix(1, nrow = nrow(W), ncol = 1))
    }else if (return == 'loglik'){
      group_loglik.k <- apply(loglik.k, MARGIN = 2, FUN=function(i){as.vector(i %*% group_mapping)})
      return(group_loglik.k)
    }else if (return == 'postpred'){
      return(NA)
    }
  }
  
  group_loglik.k <- apply(loglik.k, MARGIN = 2, FUN=function(i){as.vector(i %*% group_mapping)})
  
  #Add in pi.
  if (ncol_W == 0){
    log.pi <- log(pi)
    group_loglik.k <- logpi_adjust(group_loglik.k, log.pi)
    group_postpred_prob <- NULL
  }else{
    group_postpred_prob <- softmax_matrix(W %*% t(phi))
    group_loglik.k <- group_loglik.k + log(group_postpred_prob)
  }
  
  if (return == 'prob'){
    #Get E[z_i | -] for cluster membership
    group_E.prob <- softmax_matrix(group_loglik.k)
    
    colnames(group_E.prob) <- 1:K
    return(group_E.prob)
  }else if (return == 'loglik'){
    return(group_loglik.k)
  }else if (return == 'postpred'){
    return(group_postpred_prob)
  }else{stop('Invalid return type')}
}

update_beta <- function(X, y, E.omega, obs.E.prob, E.ridge, weights, K, global_int, blocked_X = NULL, p_X = NULL, prior_beta = NULL, cg_it = NULL, method='cpp'){
  
  #If global intercept, use dedicated function with "cg" or "cpp".  
  if (global_int){
    
    new_beta <- update_beta_global_int(blocked_X = blocked_X, p_X = p_X, 
      y = y, method = method, prior_beta = prior_beta, weights = weights,
      cg_it = cg_it, E.omega = E.omega, obs.E.prob = obs.E.prob, E.ridge = E.ridge, K = K)

    return(new_beta)
  }
  #Otherwise, do separately.
  if (method == 'cpp'){
    
    new_beta <- cpp_beta(K = K, X = X, E_ridge = E.ridge, y = y, weights = weights,
        E_omega = as.matrix(E.omega), obs_E_prob = as.matrix(obs.E.prob))

  }else if (method == 'cg'){
    
    if (is.null(cg_it)){
      cg_it <- 0
    }
    new_beta <- cg_custom(K = K, X = X, 
              list_ridge = E.ridge, 
              omega = as.matrix(E.omega), 
              s = weights * (y - 1/2), old_beta = prior_beta, 
              weights = as.matrix(obs.E.prob), tol = sqrt(.Machine$double.eps), it_max = cg_it)
    new_beta <- new_beta$beta
    
    # new_beta <- cpp_beta_CG(K = K, X = X, E_ridge = E.ridge, y = y, E_omega = E.omega, tol = sqrt(.Machine$double.eps),
    #     old_beta = prior_beta, obs_E_prob = obs.E.prob, it_max = cg_it)
    # new_beta <- new_beta$flat_beta

  }else if (method == 'base'){
    
    new_beta <- matrix(0, nrow = ncol(X), ncol = K)
    for (k in 1:K){
      #Using standard Polya-Gamma updates, do least squares with ridge.
      s.ik <- obs.E.prob[,k] * (y-1/2) * weights
      omega.ik <- E.omega[,k]
      #solve(t(X) %*% X + ridge, t(X) %*% y) for the "pure ridge" version.
      new_beta[,k] <- as.vector(solve(t(X) %*% Diagonal(x = omega.ik) %*% X + E.ridge[[k]], t(X) %*% s.ik))
    }
  }else{stop('Invalid method for update beta')}
  
  return(new_beta)
}


#' @importFrom Matrix bdiag
update_beta_global_int <- function(blocked_X, y, E.omega, obs.E.prob, weights, E.ridge, K, p_X, cg_it = NULL, prior_beta = NULL, method = 'cpp'){
  
  if (method == 'cpp'){
    #Exclude Intercept
    blocked_E <- bdiag(c(0, E.ridge))
    blocked_E <- blocked_E[-(1 + p_X * 1:K),-(1 + p_X * 1:K)]


    blocked_beta <- cpp_beta_plain(X = blocked_X,
        s = rep(weights * (y - 1/2), K) * as.vector(obs.E.prob),
        K = K, omega = sparse_diag(as.vector(E.omega)),
        ridge = blocked_E
    )

    mu <- blocked_beta[1]
    
    beta <- matrix(blocked_beta[-1], ncol = K)
    beta <- rbind(beta, mu)
  }else if (method == 'cg'){
    
    blocked_E <- bdiag(c(0, E.ridge))
    blocked_E <- blocked_E[-(1 + p_X * 1:K),-(1 + p_X * 1:K)]
    
    flat_beta <- c(prior_beta[p_X, 1], as.vector(prior_beta[-p_X,]))
    
    if (is.null(cg_it)){
      cg_it <- 0
    }
    
    cg_b <- cg_custom(K = 1, X = blocked_X,
      s = rep((y-1/2) * weights, K) * as.vector(obs.E.prob), 
      list_ridge = list(blocked_E), tol = sqrt(.Machine$double.eps),
      it_max = cg_it, old_beta = flat_beta, weights = matrix(0, nrow = 1, ncol = 0),
      omega = matrix(as.vector(E.omega))
    )

    # cg_b <- eigen_cg_with_ldlt(X = blocked_X,
    #   s = rep(y - 1/2, K) * as.vector(obs.E.prob),
    #   omega = as.vector(E.omega),
    #   ridge = blocked_E, tol = sqrt(.Machine$double.eps),
    #   it_max = cg_it, old_beta = flat_beta)
    # print(cg_b$iter)
    mu <- cg_b$beta[1]
    
    beta <- matrix(cg_b$beta[-1], ncol = K)
    beta <- rbind(beta, mu)
    
  }else{
    blocked_E <- bdiag(c(0, E.ridge))
    blocked_E <- blocked_E[-(1 + p_X * 1:K),-(1 + p_X * 1:K)]
    
    blocked_beta <- solve(t(blocked_X) %*% Diagonal(x = as.vector(E.omega)) %*% blocked_X + blocked_E, 
          t(blocked_X) %*% as.vector(rep(y - 1/2, K) * weights * as.vector(obs.E.prob)))
    
    mu <- blocked_beta[1]
    
    beta <- matrix(blocked_beta[-1], ncol = K)
    beta <- rbind(beta, mu)
    
  }
  
  return(beta)
}


calculate_null_binding <- function(Fmatrix, beta, tolerance = 1e-6){
  
  l2_beta <- mapply(list_from_cols(beta), SIMPLIFY = FALSE, FUN=function(b){
    lapply(F_EM_update(b, Fmatrix, lambda = 1, pi.gamma = 1, exo_weights = NULL), FUN=function(i){1/i[,1]})
  })
  K <- ncol(beta)
  nonbinding_ridge <- as.list(rep(NA, K))
  all_binding_null <- as.list(rep(NA, K))
  
  for (k in 1:K){
    binding_F <- mapply(l2_beta[[k]], Fmatrix, SIMPLIFY = FALSE, FUN=function(pen, F.j){
      enforce_restriction <- (pen <= tolerance)
      
      enforce_F <- do.call('rbind', F.j[enforce_restriction])
      
      ridge.j <- Reduce('+', mapply(pen[!enforce_restriction], F.j[!enforce_restriction], SIMPLIFY = FALSE, FUN=function(etau.l, l){
        etau.l * l
      }))
      if (is.null(ridge.j)){
        ridge.j <- drop0(sparseMatrix(i = 1, j = 1, x = 0, dims = dim(F.j[[1]])))
      }
      return(list(enforce = enforce_F, ridge = ridge.j))
    })
    
    null_binding <- lapply(binding_F, FUN=function(i){i$enforce})
    null_binding <- null_binding[sapply(null_binding, FUN=function(i){!is.null(i)})]
    if (length(null_binding) == 0){
      all_binding_null[[k]] <- 'No Restrictions'
    }else{
      all_binding_null[[k]] <- drop0(zapsmall(calculate_nullspace_basis(do.call('rbind', null_binding)), 15))
    }
    nonbinding_ridge[[k]] <- Reduce('+', lapply(binding_F, FUN=function(i){i$ridge}))
  }
  return(list(null_basis = all_binding_null, nonbinding_ridge = nonbinding_ridge))
}

prepare_SQUAREM <- function(SQUAREM_list, alpha_revert = -1.01, scale_alpha = 1){
  
  SQUAREM_terms <- lapply(c('beta', 'phi'), FUN=function(i){
    #Get change and change in change
    r <- SQUAREM_list[[2]][[i]] - SQUAREM_list[[1]][[i]]
    v <- SQUAREM_list[[3]][[i]] - SQUAREM_list[[2]][[i]]
    v <- v - r
    return(list(init = SQUAREM_list[[1]][[i]], r = r, v = v, norm = c(sum(r^2), sum(v^2))))
  })
  
  SQUAREM_norm <- sqrt(rowSums(sapply(SQUAREM_terms, FUN=function(i){i$norm})))
  alpha <- -SQUAREM_norm[1]/SQUAREM_norm[2] * scale_alpha
  if (!is.finite(alpha)){
    alpha <- -1
  }else if (alpha > -1){
    alpha <- alpha_revert
  }

  SQUAREM_update <- lapply(SQUAREM_terms, FUN=function(i){
    i$init - 2 * alpha * i$r + alpha^2 * i$v
  })
  names(SQUAREM_update) <- c('beta', 'phi')
  return(list(alpha = alpha, update = SQUAREM_update, terms = SQUAREM_terms))
}

simple_logit <- function(y, X, iterations, weights, obs.E.prob = NULL, binary = NULL, beta_method , beta_cg_it = NULL, prec_ridge = 1){
  
  if (!all(binary %in% c(0,1))){
    stop('simple_logit only works for a binary "in out" decision')
  }

  if (is.null(beta_cg_it) & beta_method == 'cg'){
    beta_cg_it <- 0
  }
  
  if (!is.null(binary)){
    X <- X[binary,]
    y <- y[binary]
    
  }
  if (all(X[,ncol(X)] == 1)){
    mean_init <- qlogis(mean(y))
    if (!is.finite(mean_init)){
      mean_init <- 0
    }
    beta <- c(rep(0, ncol(X) - 1), mean_init)
  }else{
    beta <- rep(0, ncol(X))
  }
  
  if (is.null(obs.E.prob)){
    obs.E.prob <- matrix(1, nrow = length(y))
  }
  
  E.ridge <- list(sparse_diag(rep(prec_ridge, ncol(X))))
  
  lagged_obj <- -Inf
  
  for (it in 1:iterations){
    xb <- as.vector(X %*% beta)
    
    ll <- sum(plogis((2 * y - 1) * xb, log.p = TRUE) * obs.E.prob)
    lprior <- -1/2 * sum(beta^2 * prec_ridge)
    obj <- sum(ll + lprior)
    
    diff <- obj - lprior
    lagged_obj <- obj
    if (diff > 0 & diff < 1e-3){
      break
    }
    
    E.omega <- obs.E.prob/(2 * xb) * tanh(xb / 2)
    E.omega[abs(xb) < 1e-10] <- 1/4    
    E.omega <- E.omega * weights
    
    beta <- update_beta(X=X,y=y,E.omega = E.omega, 
                weights = weights,
                obs.E.prob = obs.E.prob, cg_it = beta_cg_it,
                E.ridge = E.ridge, K= 1, method = beta_method, global_int = F, prior_beta = beta)
    
  }
  
  return(beta)
}
