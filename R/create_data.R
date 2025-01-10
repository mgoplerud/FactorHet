# Parse the formula
prepare_formula <- function(fmla_main, fmla_moderator = NULL, weights = NULL,
  design, group = NULL, choice_order = NULL, task = NULL, delete_response = FALSE){
  
  #
  # Parse MAIN formula:
  #
  # Parse the formula, pretending "." is character
  terms_fmla <- terms(fmla_main, allowDotAsName = T)
  if (delete_response){
    terms_fmla <- delete.response(terms_fmla)    
  }
  all_terms <- attr(terms_fmla, 'term.labels')
  # Remove the interactions and the "."
  # so only select the "base"
  factor_names <- raw_factor_names <- all_terms[!grepl(all_terms, pattern=':|^\\.$')]
  factor_names <- gsub(factor_names, pattern='^`|`$', replacement = '')
  
  # Get outcome
  if (!delete_response){
    outcome <- attr(terms_fmla, 'response')
    if (outcome == 0){stop('No response provided to formula')}
    outcome <- rownames(attr(terms_fmla, 'factors'))[outcome]
  }else{
    outcome <- NULL
  }
  
  # Parse Weights Formula
  if (inherits(weights, 'formula')){
    weights <- model.frame(update.formula(weights, ' ~ 0 + .'), design, na.action = 'na.pass')
    if (ncol(weights) != 1){stop('weights formula must contain only one term.')}
    name_weights <- colnames(weights)
    weights <- weights[,1]
  }else if (is.null(weights)){
    name_weights <- NULL
    weights <- rep(1, nrow(design))
  }else{
    stop('weights must be a formula for a column in the design.')
  }
  
  # Parse MODERATOR formula:
  if (!is.null(fmla_moderator)){
    terms_mod <- terms(fmla_moderator, data = design)
    mod_terms <- rownames(attr(terms_mod, 'factors'))
    mod_terms <- gsub(mod_terms, pattern='^`|`$', replacement = '')
    
  }else{
    terms_mod <- terms(~ 1)
    mod_terms <- NULL
  }
  if (length(base::intersect(mod_terms, factor_names)) > 0){
    stop('No factors can be used as moderators!')
  }
  
  null_or_parse <- function(x){
    if (is.null(x)){
      return(NULL)
    }else if (inherits(x, 'formula')){
      return(attr(terms(x), 'term.labels'))
    }else if (inherits(x, 'character')){
      return(x) 
    }else{
      stop('invalid type to task, group or choice_order')
    }
  }
  # Deal with group, task, choice
  name_group <- null_or_parse(group)
  name_task <- null_or_parse(task)
  name_choice <- null_or_parse(choice_order)
  
  if (any(c(name_group, name_task, name_choice) %in% c(mod_terms, factor_names))){
    stop('group, task, and choice must not be in formula or moderator')
  }
  
  if ('.sw' %in% colnames(design)){
    stop('design must not contain a column named ".sw" to ensure weights do not overwrite.')
  }
  design$.sw <- weights
  # Select only the relevant columns of the design
  design <- design[, 
    c(factor_names, mod_terms, 
      outcome, name_group, name_task, name_choice, '.sw')
  ]
  # Remove NA
  design <- na.omit(design)
  
  # Save out the weights after removing NA and remove placeholder column.
  weights <- design[[".sw"]]
  design <- design[, names(design) != ".sw", drop = FALSE]
  
  if (nrow(design) == 0){stop('design has no observations after listwise deletion.')}
  #  
  # Prepare the Interactions
  #
  # Get the full terms object (filling in ".")
  terms_fmla <- terms(fmla_main, data = design)
  # Get all interactions
  do_interactions <- attr(terms_fmla, 'term.labels')[grepl(attr(terms_fmla, 'term.labels'), pattern=':')]
  do_interactions <- strsplit(do_interactions, split=':')
  # Select only those interactions that are
  # in the main effects
  do_interactions <- lapply(do_interactions, FUN=function(i){
    if (all(i %in% raw_factor_names)){
      return(i)
    }else{
      return(NULL)
    }
  })
  do_interactions <- raw_do_interactions <- do_interactions[!sapply(do_interactions, is.null)]
  do_interactions <- lapply(do_interactions, FUN=function(i){gsub(i, pattern='^`|`$', replacement = '')})
  # Reconstruct the exact formula used
  formula_recons <- paste(outcome, ' ~ ', 
                          paste(c(raw_factor_names, sapply(raw_do_interactions, paste, collapse=':')), collapse='  + '))
  formula_mod <- paste(' ~ ', paste(attr(terms_mod, 'term.labels'), collapse = ' + '))
  if (formula_mod == ' ~  '){formula_mod <- ' ~ 1'}
  output <- list(design, outcome, factor_names, formula_recons, formula_mod, do_interactions, 
                 name_group, name_task, name_choice, weights, name_weights)
  
  names(output) <- c('design', 'outcome', 'factor_names', 
                     'formula_recons', 'formula_mod', 'do_interactions', 
                     'name_group', 'name_task', 'name_choice_order', 'weights',
                     'name_weights')
  
  # Check whether group/task/choice_order are correctly specified
  if (length(name_group) == 0){
    if (length(name_task) != 0){
      stop('task may not be provided without group.')
    }
  }else{
    if (length(name_task) != 0){
      max_per_group <- max(table(paste0(design[[name_group]], '@~@~@', design[[name_task]])))
    }else{
      max_per_group <- max(table(design[[name_group]]))
    }
    if (max_per_group != 1 & length(name_choice) == 0){
      stop('More than one observation per group/task combination! task must be provided if multiple observations per group. choice_order should be provided if conjoint design.')
    }
    if (max_per_group != 2 & length(name_choice) == 1){
      stop('Not all group/task combinations have two observations. task must be provided if multiple observations per group.')
    }
  }
  
  return(output)
}

#' Create the (Sparse) Design Matrix for Analysis
#' 
#' Create a sparse design matrix with weighted sum to zero constraints.
#' 
#' @keywords internal
#' @param design Data frame with levels of each factor assigned to observation i
#' @param penalty_for_regression An object from "create_penalty"
#' @import Matrix
create_data <- function(design, penalty_for_regression, warn_missing = TRUE, verif_row = TRUE, remove_cols = Inf){
  if (nrow(design) == 0){stop('No observations!')}
  coef_names <- penalty_for_regression$coef
  J_names <- penalty_for_regression$J_names
  J <- penalty_for_regression$J
  term_position <- penalty_for_regression$term_position
  make_interactions <- penalty_for_regression$make_interactions
  
  n_obs <- nrow(design)
  X <- matrix(nrow = 0, ncol = 2)
  
  design <- data.frame(design, check.names = FALSE)[penalty_for_regression$J_names]
  # design <- as.data.frame(select(.data = data.frame(design, check.names = FALSE), 
  #   penalty_for_regression$J_names))
  
  for (j in 1:J){
    match_main <- match(paste0(J_names[j], '(', design[, j], ')'), coef_names)
    if (any(is.na(match_main))){
      if (warn_missing){warning(paste0('Levels found in data for factor ', J_names[j], ' not found in training data.'))}
      X <- rbind(X, cbind(seq_len(n_obs)[which(!is.na(match_main))], match_main[which(!is.na(match_main))]))
    }else{
      X <- rbind(X, cbind(1:n_obs, match_main))
    }
  }
  
  for (int in make_interactions){
    sort_int <- sort(match(int, J_names))
    
    j <- sort_int[1]
    j.prime <- sort_int[2]
    
    id <- make_pad_id(data = design, j = j, j.prime = j.prime, J = J)
    #Get the IDs that should not be removed
    valid_id <- !(id %in% remove_cols)
    #Match position     
    id <- match(id, coef_names)

    if (any(is.na(id[valid_id]))){
      stop(paste0('Matching failed for factor ', J_names[j],' and ',J_names[j.prime]))
    }
    
    X <- rbind(X, cbind(seq_len(length(id))[valid_id], id[valid_id]))
  }
  
  X <- rbind(X, cbind(1:n_obs, rep(1, n_obs)))
  X <- sparseMatrix(i = X[,1], j = X[,2], x = 1,
                    dims = c(n_obs, length(coef_names)),
                    dimnames = list(NULL, coef_names))
  if (verif_row){
    rowsums_X <- rowSums(X)
    rowchecksum <- 1 + J + length(make_interactions)
    if (any(rowsums_X != rowchecksum)){
      stop('Creating design matrix failed: Not all in each row')    
    }
  }
  return(X)
}


# From a design matrix, create a "forced choice"
# between two choices.
make_forced_choice <- function(y, X, group, task, weights, unique_choice, split_character = '@@',
                               choice_order, randomize, estimate_intercept){
  if (any(is.na(choice_order))){stop('No missingness permitted for "choice_order"')}
  if (any(is.na(task))){stop('No missingness permitted for "task"')}
  if (any(is.na(group))){stop('No missingness permitted for "group"')}
  #Turn task to numeric
  choice_order <- match(choice_order, unique_choice)
  #Split the observations into person and task
  #Check whether split character is valid
  is_valid_split <- any(grepl(group, pattern=split_character)) | any(grepl(group, pattern=split_character))
  if (!is_valid_split){
    #Generate random split
    split_character <- paste(c('@', sample(c(letters, LETTERS), 25, replace = T), '@'), collapse = '')
  }
  group_task_id <- paste0(group, split_character, task)
  
  split_person_task <- split(choice_order, group_task_id)
  
  if (length(unique_choice) != 2){
    stop('Must have only two choices for forced choice tasks')
  }
  
  if (any(lengths(split_person_task) != 2)){
    easy_message(names(split_person_task)[which(lengths(split_person_task) != 2)])
    stop('Forced Choice Requires ALL group-tasks to have exactly two outcomes. "Bad" combinations listed above')
  }
  
  if (randomize){
    split_person_task <- lapply(split_person_task, FUN=function(i){sample(1:2)})
  }  
  split_X_by_row <- split(1:nrow(X), group_task_id)
  if (!is.null(y)){
    
    split_y <- split(y, group_task_id)
    
    if (any(sapply(split_y, anyDuplicated) > 0)){
      stop('Forced choice must have exactly one "yes" and "no" for each group-task combination.')
    }
    #Second is "yes"
    collapsed_y <- mapply(split_y, split_person_task, FUN=function(id_y, id_gt){
      id_y[id_gt == 2] == 1
    })
    collapsed_y <- as.numeric(collapsed_y)
    names(collapsed_y) <- NULL
    
  }else{
    collapsed_y <- NULL
  }

  collapsed_weights <- mapply(split(weights, group_task_id), names(split_person_task), FUN=function(i, p){
    if (i[1] != i[2]){
      stop(paste0('Weights for (', p, ') are not identical.'))
    }
    return(sum(i)/2)
  })

  first_task <- mapply(split_X_by_row, split_person_task, FUN=function(id_x, id_gt){
    id_x[id_gt == 1]
  })
  second_task <- mapply(split_X_by_row, split_person_task, FUN=function(id_x, id_gt){
    id_x[id_gt == 2]
  })
  
  
  #Stop if NOT paired
  stopifnot(identical(names(first_task), names(second_task)))
  #Linear predictor difference from SECOND task minus FIRST
  #If SECOND is much more positive -> vote "1".
  differenced_X <- X[second_task, ] - X[first_task,]
  if (estimate_intercept){
    differenced_X[,1] <- 1  
  }else{
    differenced_X[,1] <- NA
  }
  
  group_task <- do.call('rbind', strsplit(names(first_task), split=split_character))
  
  return(list(X = differenced_X, y = collapsed_y, weights = collapsed_weights,
              group = group_task[,1], task = group_task[,2], estimate_intercept = estimate_intercept))  
}

make_pad_id <- function(data, j, j.prime, J, design_method = 'column'){
  pre_pad_NA <- rep(NA, j - 1)
  pre_pad_NA <- paste0(pre_pad_NA, collapse = '-')
  if (pre_pad_NA != ''){
    pre_pad_NA <- paste0(pre_pad_NA, '-')
  }
  between_pad_NA <- rep(NA, j.prime - j - 1)
  between_pad_NA <- paste0(between_pad_NA, collapse = '-')
  if (between_pad_NA != ''){
    between_pad_NA <- paste0('-', between_pad_NA, '-')
  }else{
    between_pad_NA <- '-'
    
  }
  post_pad_NA <- rep(NA, J - j.prime)
  post_pad_NA <- paste0(post_pad_NA, collapse = '-')
  if (post_pad_NA != ''){
    post_pad_NA <- paste0('-', post_pad_NA)
  }
  if (design_method == 'column'){
    id <- paste(
      pre_pad_NA,
      data[,j],
      between_pad_NA,
      data[,j.prime],
      post_pad_NA,
      sep = ''
    )
  }else if (design_method == 'pairwise'){
    id <- paste(
      pre_pad_NA,
      data[,1],
      between_pad_NA,
      data[,2],
      post_pad_NA,
      sep = ''
    )
  }
  return(id)
}

# Trim a Design Matrix and Associated Penalty
# 
# Remove sufficiently rare entries from a design matrix by removing the column
# from the analysis.
remove_rare_levels <- function(X, penalty_for_regression, rare_threshold, verbose){
  fmt_rare <- rare <- NULL
  #Given a numeric rare_threshold, remove columns with fewer than rare_threshold
  #observations.
  if (is.vector(rare_threshold) & is.numeric(rare_threshold)){
    colsums_X <- colSums(X)  
    rare_columns <- which(colsums_X < rare_threshold)  
    # Only exclude rare *interactions*
    main_terms <- grep(colnames(X), pattern='-', invert = TRUE)
    rare_main <- rare_columns[rare_columns %in% main_terms]
    rare_columns <- rare_columns[!(rare_columns %in% main_terms)]
    
    if (verbose > 0 & (length(rare_main) > 0)){
      message('The following (main) factor levels are below rare_threshold. They are not excluded, but you may wish to do so or re-level the corresponding factor.', appendLF = F)
      if (verbose > 1){
        message(paste0(': ', paste(names(rare_main), collapse=', ')))    
      }
      message("\n", appendLF = F)
    }
    
  }else if (is.vector(rare_threshold) & is.character(rare_threshold)){
    
  }else{stop('rare_threshold should be number or character vector.')}  
  
  if (length(rare_columns) == 0){
    return(list(X = X, penalty_for_regression = penalty_for_regression, 
      rare = rare, fmt_rare= fmt_rare))
  }

  # Remove the rare columns from the design
  rare_X <- X[, -rare_columns, drop=F]
  # Remove the rare columns from the penalty, i.e.
  # from the list of terms
  penalty_for_regression$term_position <- penalty_for_regression$term_position[-rare_columns,,drop=F]
  penalty_for_regression$coef <- penalty_for_regression$coef[-rare_columns]
  # from the linear constraints
  penalty_for_regression$constraint <- penalty_for_regression$constraint[,-rare_columns,drop=F]
  # recalculate the basis given this reduced matrix
  penalty_for_regression$constraint_basis <- calculate_nullspace_basis(penalty_for_regression$constraint)
  # DROP from the penalty matrix
  penalty_for_regression$F <- lapply(penalty_for_regression$F, FUN=function(F_j){
    lapply(F_j, FUN=function(j){j <- j[-rare_columns,-rare_columns,drop=F]})
  }) 
  # DROP from the difference maker
  penalty_for_regression$D <- lapply(penalty_for_regression$D, FUN=function(D_j){
    lapply(D_j, FUN=function(j){j <- j[, -rare_columns, drop=F]})
  }) 
  
  fmt_rare <- strsplit(x = names(rare_columns), split='-')
  J_names <- penalty_for_regression$J_names
  
  fmt_rare <- sapply(fmt_rare, FUN=function(i){
    paste(paste0(J_names[which(i != 'NA')], '(', i[i != 'NA'], ')'), collapse = '-')
  })

  if (verbose > 0){
    message(paste0('Excluding ', length(fmt_rare), ' columns from design'), appendLF = F)
    if (verbose > 1){
      message(paste0(': ', paste(fmt_rare, collapse=', ')))    
    }
    message("\n", appendLF = F)
  }
  return(list(X = rare_X, penalty_for_regression = penalty_for_regression, rare = rare_columns, fmt_rare = fmt_rare))  
}

