
#' Create penalty matrix (list of F)
#' @keywords internal
#' @param factor_levels A list with the levels of each factor
#' @param treatment_probs A list with the levels of each factor and the probability of treatment assignment
#' @return A list of penalty matrices (F) used for sparsification
#' @import Matrix
#' @importFrom Rcpp sourceCpp
#' @useDynLib FactorHet
create_penalty <- function(factor_levels, make_interactions, 
    weight_dlist, ordered_factors, treatment_probs = NULL){
  #J is the number of factors
  J <- length(factor_levels)
  J_names <- names(factor_levels)
  
  if (is.null(J_names)){stop('Factor Levels cannot have NULL names')}
  if (anyDuplicated(J_names) != 0){
    stop('Factor Levels cannot have duplicate names')
  }
  
  if (!is.null(treatment_probs)){
    stop('Non-uniform treatment probabilities not yet implemented.')
  }else{
    #If null, assume evenly and fully randomized treatments.
    treatment_probs <- lapply(factor_levels, FUN=function(l){rep(1/length(l), length(l))})
  }
  
  #If creating interactions, do the following:
  if (is.logical(make_interactions)){
    if (make_interactions){
      make_interactions <- list_from_cols(combn(J, 2))
      make_interactions <- lapply(make_interactions, FUN=function(i){J_names[i]})
    }else{
      make_interactions <- list()
    }
  }else if (inherits(make_interactions, 'list')){
    valid_list <- all(lengths(make_interactions) == 2)
    if (!valid_list){stop('If providing list of interactions, must have two items in each.')}
    invalid_names <- setdiff(unlist(make_interactions), J_names)
    if (length(invalid_names) > 0){
      stop(paste('Some interaction names not found in factor list:', 
           paste(invalid_names, collapse=', ')))
    }
  }else{stop('Must give either T/F or list with two items')}
  
  all_interactions <- matrix(nrow = 0, ncol = J)
  factor_lengths <- lengths(factor_levels)
  
  #Create pairwise interactions between the relevant factors
  for (int in make_interactions){
    sort_int <- match(int, J_names)
    j <- sort_int[1]
    j.prime <- sort_int[2]
    int_matrix <- matrix(nrow = prod(factor_lengths[c(j, j.prime)]), ncol = J)
    int_matrix[,j] <- rep(factor_levels[[j]], each = factor_lengths[j.prime])
    int_matrix[,j.prime] <- rep(factor_levels[[j.prime]], times = factor_lengths[j])
    all_interactions <- rbind(all_interactions, int_matrix)
  }


  #Get the probability of assignment to treatment for each
  #pair of levels: List for independently randomized
  if (class(treatment_probs) != 'list'){
    stop('Non-independent arms not implemented.')
  }
  
  #Get the list of first order effects to put in front
  main_terms <- unlist(mapply(factor_levels, J_names, FUN=function(a,b){paste0(b, '(', a, ')')}))
  main_terms <- c('(Intercept)', main_terms)
  #Verify that no duplicates exist; shouldn't occur with factor_name(factor_level) phrasing
  stopifnot(anyDuplicated(main_terms) == 0)
  
  term_position <- rbind(matrix(NA, nrow = length(main_terms), ncol = J), 
    all_interactions)
  
  all_terms <- c(main_terms, apply(all_interactions, MARGIN = 1, paste, collapse = '-'))
  # all_terms <<- all_terms
  #Build the list of F matricies
  D_list <- F_list <- list()
  
  for (j in 1:J){
    levels_j <- factor_levels[[j]]
    D_j <- F_j <- list()
    #For each factor, loop over all levels.
    for (p_l in 1:length(levels_j)){
      # If levels are UNordered, fuse all pairwise combinations
      if (!ordered_factors[J_names[j]]){
        loop_range <- 1:p_l
      }else{
      # If levels are ordered, fuse adjacent categories.
        loop_range <- max(1, p_l-1):min(p_l, p_l+1)
      }
      for (p_m in loop_range){
        if (p_m != p_l){
          
          l_l <- levels_j[p_l]
          l_m <- levels_j[p_m]
          
          t_l <- which(l_l == term_position[,j])
          t_m <- which(l_m == term_position[,j])
          position_l <- term_position[t_l,,drop=F]
          position_m <- term_position[t_m,,drop=F]
          
          verify_match <- all.equal(position_l[,-j,drop=F], 
                                    position_m[,-j,drop=F])
          if (!verify_match){stop('Creating F matrix failed')}
          
          main_l <- match(paste0(J_names[j], '(', l_l, ')'), all_terms)
          main_m <- match(paste0(J_names[j], '(', l_m, ')'), all_terms)
          
          L_j <- length(levels_j)
          if (weight_dlist){
            size_weight <- length(t_l)
            if (size_weight == 0){size_weight <- 1}
            o_l <- c(length(main_l) * size_weight, rep(1, length(t_l)))
            o_m <- c(length(main_m) * size_weight, rep(1, length(t_m)))
          }else{
            o_l <- c(length(main_l) * 1, rep(1, length(t_l)))
            o_m <- c(length(main_m) * 1, rep(1, length(t_m)))
          }

          t_l <- c(main_l, t_l)
          t_m <- c(main_m, t_m)
          if (length(t_l) != length(t_m)){stop('Not aligned lengths')}
          
          Fmatrix <- rbind(
            cbind(t_l, t_l, o_l),
            cbind(t_m, t_m, o_l),
            cbind(t_m, t_l, -o_l),
            cbind(t_l, t_m, -o_l)
          )
          Fmatrix <- sparseMatrix(i = Fmatrix[,1], 
              j = Fmatrix[,2], x = Fmatrix[,3],
              dimnames = list(all_terms, all_terms),
              dims = rep(nrow(term_position), 2))
          F_j[[paste(l_m, l_l, sep = '-')]] <- Fmatrix
          
          DiffMatrix <- sparseMatrix(i = rep(1:length(t_l), 2), j = c(t_l, t_m), 
                                     x = c(rep(1, length(t_l)), rep(-1, length(t_m))),
                                     dims = c(length(t_l), nrow(term_position)))
          D_j[[paste(l_m, l_l, sep = '-')]] <- DiffMatrix
        }
      }
    }
    F_list[[j]] <- F_j
    D_list[[j]] <- D_j
  }
  names(D_list) <- names(factor_levels)
  names(F_list) <- names(factor_levels)
  #Verify that there are the right number of F
  #matrices in each factor level.
  
  stopifnot(all(mapply(F_list, factor_lengths, ordered_factors, FUN=function(Fi, fl, fo){
    if (fo){
      return(length(Fi) == (fl - 1))
    }else{
      return(length(Fi) == choose(fl, 2))
    }
  })))
  
  ##Create the constraint matrix M
  ###Create the constraints on the MAIN effects
  fmt_factor_levels <- mapply(factor_levels, J_names, SIMPLIFY = FALSE, FUN=function(a,b){paste0(b, '(', a, ')')})
  M_main <- mapply(fmt_factor_levels, treatment_probs, SIMPLIFY = FALSE, FUN=function(l, pr_l){
    sparseMatrix(i = rep(1, length(pr_l)), j = match(l, all_terms),  x = pr_l, 
                 dims = c(1,length(all_terms)))
  })
  M_main <- do.call('rbind', M_main)
  ###Create the constraints on the INTERACTIONS
  n_terms <- length(all_terms)
  M_inter <- sparseMatrix(i = 1, j = 1, x = 0, dims = c(1, n_terms))
  
  for (int in make_interactions){#Loop over all permitted interactions
    for (j in int){#loop "both" directions j and j.prime
      for (j.prime in int){
        if (j == j.prime){next} #Skip duplicated
        #Implement Equation (14) in Egami and Imai (2019)
        for (m in factor_levels[[j.prime]]){
          sumzero_position <- length(main_terms) + which(!is.na(all_interactions[,which(J_names == j)]) & 
                                                           (all_interactions[,which(J_names == j.prime)] == m))
          sumzero_weight <- treatment_probs[[j]]
          M_inter <- rbind(M_inter, sparseMatrix(i = rep(1, length(sumzero_position)), 
           j = sumzero_position, x = sumzero_weight,
           dims = c(1, n_terms)))
        }
      }
    }
  }
  
  M_inter <- M_inter[-1,,drop=F]
  dimnames(M_inter)[[2]] <- all_terms
  
  M <- rbind(M_main, M_inter)
    
  basis_M <- calculate_nullspace_basis(M)
  
  return(list("F" = F_list, "D" = D_list,
              "constraint" = M, 
              "constraint_basis" = basis_M,
              "coef" = all_terms, "J" = J, 
              "make_interactions" = make_interactions,
              'term_position' = term_position,
              'J_names' = J_names))
}


# Compare the difference between levels of factor and their interactions
prepare_fusion <- function(factor_levels, term_position, coef_names, beta, simplify = TRUE){
  F_list <- list()
  J <- length(factor_levels)
  J_names <- names(factor_levels)

  for (j in 1:J){
    levels_j <- factor_levels[[j]]
    F_j <- list()
    #For each factor, loop over all levels.
    for (p_l in 1:length(levels_j)){
      for (p_m in 1:p_l){
        if (p_m != p_l){
          
          l_l <- levels_j[p_l]
          l_m <- levels_j[p_m]
          
          t_l <- which(l_l == term_position[,j])
          t_m <- which(l_m == term_position[,j])
          position_l <- term_position[t_l,,drop=F]
          position_m <- term_position[t_m,,drop=F]
          
          position_l_negj <- apply(position_l[,-j,drop=F], MARGIN = 1, paste, collapse = ' @ ')
          position_m_negj <- apply(position_m[,-j,drop=F], MARGIN = 1, paste, collapse = ' @ ')
          
          # verify_match <- all.equal(position_l[,-j,drop=F],
          #                           position_m[,-j,drop=F])
          # if (!verify_match){stop('Creating F matrix failed')}
          # 
          # Get union of terms
          interaction_union <- base::union(position_l_negj, position_m_negj)
          # Find corresponding beta if it exists (not set to zero by fiat)
          inter_l <- beta[t_l[match(interaction_union, position_l_negj)]]
          inter_j <- beta[t_m[match(interaction_union, position_m_negj)]]
          # Set to zero
          inter_l[is.na(inter_l)] <- 0
          inter_j[is.na(inter_j)] <- 0
          diff_on_inter <- abs(inter_l - inter_j)
          
          main_l <- match(paste0(J_names[j], '(', l_l, ')'), coef_names)
          main_m <- match(paste0(J_names[j], '(', l_m, ')'), coef_names)
          
          diff_on_main <- abs(beta[main_l] - beta[main_m])
          if (length(diff_on_inter) == 0){diff_on_inter <- 0}
          if (any(is.na(diff_on_inter))){stop('')}
          F_j[[paste(l_m, l_l, sep = '-')]] <- list(main = diff_on_main, max_inter = max(diff_on_inter), all_inter = diff_on_inter)
        }
      }
    }
    F_list[[j]] <- F_j
  }
  names(F_list) <- names(factor_levels)
  if (simplify){
    
    F_simple <- mapply(F_list, names(F_list), SIMPLIFY = FALSE, FUN=function(F.j, j){
      output <- data.frame(t(sapply(F.j, FUN=function(l){cbind(l$main, l$max_inter)})))
      output$name <- names(F.j)
      output$factor <- j
      return(output)
    })
    F_simple <- do.call('rbind', F_simple)
    rownames(F_simple) <- NULL
    
    names(F_simple)[1:3] <- c('main', 'inter', 'pair')
    return(F_simple)
  }else{
    return(F_list)
  }
}

# Get fast way to calculate  binding positions
build_binding_lookup <- function(factor_levels, Fmatrix, term_position, coef_names){
  F_list <- list()
  J <- length(factor_levels)
  J_names <- names(factor_levels)
  
  for (j in 1:J){
    levels_j <- factor_levels[[j]]
    F_j <- list()
    #For each factor, loop over all levels.
    for (p_l in 1:length(levels_j)){
      for (p_m in 1:p_l){
        if (p_m != p_l){
          
          l_l <- levels_j[p_l]
          l_m <- levels_j[p_m]
          
          t_l <- which(l_l == term_position[,j])
          t_m <- which(l_m == term_position[,j])
          position_l <- term_position[t_l,,drop=F]
          position_m <- term_position[t_m,,drop=F]
          
          position_l_negj <- apply(position_l[,-j,drop=F], MARGIN = 1, paste, collapse = ' @ ')
          position_m_negj <- apply(position_m[,-j,drop=F], MARGIN = 1, paste, collapse = ' @ ')
          
          # verify_match <- all.equal(position_l[,-j,drop=F],
          #                           position_m[,-j,drop=F])
          # if (!verify_match){stop('Creating F matrix failed')}
          # 
          # Get union of terms
          interaction_union <- base::union(position_l_negj, position_m_negj)
          # Find corresponding beta if it exists (not set to zero by fiat)
          match_l <- t_l[match(interaction_union, position_l_negj)]
          match_m <- t_m[match(interaction_union, position_m_negj)]

          main_l <- match(paste0(J_names[j], '(', l_l, ')'), coef_names)
          main_m <- match(paste0(J_names[j], '(', l_m, ')'), coef_names)
          
          F_j[[paste(l_m, l_l, sep = '-')]] <- list(
            main = c(main_l, main_m),
            inter = cbind(match_l, match_m)
          )
            
        }
      }
    }
    F_list[[j]] <- F_j
  }
  #Verify integrity of name match  
  stopifnot(all(mapply(F_list, Fmatrix, FUN=function(i,j){all.equal(names(i), names(j))})))
  return(F_list)
}

calc_distance_binding <- function(b, lookup){
  lapply(lookup, FUN=function(j){
    sapply(j, FUN=function(lm){
      main_diff <- abs(diff(b[lm$main]))
      b_inter_l <- b[lm$inter[,1]]
      b_inter_m <- b[lm$inter[,2]]
      b_inter_l[is.na(b_inter_l)] <- 0
      b_inter_m[is.na(b_inter_m)] <- 0
      if (length(b_inter_l) > 0 | length(b_inter_m) > 0){
        inter_diff <- max(abs(b_inter_l - b_inter_m))
      }else{
        inter_diff <- 0
      }
      max_diff <- max(c(main_diff, inter_diff))
      return(max_diff)
    })
  })
}

standardization_weights <- function(X, D_list){

  dmat <- do.call('rbind', unlist(D_list))
  
  inv_dmat <- solve(Diagonal(n = ncol(dmat)) + crossprod(dmat))
  M_ginv <- inv_dmat %*% cbind(Diagonal(n = ncol(dmat)), t(dmat))  
  # Direct solution; sub-optimal compared to analytic inverse above
  # M_ginv <- drop0(zapsmall(ginv(as.matrix(rbind(Diagonal(n = ncol(dmat)), dmat)))))
  
  pos_id <- do.call('c', lapply(D_list, FUN=function(i){
    sapply(i, FUN=function(j){
      nrow(j)
    })
  }))
  cs_id <- cumsum(c(1, pos_id))
  
  M_ginv <- M_ginv[,-1:-ncol(X)]
  
  square_frob <- rep(NA, length(pos_id))
  for (k in 1:length(pos_id)){
    range_id <- (c(cs_id[k], cs_id[k] -1 + pos_id[k]))
    sub_X <- X %*% M_ginv[,seq(range_id[1], range_id[2])]
    square_frob[k] <- sum(sub_X^2)
  }
  
  frob_weights <- sqrt(square_frob) 
  frob_weights <- frob_weights/sqrt(nrow(X))
  
  frob_weights <- split(frob_weights, rep(1:length(D_list), sapply(D_list, length)))
  names(frob_weights) <- names(D_list)
  return(frob_weights)
}

create_standard_group <- function(D_list, weight_dlist){
  n_obs <- lapply(D_list, FUN=function(i){sapply(i, nrow)}) 
  n_obs <- sapply(n_obs, sum)  
  cum_obs <- cumsum(n_obs)
  cum_obs <- c(0, cum_obs[-length(cum_obs)])
  names(cum_obs) <- names(D_list)
  max_obs <- sum(n_obs)
  
  group_D <- mapply(cum_obs, D_list, SIMPLIFY = FALSE, FUN=function(ci, D_i){
    n_obs_i <- sapply(D_i, nrow) 
    names_i <- names(D_i)
    cum_obs_i <- cumsum(c(0, n_obs_i))
    cum_obs_i <- cum_obs_i[-length(cum_obs_i)] + ci
    names(cum_obs_i) <- names_i
    out_matrix <- mapply(cum_obs_i, D_i, SIMPLIFY = FALSE, FUN=function(ci, Di_l){
      ci_r <- ci + seq_len(nrow(Di_l))
      if (!weight_dlist){
        ci_x <- 1
      }else{
        ci_x <- rep(1, length(ci_r))
        ci_x[1] <- length(ci_r)
      }
      sparseMatrix(i = ci_r, j = ci_r, x = ci_x, dims = c(max_obs, max_obs))
    })
    return(out_matrix)
  })
  return(group_D)
}

create_log_penalty <- function(D_list, weight_dlist){
  n_obs <- lapply(D_list, FUN=function(i){sapply(i, nrow)}) 
  n_obs <- sapply(n_obs, sum)  
  cum_obs <- cumsum(n_obs)
  cum_obs <- c(0, cum_obs[-length(cum_obs)])
  names(cum_obs) <- names(D_list)
  max_obs <- sum(n_obs)
  
  group_D <- mapply(cum_obs, D_list, FUN=function(ci, D_i){
    n_obs_i <- sapply(D_i, nrow) 
    names_i <- names(D_i)
    cum_obs_i <- cumsum(c(0, n_obs_i))
    cum_obs_i <- cum_obs_i[-length(cum_obs_i)] + ci
    names(cum_obs_i) <- names_i
    out_matrix <- mapply(cum_obs_i, D_i, SIMPLIFY = FALSE, FUN=function(ci, Di_l){
      ci_r <- ci + seq_len(nrow(Di_l))
      if (!weight_dlist){
        ci_x <- 1
      }else{
        ci_x <- rep(1, length(ci_r))
        ci_x[1] <- length(ci_r)
      }
      sparseMatrix(i = ci_r, j = ci_r, x = ci_x, dims = c(max_obs, max_obs))
    })
    return(out_matrix)
  })
  return(group_D)
}
