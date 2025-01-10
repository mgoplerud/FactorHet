#' @importFrom stats aggregate
sum_mod_discrete <- function(object, type = c('bar', 'row', 'column')){
  if (!inherits(object, 'FactorHet')){
    stop('sum_mod_discrete can only be called on FactorHet object.')
  }
  type <- match.arg(type)
  phi <- coef(object, 'phi')
  
  K <- nrow(phi)
  
  if (K == 1){stop('No moderators when K = 1.')}
  
  if (is.null(object$posterior$posterior_predictive)){
    stop('No moderators provided.')
  }
  
  #Get moderator variables and their levels
  
  var_W <- attr(object$internal_parameters$W$args$terms, 'dataClasses')
  W_level <- object$internal_parameters$W$args$xlev
  var_W <- var_W[var_W %in% c("factor", "character")]
  
  if(length(var_W)==0){return(NULL)}
  
  # Get the original data for estimation  
  design <- object$internal_parameters$data$design
  # Get the unique groups in the estimation data
  group_id <- object$internal_parameters$group$unique_groups  
  group_name <- object$internal_parameters$formula$other_parameters$name_group
  
  if (length(group_name) == 0){
    unique_W <- design[,names(var_W), drop = F][group_id,, drop = F]
  }else{ 
    unique_W <- design[, names(var_W), drop = F][match(group_id, design[[group_name]]), , drop = F]
  }
  #Add survey weights
  surv_weights <- object$internal_parameters$weights$weights_W
  if (length(surv_weights) != nrow(unique_W)){
    stop('Selecting of unique individuals failed. Try doing manually.')
  }
  # Get the data by moderator for each person
  flat_data <- do.call('rbind', lapply(names(var_W), FUN=function(w){
    value_w <- unique_W[[w]]
    out <- data.frame(
      id = group_id,
      survey_weight = surv_weights,
      value = value_w,
      variable = w, stringsAsFactors = F)
    return(out)
  }))
  # Add posterior predictive probabilities
  n_merge <- nrow(flat_data)
  flat_data <- merge(flat_data, object$posterior$posterior_predictive, by.x = 'id', by.y = 'group', sort = FALSE)
  if (n_merge != nrow(flat_data)){
    stop('Merging failed for sum_mod_continuous. Try doing manually using information from "posterior".')
  }
  group_used <- setdiff(names(object$posterior$posterior_predictive), 'group')
  flat_data <- do.call('rbind', lapply(1:K, FUN=function(k){
    flat_k <- flat_data[, c(1:4, match(paste0('group_', k), names(flat_data))), drop = F]
    names(flat_k)[5] <- 'group_prob'
    flat_k$group <- k
    return(flat_k)
  }))
  if ( (n_merge * K) != nrow(flat_data)){
    stop('Merging failed for sum_mod_continuous. Try doing manually using information from "posterior".')
  }
  # Get the weight for the plot
  flat_data$plot_weight <- flat_data$survey_weight * flat_data$group_prob
  flat_data$group <- factor(paste0('Group ', flat_data$group), levels = paste0('Group ', 1:K))
  
  if (type == 'row'){
    fmt_data <- split(flat_data[, c('plot_weight', 'group', 'value', 'variable')], 
                      flat_data[, c('variable', 'value')])
    fmt_data <- fmt_data[which(sapply(fmt_data, nrow) != 0)]
    fmt_data <- do.call('rbind', lapply(fmt_data, FUN=function(i){
      agg_i <- aggregate(plot_weight ~ group, data = i, FUN = sum)
      agg_i$value <- unique(i$value)
      agg_i$variable <- unique(i$variable)
      agg_i$norm_weight <- agg_i$plot_weight/sum(agg_i$plot_weight)
      return(agg_i)
    }))    
  }else if (type %in% c('column', 'bar')){
    fmt_data <- split(flat_data[, c('plot_weight', 'group', 'value', 'variable')], 
                      flat_data[, c('group', 'variable')])
    fmt_data <- fmt_data[which(sapply(fmt_data, nrow) != 0)]
    fmt_data <- do.call('rbind', lapply(fmt_data, FUN=function(i){
      agg_i <- aggregate(plot_weight ~ value, data = i, FUN = sum)
      agg_i$group <- unique(i$group)
      agg_i$variable <- unique(i$variable)
      agg_i$norm_weight <- agg_i$plot_weight/sum(agg_i$plot_weight)
      return(agg_i)
    }))
  }else{stop('invalid type')}
  rownames(fmt_data) <- NULL
  
  fmt_data$variable <- as.factor(fmt_data$variable)
  fmt_data$group <- factor(fmt_data$group, levels = paste0('Group ', 1:K))
  
  .data <- NULL
  if (type == 'bar'){
    g <- ggplot(data=fmt_data,
        aes(y=.data[['norm_weight']], 
            x=.data[['value']], 
            fill=.data[['group']])) +
      geom_bar(stat="identity", position=position_dodge())+
      xlab("Category") +
      labs(fill="Group number")+
      facet_grid(. ~ variable, scales="free") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  }else{
    fmt_data$round_norm <- round(fmt_data$norm_weight,2)
    fmt_data$round_norm <- sprintf('%.2f', fmt_data$round_norm)
    g <- ggplot(data=fmt_data, 
        aes(
          x = .data[['group']], 
          y = .data[['value']], 
          fill= .data[['norm_weight']], 
          label=.data[['round_norm']])) +
      geom_tile()+
      geom_text(color="white")+
      xlab("Group") +
      labs(fill="Percent")+
      scale_x_discrete(position = "top") +
      scale_fill_gradient(low="#56B1F7", high="#132B43") +
      theme_bw() + 
      theme(strip.text.y.left = element_text(angle = 0), 
            panel.spacing = unit(0.1, 'lines'), strip.placement = 'outside') +
      facet_grid(variable~., switch = 'y', scales="free_y") +
      ylab('')
  }
  out <- list(plot = g, data = fmt_data)
  class(out) <- 'FactorHet_vis'
  return(out)
}

sum_mod_continuous <-function(object){
  
  if (!inherits(object, 'FactorHet')){
    stop('sum_mod_continuous can only be called on FactorHet object.')
  }
  
  phi <- coef(object, 'phi')
  K <- nrow(phi)
  
  if (K == 1){stop('No moderators when K = 1.')}
  
  if (is.null(object$posterior$posterior_predictive)){
    stop('No moderators provided.')
  }
  
  #Get moderator variables and their levels
  var_W <- attr(object$internal_parameters$W$args$terms, 'dataClasses')
  W_level <- object$internal_parameters$W$args$xlev
  var_W<-var_W[var_W=="numeric"]
  
  if(length(var_W)==0){
    return(NULL)
  }
  
  # Get the original data for estimation  
  design <- object$internal_parameters$data$design
  # Get the unique groups in the estimation data
  group_id <- object$internal_parameters$group$unique_groups  
  group_name <- object$internal_parameters$formula$other_parameters$name_group
  
  if (length(group_name) == 0){
    unique_W <- design[,names(var_W), drop = F][group_id,, drop = F]
  }else{ 
    unique_W <- design[, names(var_W), drop = F][match(group_id, design[[group_name]]), , drop = F]
  }
  #Add survey weights
  surv_weights <- object$internal_parameters$weights$weights_W
  if (length(surv_weights) != nrow(unique_W)){
    stop('Selecting of unique individuals failed. Try doing manually.')
  }
  # Get the data by moderator for each person
  flat_data <- do.call('rbind', lapply(names(var_W), FUN=function(w){
    data.frame(
      id = group_id,
      survey_weight = surv_weights,
      value = unique_W[[w]],
      variable = w, stringsAsFactors = F)
  }))
  # Add posterior predictive probabilities
  n_merge <- nrow(flat_data)
  flat_data <- merge(flat_data, object$posterior$posterior_predictive, by.x = 'id', by.y = 'group')
  if (n_merge != nrow(flat_data)){
    stop('Merging failed for sum_mod_continuous. Try doing manually using information from "posterior".')
  }
  group_used <- setdiff(names(object$posterior$posterior_predictive), 'group')
  flat_data <- do.call('rbind', lapply(1:K, FUN=function(k){
    flat_k <- flat_data[, c(1:4, match(paste0('group_', k), names(flat_data))), drop = F]
    names(flat_k)[5] <- 'group_prob'
    flat_k$group <- k
    return(flat_k)
  }))
  if ( (n_merge * K) != nrow(flat_data)){
    stop('Merging failed for sum_mod_continuous. Try doing manually using information from "posterior".')
  }
  # Get the weight for the plot
  flat_data$plot_weight <- flat_data$survey_weight * flat_data$group_prob
  flat_data$group <- factor(paste0('Group ', flat_data$group), levels = paste0('Group ', 1:K))
  #Plot
  .data <- NULL
  g <- ggplot(data=flat_data, 
      aes(x=.data[['group']], 
          y=.data[['value']], 
          weight=.data[['plot_weight']])) +
      geom_boxplot()+
      labs(title="Weighted boxplots by group", x ="Group", y = "")+
      scale_x_discrete(position = "top") +
      facet_grid(variable ~., scales="free_y", switch = 'y') +
      theme_bw() +
      theme(strip.text.y.left = element_text(angle = 0), 
        panel.spacing = unit(0.1, 'lines'), strip.placement = 'outside') 
    
  out <- list(plot = g, data = flat_data)
  class(out) <- 'FactorHet_vis'
  return(out)
}

#' Visualize the posterior by observed moderators
#'
#' Provides univariate summaries of the estimated posterior predictive
#' probabilities of group membership by the moderators. Can produce analyses
#' for continuous variables (weighted boxplot) or discrete variables (row/column
#' tables).
#' 
#' @return A list of each of the types of analyses is reported. Each element of
#'   the list contains the ggplot object and the data ("plot" and "data").
#'  
#' @details  
#'   \bold{Discrete Moderators}: Discrete moderators are shown by either a
#'   \code{"row"}, \code{"column"}, or \code{"bar"} plot. In the \code{"row"}
#'   plot, the quantity reported is, for each level of the moderator, what
#'   proportion of people fall into each group. For example, for moderator value
#'   "a", 25\% of people are in group 1 and 75\% of people are in group 2.
#'   This is estimated using a weighted average, weighting by the estimated posterior
#'   predictive probabilities of group membership and any survey weights.
#'   
#'   By contrast \code{"column"} and \code{"bar"} reports the distribution
#'   by group. For example, for Group 1, 30\% of people have moderator value "f",
#'   50\% have moderator value "g", and 20\% have moderator value "h".
#'   \code{"bar"} reports this as a bar chart whereas \code{"column"} reports as a
#'   tile plot.
#'   
#'   For all three types of plots, the data is provided in the returned output.
#'   
#'   \bold{Continuous Moderators}: Continuous moderators are shown by a
#'   histogram of the value for each group, weighted by each observation's
#'   posterior predictive probability of being in that group.
#' @param object A model fit using \code{\link{FactorHet}} or \code{\link{FactorHet_mbo}}.
#' @param visualize Specifies which types of moderators to show. Default (\code{"all"})
#'   shows all moderators. Other options include \code{"discrete"} and
#'   \code{"continuous"}.
#' @param type_discrete Show the results by \code{"row"} or \code{"column"} or
#'   \code{"all"} (i.e. both).
#' 
#' @examples 
#' data(immigration)
#' set.seed(15)
#' # Estimate model with arbitrary choice of lambda
#' fit <- FactorHet(Chosen_Immigrant ~ Plans + Ed + Country,
#'   design = immigration, lambda = 1e-2,
#'   moderator = ~ party_ID,
#'   K = 2, group = ~ CaseID,
#'   control = FactorHet_control(init = 'mclust'),
#'   task = ~ contest_no, choice_order = ~ choice_id)
#' posterior_by_moderators(fit)
#' @importFrom stats quantile
#' @export
posterior_by_moderators <- function(object, 
  visualize = c('all', 'discrete', 'continuous'),
  type_discrete = c('bar', 'row', 'column', 'all')){
  
  type_discrete <- match.arg(type_discrete)
  visualize <- match.arg(visualize)
  
  plot_discrete_row <- plot_discrete_col <- plot_continuous <- NULL
  
  if (visualize %in% c('all', 'discrete')){
    if (type_discrete %in% c('all', 'bar')){
      plot_discrete_bar <- sum_mod_discrete(object = object, type = 'bar')
    }
    if (type_discrete %in% c('all', 'row')){
      plot_discrete_row <- sum_mod_discrete(object = object, type = 'row')
    }
    if (type_discrete %in% c('all', 'column')){
      plot_discrete_col <- sum_mod_discrete(object = object, type = 'column')
    }
  }
  if (visualize %in% c('all', 'continuous')){
    plot_continuous <- sum_mod_continuous(object = object)
  }
  
  out <- list(
    'discrete_bar' = plot_discrete_bar,
    'discrete_column' = plot_discrete_col,
    'discrete_row' = plot_discrete_row,
    'continuous' = plot_continuous
  )
  out <- out[!sapply(out, is.null)]
  if (length(out) == 1){
    out <- out[[1]]
  }
  return(out)
}

#Create table of joint posterior predictive prob of two groupings
#data_obj is a data frame with the posterior pred prob for each individual for each grouping
#clus_var_nam1 is a vector of the column names corresponding to the groups in the first grouping
#clus_var_nam2 is the same for the second grouping
group_table_fun<-function(data_obj, clus_var_nam1, clus_var_nam2){
  table.dat <- data.frame(matrix(ncol = length(clus_var_nam1), nrow = length(clus_var_nam2)), row.names=clus_var_nam2)
  colnames(table.dat) <- clus_var_nam1
  for(i in clus_var_nam1){
    for(j in clus_var_nam2){
      table.dat[j,i]<-mean(data_obj[,i]*data_obj[,j])
    }
  }
  return(table.dat)
}


#' Compute association between moderators and group membership
#' 
#' @description This function computes the impact of changing a
#'   moderator on the group membership probabilities.
#' @details This function computes the change in \eqn{\pi_k(X_i)} for the change
#'   in one of the moderators in \eqn{X_i}. The change is averaged across the
#'   distribution of the other moderators found in \code{newdata} (or, by
#'   default, the estimation data). It thus can be thought of as the "marginal
#'   effect" of changing one moderator on the probability of group memberships,
#'   holding all other moderators constant. It returns a data.frame of the
#'   estimated effects as well as a plot to visualize the changes in
#'   \eqn{\pi_k(X_i)}. Goplerud et al. (2025) provides more discussion of this
#'   method.
#' 
#' @param object An object from \code{\link{FactorHet}} or
#'   \code{\link{FactorHet_mbo}}.
#' @param newdata An optional argument that provides the data over which to
#'   average the distribution of the other moderators. The default is
#'   \code{NULL} which uses the estimation data.
#' @param vcov A logical value indicating whether the standard errors should be
#'   computed. The default is \code{TRUE}.
#' @param se.method An optional argument as to the type of standard errors used.
#'   The default is \code{NULL} uses estimated standard errors.
#'   \code{\link{vcov.FactorHet}} provides more information.
#' @param quant_continuous A numeric vector consisting of two values between 0
#'   and 1. For continuous moderators, it sets two quantiles of the moderator's
#'   distribution to show the difference between. The default \code{c(0.25,
#'   0.75)} compares the effect of changing the moderator from its 25th
#'   percentile to its 75th percentile.
#' @param abs_diff A logical value as to whether the difference or absolute
#'   difference in the change in \eqn{\pi_k(X_i)} should be shown. The default
#'   is \code{FALSE} which returns the standard "marginal effect" of changing
#'   the moderators with a standard error computed via the delta method. The
#'   value \code{TRUE} draws 10,000 samples from the asymptotic distribution of
#'   the moderators and computes the average the \bold{absolute values} of the
#'   marginal effects for each observation in \code{newdata} using those
#'   samples. This is considerably slower than the default setting. The appendix
#'   of Goplerud et al. (2025) illustrates one use of this argument.
#'  
#' @examples
#' # Estimate model with arbitrary choice of lambda
#' data(immigration)
#' set.seed(15)
#' # Estimate model with arbitrary choice of lambda
#' fit <- FactorHet(Chosen_Immigrant ~ Plans + Ed + Country,
#'   design = immigration, lambda = 1e-2,
#'   moderator = ~ party_ID,
#'   K = 2, group = ~ CaseID,
#'   control = FactorHet_control(init = 'mclust'),
#'   task = ~ contest_no, choice_order = ~ choice_id)
#' margeff_moderators(fit)
#' @return Returns a named list with the underlying data (\code{"data"}) and the
#'   plot (\code{"plot"}).
#' @import ggplot2
#' @export
margeff_moderators <- function(object, newdata = NULL, vcov = TRUE,
                          se.method = NULL, 
                          quant_continuous = c(0.25, 0.75),
                          abs_diff = FALSE){
  
  if (!inherits(object, 'FactorHet')){
    stop('object must be from FactorHet')
  }
  
  phi <- coef(object, 'phi')
  K <- nrow(phi)
  if (K == 1){stop('No moderators when K = 1.')}
  
  if (is.null(object$posterior$posterior_predictive)){
    stop('No moderators provided.')
  }
  
  if (vcov){
    object_vcov <- vcov.FactorHet(object, phi = TRUE, se.method = se.method)
    if (is.null(object_vcov)){
      vcov <- FALSE
    }else{
      n_p <- nrow(coef(object))
      object_vcov <- object_vcov[-seq_len(n_p * K), -seq_len(n_p * K)]
    }
  }else{
    object_vcov <- NULL
  }
  
  if (length(quant_continuous) != 2){stop('Provide two quantiles to compare.')}
  
  #Get the "base" variables used in moderators
  var_W <- attr(object$internal_parameters$W$args$terms, 'dataClasses')
  W_level <- object$internal_parameters$W$args$xlev
  #For each variable, do the prediction
  all_mfx <- data.frame()
  
  if (abs_diff){
    if (!vcov){stop('Cannot do "abs_diff" if vcov=FALSE')}
    # Draw samples from N(mu, V)    
    # L %*% t(L) = V
    # Cholesky decompose the variance
    chol_vcov <- t(chol(object_vcov))
    stopifnot(isTRUE(all.equal(as.matrix(chol_vcov %*% t(chol_vcov)) , object_vcov)))
    # Draws samples from N(0, V)
    sim_phi <- matrix(rnorm(10^4 * ncol(object_vcov)), ncol = ncol(object_vcov)) %*% t(chol_vcov)
    # Add in mean mu to each sample
    sim_phi <- sweep(sim_phi, MARGIN = 2, FUN = '+', STATS = as.vector(t(coef(object, 'phi')[-1,])))
  }
  
  for (v in names(var_W)){
    
    copy_design <- object$internal_parameters$data$design
    
    #Get the levels to predict over
    if (var_W[v] %in% c('numeric', 'matrix')){
      unique_levels_v <- length(unique(copy_design[[v]]))
      if (unique_levels_v > 2){
        unique_levels_v <- quantile(copy_design[[v]], quant_continuous)
        var_type <- 'continuous'
      }else{
        if (unique_levels_v != 2){stop(paste0(v, ' is constant?'))}
        unique_levels_v <- sort(unique(copy_design[[v]]))
        var_type <- 'binary'
      }
    }else{
      unique_levels_v <- W_level[[v]]
      baseline <- unique_levels_v[1]
      unique_levels_v <- unique_levels_v[-1]
      if (length(unique_levels_v) == 0){stop(paste0(v, ' has only one level?'))}
      var_type <- 'factor'
    }
    
    if (!is.null(newdata)){
      copy_design <- newdata
    }
    
    
    if (var_type == 'factor'){
      copy_design[,v] <- baseline
      pred_baseline <- predict(object, newdata = copy_design, calc_gradient = vcov, 
                               override_weights = FALSE,
                               return = 'postpred_only')
      if (vcov){
        W_baseline <- attributes(pred_baseline)$W
      }
      norm_weights <- attributes(pred_baseline)$norm_weights
      
      mean_baseline <- colSums(Diagonal(x = norm_weights) %*% pred_baseline)
      
      for (u in unique_levels_v){
        copy_design[,v] <- u
        pred_u <- predict(object, newdata = copy_design, 
          calc_gradient = vcov, 
          return = 'postpred_only', override_weights = FALSE)
        
        stopifnot(identical(attributes(pred_u)$norm_weights,
                            norm_weights))
        
        if (vcov){
          W_u <- attributes(pred_u)$W
          grad_diff <- delta_ME_multinom(K = K, prob_high = pred_u, W_high = W_u,
             prob_low = pred_baseline, W_low = W_baseline, vcov = object_vcov,
             weights = norm_weights)              
        }else{
          grad_diff <- rep(NA, K)
        }
        
        mean_u <- colSums(Diagonal(x = norm_weights) %*% pred_u)
        #Average Posterior Predictive After Change
        pred_u <- colSums(Diagonal(x = norm_weights) %*% pred_u)
        
        if (abs_diff){
          # Simulate the values
          pred_sim <- do.call('rbind', lapply(1:nrow(sim_phi), FUN=function(i){
            sphi <- rbind(0, matrix(sim_phi[i,,drop=FALSE], byrow = TRUE, nrow = K - 1))
            diffphi <- softmax_matrix(W_u %*% t(sphi)) - softmax_matrix(W_baseline %*% t(sphi))
            data.frame(
              abs_diff = colSums(Diagonal(x=norm_weights) %*% abs(diffphi)),
              K = 1:K
            )
          }))
          out_mfx <- do.call('rbind', lapply(c('abs_diff'), FUN=function(v){
            out_v <- sapply(split(pred_sim[[v]], pred_sim$K), FUN=function(i){
              c('mean' = mean(i), 'median' = median(i), quantile(i, c(0.025, 0.975)))
            })
            out_v <- data.frame(t(out_v))
            names(out_v) <- c('mean', 'median', 'll', 'ul')
            out_v$K <- 1:nrow(out_v)
            out_v$qoi <- v
            return(out_v)
          }))
          out_mfx$point_estimate_average <- mean_u - mean_baseline
          out_mfx$variable <- paste0(v, '(', u, ')')
          out_mfx$type <- var_type
          all_mfx <- rbind(all_mfx, out_mfx)
        }else{
          all_mfx <- rbind(all_mfx, data.frame(variable = paste0(v, '(', u, ')'), 
           type = var_type, 'diff' = t(mean_u - mean_baseline), 
           'var' = t(grad_diff),
           'changed' = t(mean_u), 'baseline' = t(mean_baseline), 
           stringsAsFactors = F))
          
        }
        
      }
    }else{
      
      copy_design[,v] <- unique_levels_v[1]
      predict_low <- predict(object, newdata = copy_design, calc_gradient = vcov, return = 'postpred_only', override_weights = FALSE)
      copy_design[,v] <- unique_levels_v[2]
      predict_high <- predict(object, newdata = copy_design, calc_gradient = vcov, return = 'postpred_only', override_weights = FALSE)
      
      #check no issues with weights
      stopifnot(identical(attributes(predict_high)$norm_weights,
                          attributes(predict_low)$norm_weights))
      norm_weights <- attributes(predict_high)$norm_weights

      if (vcov){

        W_low <- attributes(predict_low)$W
        W_high <- attributes(predict_high)$W
        
        grad_diff <- delta_ME_multinom(K = K, prob_high = predict_high, W_high = W_high,
                                       prob_low = predict_low, W_low = W_low, vcov = object_vcov,
                                       weights = norm_weights)        
      }else{
        grad_diff <- rep(NA, K)
      }
      
      predict_low <- colSums(Diagonal(x = norm_weights) %*% predict_low)
      predict_high <- colSums(Diagonal(x = norm_weights) %*% predict_high)
      
      if (abs_diff){
        # Simulate the average of the absolute values
        pred_sim <- do.call('rbind', lapply(1:nrow(sim_phi), FUN=function(i){
          sphi <- rbind(0, matrix(sim_phi[i,,drop=FALSE], byrow = TRUE, nrow = K-1))
          diffphi <- softmax_matrix(W_high %*% t(sphi)) - softmax_matrix(W_low %*% t(sphi))
          data.frame(
                     abs_diff = colSums(Diagonal(x=norm_weights) %*% abs(diffphi)),
                     K = 1:K
          )
        }))
        out_mfx <- do.call('rbind', lapply(c('abs_diff'), FUN=function(v){
          out_v <- sapply(split(pred_sim[[v]], pred_sim$K), FUN=function(i){
            c('mean' = mean(i), 'median' = median(i), quantile(i, c(0.025, 0.975)))
          })
          out_v <- data.frame(t(out_v))
          names(out_v) <- c('mean', 'median', 'll', 'ul')
          out_v$K <- 1:nrow(out_v)
          out_v$qoi <- v
          return(out_v)
        }))
        out_mfx$variable <- v
        out_mfx$type <- var_type
        out_mfx$point_estimate_average <- predict_high - predict_low
        all_mfx <- rbind(all_mfx, out_mfx)
        # Average absolute Change in Posterior Predictive
        # pred_absolute <- colSums(Diagonal(x=norm_weights) %*% abs(predict_high - predict_low))
        # all_mfx <- rbind(all_mfx, data.frame(variable = v, type = var_type, 
        #                                      'diff' = t(pred_absolute), var = t(grad_diff),
        #                                      'signed_change' = t(colSums(Diagonal(x=norm_weights) %*% sign(predict_high - predict_low))),
        #                                      'changed' = NA, 'baseline' = NA,  stringsAsFactors = F))
      }else{
        #Average Change in Posterior Predictive
        all_mfx <- rbind(all_mfx, data.frame(variable = v, type = var_type, 
                                             'diff' = t(predict_high - predict_low), var = t(grad_diff),
                                             'changed' = t(predict_high), 'baseline' = t(predict_low),  stringsAsFactors = F))
      }
      
    }
  }  
  mfx_order <- as.vector(unlist(mapply(names(W_level), W_level, 
                                       FUN=function(i,j){paste0(i,'(', j, ')')})))
  
  mfx_order <- c(mfx_order, names(var_W[var_W == 'numeric']))

  if (abs_diff){
    # Already in long format...
  }else{
    if (vcov){
      for (k in 1:K){
        t_stat <- all_mfx[[paste0('diff.', k)]] / sqrt(all_mfx[[paste0('var.', k)]])
        sig <- abs(t_stat) > 1.96
        all_mfx[[paste0('t.', k)]] <- t_stat
        all_mfx[[paste0('sig.', k)]] <- as.numeric(sig)
      }
    }else{
      for (k in 1:K){
        all_mfx[[paste0('t.', k)]] <- NA
        all_mfx[[paste0('sig.', k)]] <- NA
      }
    }
    # Old version that uses reshape2
    # vis_mfx <- melt(all_mfx, id.vars = c('type', 'variable'),
    #   measure.vars = grep(names(all_mfx), pattern='^(t$|sig|changed|baseline)'),
    #   variable.name = 'mfx_type')
    # vis_mfx$group <- gsub(vis_mfx$mfx_type, pattern='[^0-9]+', perl = T, replacement = 'Group ')
    # vis_mfx$mfx_type <- gsub(vis_mfx$mfx_type, pattern='[\\.0-9]+', replacement ='')
    # vis_mfx <- dcast(vis_mfx, group + variable ~ mfx_type, value.var = 'value')
    vis_mfx <- do.call('rbind', lapply(grep(names(all_mfx), pattern='^(t$|sig|changed|baseline)'), FUN=function(i){
      data.frame(type = all_mfx$type,
                 variable = all_mfx$variable,
                 mfx_type = names(all_mfx)[i],
                 value = all_mfx[,i],
                 stringsAsFactors = FALSE)
    }))
    vis_mfx$group <- gsub(vis_mfx$mfx_type, pattern='[^0-9]+', perl = T, replacement = 'Group ')
    vis_mfx$mfx_type <- gsub(vis_mfx$mfx_type, pattern='[\\.0-9]+', replacement ='')
    joint_id <- paste(vis_mfx$group, '@@@@', vis_mfx$variable)
    int_vis <- split(vis_mfx[, c('value', 'mfx_type')], joint_id)
    vis_mfx <- do.call('rbind', lapply(int_vis, FUN=function(i){
      dat <- data.frame(x = t(i[,1]))
      names(dat) <- i[,2]
      return(dat)
    }))
    vis_mfx[,c('group', 'variable')] <- do.call('rbind', strsplit(names(int_vis), ' @@@@ '))
    rownames(vis_mfx) <- NULL
    vis_mfx <- vis_mfx[, c('group', 'variable', 'baseline', 'changed', 'sig')]
    vis_mfx$fmt_name <- factor(vis_mfx$variable, levels = mfx_order)
    vis_mfx$sig <- factor(vis_mfx$sig)
  }
  
  .data <- NULL
  if (abs_diff){
    all_mfx$fmt_name <- factor(all_mfx$variable, levels = mfx_order)
    vis_mfx <- all_mfx[all_mfx$qoi == 'abs_diff',]
    vis_mfx$group <- paste0('Group ', vis_mfx$K)
    vis_mfx$abs_orig <- abs(vis_mfx$point_estimate_average)
    g <- ggplot(vis_mfx) + 
      geom_point(aes(x=.data[['fmt_name']],y=.data[['mean']])) +
      geom_errorbar(aes(x=.data[['fmt_name']],
                               ymin=.data[['ll']],
                               ymax=.data[['ul']])) +
      coord_flip() + facet_wrap(~group) +
      theme_bw() + ylab('Change in Posterior Predictive Probability of Group Membership') + 
      xlab('Covariate')  +
      geom_hline(aes(yintercept=0), linetype='dashed') +
      geom_point(aes(x=.data[['fmt_name']], y=.data[['abs_orig']]), col = 'red', pch = 8)
  }else{
    g <- ggplot(vis_mfx, aes(alpha = .data[['sig']])) + 
      geom_point(aes(x=.data[['fmt_name']],y=.data[['baseline']])) +
      geom_segment(aes(x=.data[['fmt_name']],
                       xend=.data[['variable']],
                       y=.data[['baseline']],
                       yend=.data[['changed']]),
                   arrow = arrow(length =unit(0.03, 'npc'))) + 
      coord_flip() + facet_wrap(~group) +
      theme_bw() + ylab('Posterior Predictive Probability of Group Membership') + 
      xlab('Covariate') +
      scale_alpha_manual(values = c(0.25, 1), guide = 'none')
  }
  output <- list(plot = g, data = all_mfx)
  class(output) <- 'FactorHet_vis'
  return(output)
}

#' @rdname deprecated
#' @keywords internal
#' @export
moderator_AME <- function(...){
  .Deprecated(new = 'margeff_moderators', old = 'moderator_AME')
  return(margeff_moderators(...))
}

delta_ME_multinom <- function(K, 
  prob_high, prob_low, W_high, W_low, 
  vcov, weights, absolute = FALSE){
  
  delta_out <- sapply(1:K, FUN=function(k){
    if (absolute){
      grad_delta <- do.call('c', lapply(2:K, FUN=function(l){
        g_high <- (Diagonal(x = weights) %*% Diagonal(x = prob_high[,k] * ((l == k) - prob_high[,l])) %*% W_high)
        g_low <- (Diagonal(x = weights) %*% Diagonal(x = prob_low[,k] * ((l == k) - prob_low[,l])) %*% W_low)
        g_diff <- colSums(Diagonal(x=sign(prob_high[,k] - prob_low[,k])) %*% (g_high - g_low))
        return(g_diff)
      }))
      
    }else{
      grad_delta <- do.call('c', lapply(2:K, FUN=function(l){
        g_high <- colSums(Diagonal(x = weights) %*% Diagonal(x = prob_high[,k] * ((l == k) - prob_high[,l])) %*% W_high)
        g_low <- colSums(Diagonal(x = weights) %*% Diagonal(x = prob_low[,k] * ((l == k) - prob_low[,l])) %*% W_low)
        return(g_high - g_low)
      }))
    }
    delta_var <- as.numeric(t(grad_delta) %*% vcov %*% grad_delta)
    return(delta_var)
  })
  return(delta_out)
}
