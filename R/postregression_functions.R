extract_coefficient <- function(object, vcov.type = NULL, baseline = NA, verbose = TRUE){
  
  object_factor_levels <- object$internal_parameters$factor_levels
  J <- length(object_factor_levels)
  beta <- object$parameters$beta
  K <- ncol(beta)
  if (is.null(vcov.type)){
    vcov_object <- vcov.FactorHet(object, phi = FALSE)
  }else{
    vcov_object <- vcov.FactorHet(object, se.method = vcov.type, phi = FALSE)
  }
  vcov_object <- as.matrix(vcov_object)

  if (is.null(baseline)){
    baseline <- sapply(object_factor_levels, FUN=function(i){i[1]})
    use_baseline <- TRUE
    if (verbose){
      message('baseline = NULL gives the following baseline')
      print(unlist(baseline))
    }
  }else{
    if (is.na(baseline[1])){
      if (verbose){
        message('baseline = NA implies showing coefficients. Interpret by looking at differences. Set "null" to choose baseline automatically.')
      }
      use_baseline <- FALSE
    }else{
      if (length(baseline) != J){stop('Must provide J levels as baseline')}
      
      verify_baseline <- all.equal(sort(names(baseline)), sort(names(object_factor_levels)))
      if (!isTRUE(verify_baseline)){
        stop('baseline must be a named list with one entry for each factor_level. Some factors in "baseline" are not in the estimated model.')
      }
      baseline <- baseline[match(names(object_factor_levels), names(baseline))]
      
      verify_baseline <- all(mapply(baseline, object_factor_levels, FUN=function(b, l){b %in% l}))
      if (!verify_baseline){stop('All baseline levels must be in the associated factor.')}
      use_baseline <- TRUE      
    }
  }  
  
  out <- data.frame()
  
  for (j in 1:J){
    
    main_effect <- paste0(names(object_factor_levels)[j], '(', object_factor_levels[[j]], ')')
    coef_AME <- beta[match(main_effect, rownames(beta)),,drop=F]
    
    if (use_baseline){
      which_baseline <- match(baseline[j], object_factor_levels[[j]])
      
      baseline_AME <- beta[match(main_effect[which_baseline], rownames(beta)),,drop=F]
      
      est_AME <- sapply(1:K, FUN=function(k){
        coef_AME[,k] - baseline_AME[k] 
      })
      
      var_AME <- sapply(1:K, FUN=function(k){
        sapply(rownames(est_AME), FUN=function(j){
          fmt_names <- paste0('beta', k, '_', c(j, rownames(baseline_AME)))
          c(1,-1) %*% vcov_object[fmt_names, fmt_names] %*% c(1,-1)
        })
      })
      
    }else{
      est_AME <- as.matrix(coef_AME)
      
      var_AME <- sapply(1:K, 
             FUN=function(i){diag(vcov_object)[
               paste0('beta', i, '_', rownames(est_AME))]})
      
    }

    stopifnot(all(!is.na(var_AME)))  
    
    AME_flat <- data.frame(variable = rep(rownames(est_AME), K),
     group = rep(1:K, each = nrow(est_AME)),
     coef = as.vector(as.matrix(est_AME)), 
     var = as.vector(var_AME),
     id = rep(1:nrow(est_AME), K),
     stringsAsFactors = F)
    
    
    AME_flat$level <- object_factor_levels[[j]][AME_flat$id]
    AME_flat$factor <- names(object_factor_levels)[j]
    stopifnot(all.equal(AME_flat$variable, paste0(AME_flat$factor, '(', AME_flat$level, ')')))
    
    out <- rbind(out, AME_flat)
  }
  
  return(out)
}

#' Plot a FactorHet object
#' 
#' Plots the coefficients \eqn{\beta_k} from a fitted FactorHet object for the
#' main effects only. Use \code{\link{AME}} to calculate average marginal effects. This
#' provides a fast method to examine the impact of a factor.
#' 
#' @param object An object from \code{\link{FactorHet}} or \code{\link{FactorHet_mbo}}.
#' @param baseline A specification of baseline levels of each factor. The
#'   default is \code{NA} and shows all coefficients. The documentation for
#'   \code{\link{AME}} discusses how a named list could be provided to show the
#'   difference between coefficients.
#' @param plot A logical value as to whether the function should print the plot
#'   immediately or quietly provide an object containing the plot and data. The
#'   default is \code{TRUE}.
#' @seealso \link{AME}
#' @return \item{plot}{The ggplot2 visualization of the coefficients.}
#' \item{data}{The data used to create the ggplot2 visualization.}
#' 
#' @examples 
#' # Fit with one group and limited regularization for example only
#' # Ignore conjoint structure for simplicity
#' fit <- FactorHet(Chosen_Immigrant ~ Plans + Ed + Country,
#'   design = immigration, lambda = 1e-2,
#'   K = 1, group = ~ CaseID, task = ~ contest_no, choice_order = ~ choice_id)
#' # Plot the raw coefficients
#' cjoint_plot(fit)
#' 
#' @import ggplot2
#' @export
cjoint_plot <- function(object, baseline = NA, plot = TRUE) {
  return(do.call('internal_cjoint', as.list(environment())))
  # if (inherits(object, 'FactorHet_crossfit')){
  #   do.call('internal_cjoint_crossfit', as.list(environment()))
  # }else{
  #   do.call('internal_cjoint', as.list(environment()))
  # }
}

internal_cjoint <- function(object, baseline, plot){
  
  if (inherits(object, 'FactorHet')){
    
    object_factor_levels <- object$internal_parameters$factor_levels
    opts <- list()
    opts$verbose <- FALSE
    opts <- c(opts, list('baseline' = baseline))
    opts$object <- object
    object <- do.call('extract_coefficient', args = opts)
    rm(opts)
  }else{stop('Should pass FactorHet object to cjoint_plot')}
  
  joint_levels <- unlist(mapply(names(object_factor_levels), object_factor_levels, FUN=function(i,j){paste(i,j, sep ='_@_')}))
  object$joint_level <- paste(object$factor, object$level, sep='_@_')
  object$joint_level <- factor(object$joint_level, joint_levels)
  object$fmt_group <- paste0('Group ', object$group)
  parse_lbl <- function(x){sapply(strsplit(x, split='_@_'), FUN=function(i){i[2]})}
  
  object$ll <- object$coef - 1.96 * sqrt(object$var)
  object$ul <- object$coef + 1.96 * sqrt(object$var)
  g <- make_plot_cjoint(object = object, parse_lbl = parse_lbl)
  # labs(col = 'Group:', shape = 'Group:')
  if (plot){
    print(g)
  }
  out <- (list(plot = g, data = object[, c('factor', 'level', 'group', 'coef', 'var')]))
  class(out) <- 'FactorHet_vis'
  return(out)
}

make_plot_cjoint <- function(object, parse_lbl){
  
  .data <- NULL
  # New version without aes_string...
  g <- ggplot(object) + 
    geom_hline(aes(yintercept=0)) +
    geom_point(aes(x=.data[['joint_level']], 
                   y=.data[['coef']], 
                   col =.data[['fmt_group']], 
                   pch =.data[['fmt_group']])) +
    geom_errorbar(aes(x=.data[['joint_level']],
                             ymin=.data[['ll']],
                             ymax=.data[['ul']],
                             col = .data[['fmt_group']])) +
    coord_flip() + facet_grid(factor ~ fmt_group, 
        scales = 'free_y', space = 'free_y', switch = 'y', 
        labeller = label_wrap_gen()) +
    theme_bw(base_size = 8) +
    theme(strip.text.y.left = element_text(angle = 0), 
          panel.spacing = unit(0.1, 'lines'), strip.placement = 'outside') +
    scale_color_discrete(guide = 'none') +
    scale_x_discrete(labels = parse_lbl) +
    xlab('Factor') + ylab('Coefficient') +
    theme(legend.position = 'none') 
  return(g)
  # Old Verison using aes_string...
  # ggplot(object) + 
  #   geom_hline(aes(yintercept=0)) +
  #   geom_point(aes_string(x='joint_level', y='coef', col = 'fmt_group', pch = 'fmt_group')) +
  #   geom_errorbar(aes_string(x='joint_level',ymin='ll',ymax='ul',
  #                            col = 'fmt_group')) +
  #   coord_flip() + facet_grid(factor ~ fmt_group, scales = 'free_y', space = 'free_y', switch = 'y', 
  #                             labeller = label_wrap_gen()) +
  #   theme_bw(base_size = 8) +
  #   theme(strip.text.y.left = element_text(angle = 0), 
  #         panel.spacing = unit(0.1, 'lines'), strip.placement = 'outside') +
  #   scale_color_discrete(guide = 'none') +
  #   scale_x_discrete(labels = parse_lbl) +
  #   xlab('Factor') + ylab('Coefficient') +
  #   theme(legend.position = 'none') 
}

make_plot_AME <- function(store_mfx, ignore_bl = FALSE, bl_shape = 16, facet = TRUE){
  
  plot_bl <- store_mfx[store_mfx$baseline,]
  if (ignore_bl){
    plot_bl <- NULL
    bfun <- NULL
    alpha_type <- NULL
    plot_alpha_scale <- NULL
    plot_vertical <- geom_hline(aes(yintercept=0), linetype = 'dashed') 
  }else if (nrow(plot_bl) == 0){
    plot_bl <- NULL
    bfun <- NULL
    alpha_type <- NULL
    plot_alpha_scale <- NULL
    plot_vertical <- geom_hline(aes(yintercept=0.5), linetype = 'dashed') 
  }else{
    bfun <- function(x){x[!x$baseline,]}
    opp.bfun <- function(x){x[x$baseline,]}
    plot_bl <- geom_point(data = opp.bfun, pch = bl_shape)
    plot_vertical <- geom_hline(aes(yintercept=0), linetype = 'dashed') 
    alpha_type <- 'not_baseline'
    store_mfx[['not_baseline']] <- ! store_mfx[['baseline']]
    plot_alpha_scale <- scale_alpha_manual(values = c(0, 1), guide = 'none')
  }
  
  if (!is.null(alpha_type)){
    plot_point <- geom_point(aes(
      col = .data[['group']], 
      pch = .data[['group']], 
      alpha = .data[[alpha_type]]))
  }else{
    plot_point <- geom_point(aes(
      col = .data[['group']], 
      pch = .data[['group']]
    ))
  }
  # Using .data instead of aes_string
  g <- ggplot(store_mfx,
      aes(
        x=.data[['fmt_level']], 
        ymin=.data[['ll']],
        ymax=.data[['ul']],
        y=.data[['marginal_effect']])
    ) +
    plot_point + 
    plot_alpha_scale +
    plot_bl +
    geom_errorbar(data = bfun,
                  aes(col = .data[['group']])) +
    plot_vertical +
    coord_flip() +    
    theme_bw(base_size = 8) +
    theme(strip.text.y.left = element_text(angle = 0), 
          panel.spacing = unit(0.1, 'lines'), strip.placement = 'outside') +
    xlab('Factor') + ylab('Effect') +
    theme(legend.position = 'none')
  
  # Previous ggplot2 using aes_string...    
  # g <- ggplot(store_mfx,
  #             aes_string(x='fmt_level', ymin='ll',
  #                        ymax='ul',
  #                        y='marginal_effect')) +
  #   geom_point(aes_string(col = 'group', pch = 'group', 
  #                         alpha = alpha_type)) +
  #   plot_alpha_scale +
  #   plot_bl +
  #   geom_errorbar(data = bfun, aes_string(col = 'group')) +
  #   plot_vertical +
  #   coord_flip() +    
  #   theme_bw(base_size = 8) +
  #   theme(strip.text.y.left = element_text(angle = 0), 
  #         panel.spacing = unit(0.1, 'lines'), strip.placement = 'outside') +
  #   xlab('Factor') + ylab('Effect') +
  #   theme(legend.position = 'none')
  if (facet){
    g <- g + facet_grid(factor ~ paste0('Group ', group), 
     scales = 'free_y', space = 'free_y', switch = 'y',
     labeller = label_wrap_gen())
  }else{
    g <- g + facet_grid(factor ~ ., 
      scales = 'free_y', space = 'free_y', switch = 'y', 
      labeller = label_wrap_gen()) +
      theme(legend.position = 'bottom')
  }
  return(g)
}
#' Calculate marginal effects
#' 
#' Calculate the average marginal (component) effect (AME or AMCE), the average
#' combination effect (ACE), or the average marginal interaction effect (AMIE)
#' with a FactorHet model.
#' 
#' @rdname calc_effects
#' 
#' @param object An object from \code{\link{FactorHet}} or \code{\link{FactorHet_mbo}}.
#' @param design A dataset used to estimate the marginal effects. The default,
#'   \code{NULL}, uses the estimation data.
#' @param baseline A named list consisting of each factor and a corresponding
#'   baseline level. The default (\code{NULL}) computes the effect for all
#'   factors using the first level as the baseline. \code{NA} uses no baseline
#'   and approximates the "marginal means" from Leeper et al. (2020).
#' @param vcov A logical value indicating whether the standard errors for the
#'   AME should be computed. The default is \code{TRUE}. Standard errors are not
#'   implemented for the AMIE.
#' @param vcov.type A string indicating the type of standard errors to be
#'   computed. The default is \code{NULL} and uses the default settings in
#'   \code{\link{vcov.FactorHet}}; options are specified by that function's
#'   \code{se.method} argument.
#' @param average_position A logical value indicating whether, for forced choice
#'   designs, the "left" and "right" profiles should be averaged. The default is
#'   \code{TRUE}. Goplerud et al. (2025) provide discussion of this point.
#' @param plot A logical value as to whether the function should print the plot
#'   immediately or quietly provide an object containing the plot and data. The
#'   default is \code{TRUE}.
#' @param verbose A logical value as to whether more information should be
#'   provided on the progress of estimating the effects. The default is
#'   \code{TRUE}.
#' @param ignore_restrictions A logical value about whether to ignore
#'   randomization restrictions when calculating the marginal effects. The
#'   default is \code{FALSE}. "Details" provides more information.
#' @param extra_restriction A list of additional restrictions to include when
#'   computing the marginal effects. The default is \code{NULL}, i.e. no
#'   additional restrictions. "Details" provides more information about the use
#'   of this function.
#' 
#' @references 
#'   
#'   Egami, Naoki and Kosuke Imai. 2019. "Causal Interaction in
#'   Factorial Experiments: Application to Conjoint Analysis." \emph{Journal of the
#'   American Statistical Association}. 114(526):529-540.
#'
#'   Goplerud, Max, Kosuke Imai, and Nicole E. Pashley. 2025. "Estimating
#'   Heterogeneous Causal Effects of High-Dimensional Treatments: Application to
#'   Conjoint Analysis." arxiv preprint: \url{https://arxiv.org/abs/2201.01357}
#'  
#'   Leeper, Thomas J., Sara B. Hobolt, and James Tilley. 2020. "Measuring Subgroup
#'   Preferences in Conjoint Experiments." \emph{Political Analysis}. 28(2):207-221.
#'
#' @return \item{plot}{The ggplot2 visualization of the marginal effects.}
#' \item{data}{The data used to create the ggplot2 visualization.}
#' 
#' @details 
#' 
#' \bold{Choice of Baseline}: For ACE and AMIE, a choice of baseline is
#' required. See Egami and Imai (2019) for details. For AME, a choice of
#' baseline corresponds to a "standard" AME (see Egami and Imai 2019). The
#' option \code{NULL} chooses the first level of each factor. It can be manually
#' specified using a named list. If a named list is provided, only AMEs for
#' those named factors are calculated. This can be helpful if there are many
#' factors. 
#' 
#' If \code{NA} is provided as the baseline level, the AME is calculated without
#' a baseline; while this does not correspond to a "proper" AME, it is designed
#' to approximate the "marginal means" discussed in Leeper et al. (2020). Note
#' that in the presence of randomization restrictions, the quantity estimated
#' with a \code{NA} baseline may not be centered around 0.5. Ignoring the
#' randomization restrictions may be useful in this scenario. Supporting information
#' from Goplerud et al. (2025) provides more discussion of this point.
#' 
#' \bold{Randomization Restrictions}: Randomization restrictions can be set in
#' one of two ways. By default, FactorHet checks whether for each pairwise
#' combinations of factors, some combination of levels do not occur at all (e.g.
#' "doctor" and "high school") or whether some included interactions are
#' extremely rare (see \code{rare_threshold} in \code{\link{FactorHet_control}}). Those
#' are assumed to be the randomization restrictions implied by the design.
#' Setting \code{rare_threshold = 0} forces the inclusion of all interaction
#' terms.
#' 
#' If this procedure for automatically generating randomization restrictions is
#' inappropriate for a specific dataset, randomization restrictions can be set
#' manually as follows using the \code{manual_AME} function. First, set
#' \code{ignore_restrictions = TRUE}. This will ignore all "data-driven"
#' estimates of randomization restrictions. Second, the argument
#' \code{extra_restriction} should be a named list where the name of each
#' element corresponds to a factor (e.g. "Job") and each element is a vector of
#' the levels that \emph{cannot} be used. When using this approach, \code{AME}
#' should be used only for one factor at a time. An example is shown below.
#'
#' \bold{Plots}: Note that for the ggplot2 visualizations of the ACE and AMIE,
#' gray squares indicate combinations that are excluded due to randomization
#' restrictions. White indicates baseline levels.
#' @examples 
#' data(immigration)
#' # Induce "fake" randomization restriction
#' immigration$joint_id <- paste(immigration$CaseID, immigration$contest_no)
#' remove_profiles <- subset(immigration, Plans == 'No plans' & Ed == 'GradDeg')
#' immigration <- subset(immigration, !(joint_id %in% remove_profiles$joint_id))
#' # Fit with one group and limited regularization for example only
#' fit <- FactorHet(Chosen_Immigrant ~ Plans + Ed + Country,
#'   design = immigration, lambda = 1e-4,
#'   K = 1, group = ~ CaseID, task = ~ contest_no, choice_order = ~ choice_id)
#' # Estimate AME of "Ed" with randomization restriction
#' est_AME <- AME(fit, baseline = list('Ed' = 'GradDeg'))
#' \donttest{
#' # Estimate AME ignoring randomization restriction
#' est_AME_norr <- AME(fit, 
#'   baseline = list('Ed' = 'GradDeg'), ignore_restrictions = TRUE)
#' # Estimate AME by manually specifying randomization restrictions
#' # this uses the 'manual_AME' function
#' est_AME_rr_manual <- manual_AME(fit,
#'   baseline = list('Ed' = 'GradDeg'), ignore_restrictions = TRUE,
#'   extra_restriction = list('Plans' = 'No plans'))
#' stopifnot(isTRUE(all.equal(est_AME$data, est_AME_rr_manual$data)))
#' # Estimate without baseline
#' est_MM <- AME(fit, baseline = list('Ed' = NA))
#' }
#' # Estimate ACE and AMIE
#' \donttest{
#' est_ACE <- ACE(fit, baseline = list('Ed' = 'GradDeg', 'Plans' = 'Has contract'))
#' }
#' est_AMIE <- AMIE(fit, baseline = list('Ed' = 'GradDeg', 'Plans' = 'Has contract'))
#' @export
AME <- function(object, baseline = NULL, vcov = TRUE, 
                         design = NULL,
                         ignore_restrictions = FALSE, vcov.type = NULL,
                         average_position = TRUE, verbose = TRUE, plot = TRUE){
  # if (inherits(object, 'FactorHet_crossfit')){
  #   do.call('internal_AME_crossfit', as.list(environment()))
  # }else{
  #   do.call('internal_AME', as.list(environment()))
  # }
  args <- as.list(environment())
  return(do.call('internal_AME', args))
}

#' @rdname calc_effects
#' @export
manual_AME <- function(object, baseline, vcov = TRUE, 
                       design = NULL, extra_restriction = NULL,
                       ignore_restrictions = FALSE, vcov.type = NULL,
                       average_position = TRUE, verbose = TRUE, plot = TRUE){
  args <- as.list(environment())
  if (length(args$baseline) > 1){
    warning('manual_AME should generally only be used to compute one AME at a time.')
  }
  return(do.call('internal_AME', args))
}

#' Deprecated Functions
#'
#' The names of these functions have changed, e.g., \code{AME} to
#' \code{AME}, for clarity. Please adjust your code to use the current names.
#' 
#' @rdname deprecated
#' @keywords internal
#' @export
marginal_AME <- function(...){
  .Deprecated(new = 'AME', old = 'marginal_AME')
  return(AME(...))
}

internal_AME <- function(object, baseline = NULL, vcov = TRUE, 
                         design = NULL,
                         ignore_restrictions = FALSE, vcov.type = NULL,
                         average_position = TRUE, verbose = TRUE, plot = TRUE,
                         extra_restriction = NULL){
  
  if (!inherits(object, 'FactorHet')){
    stop('Should pass FactorHet object to AME')
  }
  
  object_factor_levels <- object$internal_parameters$factor_levels
  J <- length(object_factor_levels)
  K <- ncol(coef.FactorHet(object))
  
  rare_cols <- object$internal_parameters$rare$rare_fmt_col

  restriction_list <- build_restrictions(removed_cols = rare_cols,
                     factor_names = names(object_factor_levels))
  
  if (ignore_restrictions){
    warning('Overriding restrictions implied by restricted randomization')
    restriction_list <- lapply(names(object_factor_levels), FUN=function(i){NULL})
    names(restriction_list) <- names(object_factor_levels)
  }
  
  if (is.null(baseline)){
    baseline <- lapply(object_factor_levels, FUN=function(i){i[1]})
    if (verbose){
      message('baseline = NULL gives the following baseline')
      print(unlist(baseline))
    }
  } else if (identical(baseline, NA)) {
    baseline <- lapply(object_factor_levels, FUN=function(i){NA})
    if (verbose){
      message('baseline = NA approximates "marginal means".')
    }
  } else{
    
    if (!isTRUE(all(names(baseline) %in% names(object_factor_levels)))){
      stop("baseline must be a named list with one entry for each factor.")
    }
    if (!isTRUE(all(lengths(baseline) == 1))){
      stop('baseline must be a named list with one entry for each factor.')
    }
    
    verify_baseline <- all(mapply(baseline, names(baseline), FUN=function(b, l){b %in% c(NA, object_factor_levels[[l]])}))
    if (!verify_baseline){stop('All baseline levels must be in the associated factor or "NA".')}
    
  }  
  
  # For custom training data:
  if (is.null(design)){
    baseline_design <- data.frame(object$internal_parameters$data$design, check.names = FALSE)
  }else{
    select_cols <- c(
      names(object_factor_levels), unlist(formula(object)$other_parameters)
    )
    diff_cols <- setdiff(select_cols, names(design))
    if (length(diff_cols) != 0){
      stop(paste('Missing columns from predict data: ', paste(diff_cols, collapse=', ')))
    }
    baseline_design <- na.omit(design[, select_cols])
  }
  
  unique_choice <- object$internal_parameters$unique_choice
  if (is.null(unique_choice)){
    pure_factorial <- TRUE
    unique_choice <- 'factorial'
  }else{
    pure_factorial <- FALSE
  }
  
  choice_name <- object$internal_parameters$formula$other_parameters$name_choice_order
  raw_design <- baseline_design
  
  if (vcov){
    if (is.null(vcov.type)){
      full_object_vcov <- vcov.FactorHet(object, phi = FALSE)
    }else{
      full_object_vcov <- vcov.FactorHet(object, phi = FALSE, se.method = vcov.type)
    }
    if (is.null(full_object_vcov)){
      vcov <- FALSE
    }else{
      n_p <- nrow(coef(object))
      stopifnot(ncol(full_object_vcov) == (n_p *K))
      
      object_vcov <- lapply(1:K, FUN=function(k){
        ind_k <- 1:n_p + (k-1) * n_p
        return(full_object_vcov[ind_k, ind_k])
      })
      
      full_object_vcov <- full_object_vcov[1:(n_p * K), 1:(n_p * K)]
      
    }
  }else{
    object_vcov <- NULL
    full_object_vcov <- NULL
  }

  
  store_mfx <- data.frame()
  vcov_counter <- 1
  store_grad <- list()
  for (fac in names(baseline)){
    restricted_fac <- restriction_list[[fac]]
    # if (length(restricted_fac) > 1){stop()}
    
    if (!is.null(extra_restriction)){
      if (is.null(restricted_fac)){restricted_fac <- list()}
      restricted_fac <- c(restricted_fac, extra_restriction)
      un <- unique(names(restricted_fac))
      restricted_fac <- lapply(un, FUN=function(i){unique(unlist(restricted_fac[names(restricted_fac) == i]))})
      names(restricted_fac) <- un
      restricted_fac <- restricted_fac[unique(names(restricted_fac))]
    }
    no_restrictions <- is.null(restricted_fac) | (length(restricted_fac) == 0)
    
    #Set each profile to "l" - what is the probability of choosing "baseline"
    #if no effect of factor j
    baseline_fac <- baseline[[fac]] #Get the baseline for j

    baseline_data <- prune_empirical_dist(baseline_design = raw_design, 
       cjoint_names = formula(object)$other_parameters,
       unique_choice = unique_choice, no_restrictions = no_restrictions,
       pure_factorial = pure_factorial, restricted_fac = restricted_fac,
       verbose = verbose
      )

    if (any(sapply(baseline_data, is.na))){
      stop(paste0('Randomization restrictions removed all observations for "', fac, '". Either adjust "rare_threshold" and refit model or allow AME to ignore restrictions.'))
    }
    raw_treat_design <- baseline_data

    if (!is.na(baseline_fac)){
      
      if (pure_factorial){
        
        baseline_data[['factorial']]$data[,fac] <- baseline_fac
        
        predict_baseline <- predict(object, newdata = baseline_data[['factorial']]$data, 
                                    return_task = TRUE, by_group = TRUE,
                                    ignore_moderator = TRUE,
                                    calc_gradient = vcov)
        if (vcov){
          grad_baseline <- delta_grad_AME(predict_baseline)
        }
        
        predict_baseline <- colMeans(predict_baseline)
        
        baseline_data$factorial <- c(baseline_data$factorial, 
                                     list('gradient' = grad_baseline, 'predict' = predict_baseline))
        rm(predict_baseline, grad_baseline)
        
      }else{
        
        baseline_data <- lapply(unique_choice, FUN=function(p){
          
          bd_p <- baseline_data[[p]]$data
          
          # Set the p (e.g. left) to the baseline level
          bd_p[bd_p[,choice_name] == p, fac] <- baseline_fac
          
          predict_bd_p <- predict(object, newdata = bd_p, ignore_moderator = TRUE,
                                  return_task = TRUE, by_group = TRUE,
                                  calc_gradient = vcov)
          if (vcov){
            grad_bd_p <- delta_grad_AME(predict_bd_p)
          }else{
            grad_bd_p <- NULL
          }
          predict_bd_p <- colMeans(predict_bd_p)
          
          return(list('data' = bd_p,
                      'predict' = predict_bd_p,
                      'gradient' = grad_bd_p))
        })
        names(baseline_data) <- unique_choice      
      }
      
    }else{
      baseline_data <- NA
    }
    
    
    #Loop over levels in factor j that are not baseline:
    if (verbose){message('.', appendLF = FALSE)}
    
    for (l in setdiff(object_factor_levels[[fac]], baseline_fac)){
      #Loop over "left and right" configurations
      for (p in unique_choice){
        treat_design <- raw_treat_design[[p]]$data
        if (pure_factorial){
          treat_design[,fac] <- l
        }else{
          #Set profile "p" (e.g. left) to level j
          treat_design[treat_design[,choice_name] == p, fac] <- l
        }
        predict_treat <- predict(object, newdata = treat_design, 
                                 by_group = TRUE, ignore_moderator = TRUE,
                                 calc_gradient = vcov)
        if (vcov){
          grad_treat <- delta_grad_AME(predict_treat)
          if (!is.na(baseline_fac)){
            grad_diff <- grad_treat - baseline_data[[p]]$gradient
            if (!pure_factorial){
              grad_diff <- ifelse(p == unique_choice[2], 1, -1) * grad_diff
            }
          }else{
            grad_diff <- grad_treat
            if (!pure_factorial){
              grad_diff <- ifelse(p == unique_choice[2], 1, -1) * grad_diff
            }
            
          }
          store_grad[vcov_counter + 0:(K-1)] <- list_from_cols(grad_diff)
          vcov_counter <- vcov_counter + K
        }
        
        predict_treat <- colMeans(predict_treat)
        if (!is.na(baseline_fac)){
          diff_predict <- predict_treat - baseline_data[[p]]$predict
          if (!pure_factorial){
            diff_predict <- ifelse(p == unique_choice[2], 1, -1) * diff_predict
          }
        }else{
          diff_predict <- predict_treat
          if (!pure_factorial){
            if (p == unique_choice[2]){
              diff_predict <- diff_predict
            }else{
              diff_predict <- 1 - diff_predict
            }
          }
        }
        
        store_mfx <- rbind(store_mfx, 
          data.frame(marginal_effect = diff_predict, 
           factor = fac, level = l, stringsAsFactors = FALSE, switch = p, group = 1:K))
      }
    }
    
  }
  if (!pure_factorial){
    store_mfx$group <- as.numeric(store_mfx$group)
  }
  if (average_position){
    joint_id <- with(store_mfx, paste(factor, level, group, sep = '@@@'))

    if (vcov){
      store_grad <- lapply(split(1:nrow(store_mfx), joint_id), FUN=function(i){
        grad_i <- Reduce('+', store_grad[i]) * 1/2
      })
    }
    
    store_mfx <- data.frame(marginal_effect = sapply(split(store_mfx$marginal_effect, joint_id), mean))
    store_mfx[,c('factor', 'level', 'group')] <- do.call('rbind', strsplit(rownames(store_mfx), split = '@@@'))
    store_mfx$group <- as.numeric(store_mfx$group)
    
    if (vcov){#Get the standard errors
      stopifnot(identical(rownames(store_mfx), names(store_grad)))
      store_mfx$var <- mapply(store_grad, store_mfx$group, FUN=function(i,k){
        as.numeric(t(i) %*% object_vcov[[k]] %*% i)
      })
    }else{
      store_mfx$var <- NA
    }
    rownames(store_mfx) <- NULL
  }else{
    if (vcov){#Get the standard errors
      store_mfx$var <- mapply(store_grad, store_mfx$group, FUN=function(i,k){
        as.numeric(t(i) %*% object_vcov[[k]] %*% i)
      })
    }else{
      store_mfx$var <- NA
    }
  }
  store_mfx$baseline <- FALSE
  copy_mfx <- store_mfx
  # Replace very small negative numbers with "0"
  numerical_zero <- (sign(store_mfx$var) == -1) & (sqrt(abs(store_mfx$var)) < sqrt(.Machine$double.eps))
  store_mfx$var[numerical_zero] <- 0

  #Add in baseline
  baseline_mfx <- data.frame(marginal_effect = 0, baseline = TRUE, level = unlist(baseline), factor = names(baseline), stringsAsFactors = F, row.names = NULL)
  baseline_mfx <- do.call('rbind', lapply(1:K, FUN=function(i){baseline_mfx$group <- i; return(baseline_mfx)}))
  baseline_mfx$var <- NA
  baseline_mfx <- baseline_mfx[!is.na(baseline_mfx$level),]
  store_mfx <- rbind(store_mfx, baseline_mfx)
  
  fmt_order <- c(unlist(baseline), setdiff(unlist(object_factor_levels), unlist(baseline)))
  store_mfx$fmt_level <- factor(store_mfx$level, levels = unique(fmt_order), ordered = TRUE)
  store_mfx$ll <- store_mfx$marginal_effect - 1.96 * sqrt(store_mfx$var)
  store_mfx$ul <- store_mfx$marginal_effect + 1.96 * sqrt(store_mfx$var)
  store_mfx$group <- factor(store_mfx$group)
  
  plot_bl <- store_mfx[store_mfx$baseline,]
  if (nrow(plot_bl) == 0){
    plot_bl <- NULL
    plot_vertical <- geom_hline(aes(yintercept=0.5), linetype = 'dashed') 
  }else{
    plot_bl <- geom_point(data = store_mfx[store_mfx$baseline,])
    plot_vertical <- geom_hline(aes(yintercept=0), linetype = 'dashed') 
  }
  g <- make_plot_AME(store_mfx)
  g_nofacet <- make_plot_AME(store_mfx, facet = FALSE)
  if (plot){print(g)}
  
  out <- list(plot = g, plot_nofacet = g_nofacet,
    data = store_mfx,
    gradient = store_grad)
  class(out) <- 'FactorHet_vis'
  invisible(out)
}

#' Difference between AMEs in each group
#' 
#' @description Computes the differences between AME between two groups with an
#'   accompanying standard error.
#'   
#' @rdname diff_AME
#' @param object An object from \code{\link{FactorHet}} or
#'   \code{\link{FactorHet_mbo}}.
#' @param AME An object containing the average marginal effects estimated using
#'   \code{\link{AME}}.
#' @param baseline_group An integer that denotes the baseline group. The
#'   function will show the difference in AMEs from this group, i.e. Group
#'   \code{k} - Group \code{baseline_group}.
#' @param plot A logical value as to whether the function should print the plot
#'   immediately or quietly provide an object containing the plot and data. The
#'   default is \code{TRUE}.
#' 
#' @examples
#'   # Fit with two groups and fixed lambda for quick illustration
#'   \donttest{
#'   data(immigration)
#'   set.seed(15)
#'   fit <- FactorHet(formula = Chosen_Immigrant ~ Country + Ed,
#'     design = immigration, group = ~ CaseID, task = ~ contest_no,
#'     choice_order = ~ choice_id, lambda = 1e-2,
#'     control = FactorHet_control(init_method = 'random_member'),
#'     K = 2)
#'   fit_AME <- AME(fit)
#'   diff_AME(fit, fit_AME, baseline_group = 2)
#'   }
#' @export
diff_AME <- function(object, AME, baseline_group,
                     plot = TRUE){
  
  if (!inherits(object, 'FactorHet')){
    stop('Must be object estimated using FactorHet or FactorHet_mbo')
  }
  if (!inherits(AME, 'FactorHet_vis')){
    stop('Must be objected estimated using AME')
  }
  # Get unique combinations
  estimates <- AME$data
  gradients <- AME$gradient
  vcov <- vcov(object, phi = FALSE)

  unique_combinations <- unique(estimates[estimates$baseline == FALSE, c('factor', 'level')])
  fmt_order <- levels(estimates$fmt_level)
  K <- length(unique(estimates$group))
  seq_K <- 1:K
  nonbaseline_K <- setdiff(seq_K, baseline_group)
  
  out_diff <- lapply(1:nrow(unique_combinations), FUN=function(v){
    combo_v <- unique_combinations[v,]
    pos_v <- which((estimates$factor == combo_v$factor) & (estimates$level == combo_v$level))
    mfx_v <- estimates[pos_v, ]  
    est_v <- mfx_v$marginal_effect
    var_v <- mfx_v$var
    group_v <- as.numeric(mfx_v$group)
    grad_v <- gradients[pos_v]
    grad_v <- do.call('cbind', grad_v)
    stopifnot(isTRUE(all.equal(group_v, seq_K)))
    
    out_v <- do.call('rbind', lapply(nonbaseline_K, FUN=function(k){
      
      sign_k <- rep(0, K)
      sign_k[k] <- 1
      sign_k[baseline_group] <- -1
      transf_grad_v_k <- as.vector(grad_v %*% Diagonal(x = sign_k))
      
      diff_est_v <- sum(est_v * sign_k)
      diff_var_v <- as.numeric(t(transf_grad_v_k) %*% vcov %*% transf_grad_v_k)
      return(data.frame(est = diff_est_v,
                        var = diff_var_v, 
                        group = k, 
                        stringsAsFactors = FALSE))  
    }))
    out_v$level <- combo_v$level
    out_v$factor <- combo_v$factor
    # out_v$group <- paste0(out_v$group, '-', baseline_group)
    return(out_v)
  })
  out_diff <- do.call('rbind', out_diff)
  out_diff$fmt_level <- factor(out_diff$level, levels = unique(fmt_order), ordered = TRUE)
  out_diff$ll <- out_diff$est - 1.96 * sqrt(out_diff$var)
  out_diff$ul <- out_diff$est + 1.96 * sqrt(out_diff$var)
  out_diff$group <- factor(out_diff$group)
  out_diff$marginal_effect <- out_diff$est
  out_diff$baseline <- FALSE
  g_diff <- make_plot_AME(out_diff, ignore_bl = TRUE)
  g_diff <- g_diff + ylab(paste0('Difference in Effect from Group ', baseline_group))
  if (plot){print(g_diff)}
  out <- list(plot = g_diff, 
              data = out_diff)
  class(out) <- 'FactorHet_vis'
  invisible(out)
}

#' @rdname deprecated
#' @keywords internal
#' @export
marginal_ACE <- function(...){
  .Deprecated(new = 'ACE', old = 'marginal_ACE')
  return(ACE(...))
}

#' @rdname calc_effects
#' @export
ACE <- function(object, baseline, design = NULL, average_position = TRUE, 
                         ignore_restrictions = FALSE, extra_restriction = NULL,
                         verbose = TRUE, plot = TRUE){
  
  object_factor_levels <- object$internal_parameters$factor_levels
  J <- length(object_factor_levels)
  K <- ncol(coef(object))
  
  rare_cols <- object$internal_parameters$rare$rare_fmt_col
  restriction_list <- build_restrictions(removed_cols = rare_cols,
                                         factor_names = names(object_factor_levels))
  
  if (ignore_restrictions){
    warning('Overriding restrictions implied by restricted randomization')
    restriction_list <- lapply(names(object_factor_levels), FUN=function(i){NULL})
    names(restriction_list) <- names(object_factor_levels)
  }
  
  if (length(baseline) != 2){stop('Must provide exactly 2 factors as baseline for ACE.')}

  if (!isTRUE(all(names(baseline) %in% names(object_factor_levels)))){
    stop('baseline must be a named list with one entry for each of the two factor_levels. Some of the factors in "baseline" are not in the original model.')
  }

  verify_baseline <- all(mapply(baseline, names(baseline), FUN=function(b, l){b %in% object_factor_levels[[l]]}))
  if (!verify_baseline){stop('All baseline levels must be in the associated factor.')}

  
  unique_choice <- object$internal_parameters$unique_choice
  if (is.null(unique_choice)){
    pure_factorial <- TRUE
    unique_choice <- 'factorial'
  }else{
    pure_factorial <- FALSE
  }
  
  choice_name <- object$internal_parameters$formula$other_parameters$name_choice_order
  
  fac_l <- names(baseline)[1]
  fac_m <- names(baseline)[2]
  
  baseline_fac_l <- baseline[[fac_l]] #Get the baseline for j
  baseline_fac_m <- baseline[[fac_m]] #Get the baseline for j
  
  # For custom training data:
  if (is.null(design)){
    baseline_design <- data.frame(object$internal_parameters$data$design, check.names = FALSE)
  }else{
    baseline_design <- na.omit(design[, names(object_factor_levels)])
  }

  restricted_fac <- list()

  for (f in c(fac_l, fac_m)){
    restricted_fac <- c(restricted_fac, restriction_list[[f]])
  }
  if (!is.null(extra_restriction)){
    restricted_fac <- c(restricted_fac, extra_restriction)
  }
  un <- unique(names(restricted_fac))
  restricted_fac <- lapply(un, FUN=function(i){unique(unlist(restricted_fac[names(restricted_fac) == i]))})
  names(restricted_fac) <- un
  # Get all factors and levels that must be excluded from ACE
  restricted_fac <- restricted_fac[unique(names(restricted_fac))]
  if (length(restricted_fac) == 0){
    no_restrictions <- TRUE
    restricted_fac <- NULL
  }else{
    no_restrictions <- FALSE
  }
  
  baseline_data <- prune_empirical_dist(baseline_design = baseline_design, 
      cjoint_names = formula(object)$other_parameters,
      unique_choice = unique_choice, no_restrictions = no_restrictions,
      pure_factorial = pure_factorial, restricted_fac = restricted_fac,
      verbose = verbose
  )
  
  if (any(sapply(baseline_data, is.na))){
    stop(paste0('Randomization restrictions removed all observations for ', paste0(fac_l, ',', fac_m), '. Either adjust "rare_threshold" and refit model or allow AME to ignore restrictions.'))
  }
  
  if (pure_factorial){
    baseline_data$factorial$data[, fac_l] <- baseline_fac_l
    baseline_data$factorial$data[, fac_m] <- baseline_fac_m
    predict_baseline <- predict(object, 
      newdata = baseline_data$factorial$data, by_group = TRUE)
    baseline_data$factorial$predict <- colMeans(predict_baseline)
  }else{
    baseline_data <- lapply(unique_choice, FUN=function(p){
      bd_p <- baseline_data[[p]]$data
      # Set the p (e.g. left) to the baseline level
      bd_p[bd_p[,choice_name] == p, fac_l] <- baseline_fac_l
      bd_p[bd_p[,choice_name] == p, fac_m] <- baseline_fac_m
      
      predict_bd_p <- predict(object, newdata = bd_p, by_group = TRUE)
      predict_bd_p <- colMeans(predict_bd_p)
      return(list('data' = bd_p, 'predict' = predict_bd_p))
    })
    names(baseline_data) <- unique_choice
  }
  
  store_mfx <- data.frame()
  store_res <- data.frame()
  # Loop over levels

  for (l in object_factor_levels[[fac_l]]){
    # Skip restricted levels for j
    if (verbose){message('.', appendLF = FALSE)}
    for (m in object_factor_levels[[fac_m]]){
      # Skip baseline
      if ((l == baseline_fac_l) & (m == baseline_fac_m)){
        next
      }
      # Skip restricted levels for (l,m)
      candidate_restriction <- c(paste0(fac_l, '(', l, ')-', fac_m, '(', m, ')'), 
        paste0(fac_m, '(', m, ')-', fac_l, '(', l, ')'))
      if (any(candidate_restriction %in% rare_cols)){
        store_res <- rbind(store_res, 
          data.frame(factor_l = l, factor_m = m, stringsAsFactors = F))
        next
      }

      for (p in unique_choice){
        
        treat_design <- baseline_data[[p]]$data
        
        if (pure_factorial){
          treat_design[,fac_l] <- l
          treat_design[,fac_m] <- m
        }else{
          #Set profile "p" (e.g. left) to level j
          treat_design[treat_design[,choice_name] == p, fac_l] <- l
          treat_design[treat_design[,choice_name] == p, fac_m] <- m
        }
        #Get the probability of choosing left XX
        predict_treat <- predict(object, newdata = treat_design, by_group = TRUE)
        
        predict_treat <- colMeans(predict_treat)
        
        store_mfx <- rbind(store_mfx, data.frame(
          marginal_effect = predict_treat - baseline_data[[p]]$predict, 
          factor_l = l, factor_m = m, stringsAsFactors = FALSE, switch = p, group = 1:K))
      }
    }
  }
  if (!pure_factorial){
    store_mfx$marginal_effect <- with(store_mfx, 
      ifelse(switch == unique_choice[2], marginal_effect, -marginal_effect))
  }
  if (average_position){
    store_mfx <- data.frame(marginal_effect = sapply(with(store_mfx,
      split(marginal_effect, paste(factor_l, factor_m, group, sep = '@@@'))), mean))
    store_mfx[,c('factor_l', 'factor_m', 'group')] <- do.call('rbind', strsplit(rownames(store_mfx), split = '@@@'))
    store_mfx$group <- as.numeric(store_mfx$group)
    rownames(store_mfx) <- NULL
  }
  store_mfx$baseline <- FALSE
  #Add in baseline
  
  bl <- do.call('rbind', mapply(baseline, names(baseline), SIMPLIFY = FALSE, FUN=function(i,j){
    v <- as.data.frame(baseline, stringsAsFactors = F)
    v[,] <- NA
    v[[j]] <- i
    return(v)
  }))
  bl$marginal_effect <- NA
  bl$baseline <- TRUE
  bl$group <- NA
  rownames(bl) <- NULL
  
  store_mfx$factor_l <- factor(store_mfx$factor_l, levels = object_factor_levels[[fac_l]])
  store_mfx$factor_m <- factor(store_mfx$factor_m, levels = object_factor_levels[[fac_m]])
  
  fmt_order <- c(unlist(baseline), setdiff(unlist(object_factor_levels), unlist(baseline)))
  
  if (nrow(store_res) == 0){
    store_res <- NULL
  }
  
  g <- ggplot(subset(store_mfx, !baseline)) +
    geom_tile(aes(x=.data[['factor_l']],
                  y=.data[['factor_m']], 
                  fill = .data[['marginal_effect']])) +
    geom_tile(aes(x=.data[['factor_l']], 
                  y=.data[['factor_m']]), 
              fill = 'gray', alpha = 0.5, data = store_res) +
    facet_wrap(~paste0('Group ', group)) + theme_bw() +
    theme(legend.position = 'bottom', panel.grid = element_blank(),
          axis.text.x = element_text(hjust=1,vjust=0, angle = 90)) +
    labs(fill = 'ACE') +
    xlab(fac_l) + ylab(fac_m) +scale_y_discrete(drop = F) +
    scale_x_discrete(drop = F)
  
  if (plot){print(g)}
  names(store_mfx)[names(store_mfx) == 'factor_l'] <- fac_l
  names(store_mfx)[names(store_mfx) == 'factor_m'] <- fac_m
  
  store_mfx <- rbind(store_mfx, bl)

  out <- (list(plot = g, data = store_mfx, res = store_res))
  class(out) <- 'FactorHet_vis'
  return(out)
}


#' @rdname deprecated
#' @keywords internal
#' @export
marginal_AMIE <- function(...){
  .Deprecated(new = 'AMIE', old = 'marginal_AMIE')
  return(AMIE(...))
}

#' @rdname calc_effects
#' @export
AMIE <- function(object, design = NULL, baseline = NULL, average_position = TRUE, 
              ignore_restrictions = FALSE, verbose = FALSE, plot = TRUE){
  
  object_factor_levels <- object$internal_parameters$factor_levels
  J <- length(object_factor_levels)
  K <- ncol(coef(object))

  
  rare_cols <- object$internal_parameters$rare$rare_fmt_col
  restriction_list <- build_restrictions(removed_cols = rare_cols,
                                         factor_names = names(object_factor_levels))

  if (ignore_restrictions){
    warning('Overriding restrictions implied by restricted randomization')
    restriction_list <- lapply(names(object_factor_levels), FUN=function(i){NULL})
    names(restriction_list) <- names(object_factor_levels)
  }
  
  if (is.null(baseline)){
    baseline <- lapply(object_factor_levels, FUN=function(i){i[1]})
    if (verbose){
      message('baseline = NULL gives the following baseline')
      print(unlist(baseline))
    }
  }else{
    if (length(baseline) < 2){
      stop('Must provide at least two levels for AMIE; use "NULL" for all pairs.')
    }
    if (!isTRUE(all(names(baseline) %in% names(object_factor_levels)))){
      stop("baseline must be a named list with one entry for each factor.")
    }
    if (!isTRUE(all(lengths(baseline) == 1))){
      stop('baseline must be a named list with one entry for each factor.')
    }
    
    verify_baseline <- all(mapply(baseline, names(baseline), FUN=function(b, l){b %in% object_factor_levels[[l]]}))
    if (!verify_baseline){stop('All baseline levels must be in the associated factor.')}
    
  }  
  
  verify_baseline <- all(mapply(baseline, names(baseline), FUN=function(b, l){b %in% object_factor_levels[[l]]}))
  if (!verify_baseline){stop('All baseline levels must be in the associated factor.')}

  baseline_combinations <- combn(names(baseline), 2)
  
  all_AMIE <- lapply(1:ncol(baseline_combinations), FUN=function(i){
    # Get two factors
    combo_i <- baseline_combinations[,i]
    if (verbose){print(combo_i)}
    extra_res <- c()
    for (v in restriction_list[combo_i]){
      extra_res <- c(extra_res, v)
    }
    extra_res <- extra_res[unique(names(extra_res))]
    est_AME <- suppressWarnings(manual_AME(object = object, design = design,
                            extra_restriction = extra_res,
                            baseline = baseline[combo_i], plot = FALSE)$data)
    est_ACE <- ACE(object = object, design = design,
                   extra_restriction = extra_res,
                   baseline = baseline[combo_i], plot = FALSE)
    est_ACE_res <- est_ACE$res
    est_ACE <- est_ACE$data
    est_ACE$combination_effect <- est_ACE$marginal_effect
    est_ACE <- est_ACE[,!(names(est_ACE) %in% c('marginal_effect', 'level', 'factor'))]
    #Merge in the AMEs to subtract
    merge_l <- merge_m <- c('group', 'baseline', 'level')
    names(merge_m) <- names(merge_l) <- merge_l
    names(merge_l)[3] <- combo_i[1]
    names(merge_m)[3] <- combo_i[2]
    
    est_ACE <- merge(est_ACE, est_AME[est_AME$factor == combo_i[1],], 
                     by.x = names(merge_l), by.y = merge_l)
    est_ACE <- est_ACE[, !(names(est_ACE) %in% 'factor')]
    est_ACE <- merge(est_ACE, est_AME[est_AME$factor == combo_i[2],], 
                    by.x = names(merge_m), by.y = merge_m)
    est_ACE <- est_ACE[, !(names(est_ACE) %in% 'factor')]
    est_ACE <- est_ACE[!is.na(est_ACE$group),]
    est_ACE$marginal_effect.x <- with(est_ACE, ifelse(is.na(marginal_effect.x), 0, marginal_effect.x))
    est_ACE$marginal_effect.y <- with(est_ACE, ifelse(is.na(marginal_effect.y), 0, marginal_effect.y))
    
    est_ACE$AMIE <- with(est_ACE, combination_effect - marginal_effect.x - marginal_effect.y)
    est_ACE[[combo_i[1]]] <- factor(est_ACE[[combo_i[1]]], levels = object_factor_levels[[combo_i[1]]])
    est_ACE[[combo_i[2]]] <- factor(est_ACE[[combo_i[2]]], levels = object_factor_levels[[combo_i[2]]])
    
    if (!is.null(est_ACE_res)){
      extra_res_geom <- geom_tile(aes(x=.data[['factor_l']], 
                                      y=.data[['factor_m']]), 
                fill = 'gray', alpha = 0.5, data = est_ACE_res)
    }else{
      extra_res_geom <- NULL
    }

    g <- ggplot(est_ACE[!est_ACE$baseline,]) +
      geom_tile(aes(x=.data[[combo_i[1]]],
                    y=.data[[combo_i[2]]], 
                           fill = .data[['AMIE']])) +
      extra_res_geom +
      facet_wrap(~paste0('Group ', group)) + theme_bw() +
      theme(legend.position = 'bottom', panel.grid = element_blank(),
            axis.text.x = element_text(hjust=1,vjust=0, angle = 90)) +
      labs(fill = 'AMIE') +
      guides(fill = guide_colourbar(label.theme = element_text(size = 8, 
                                                               hjust = 1, angle = 90))) +
      scale_x_discrete(drop = F) + scale_y_discrete(drop = F)
    
    return(list(plot = g, data = est_ACE, factors = combo_i))
  })
  
  names(all_AMIE) <- apply(baseline_combinations, MARGIN = 2, FUN=function(i){paste(i, collapse =' ')})

  if (plot){
    all_plot <- lapply(all_AMIE, FUN=function(i){i$plot})
    print(all_plot)
  }else{
    all_plot <- NULL
  }
  out_AMIE <- list('data' = lapply(all_AMIE, FUN=function(i){i$data}),
       'plot' = all_plot)
  class(out_AMIE) <- 'FactorHet_vis'
  return(out_AMIE)
}

build_restrictions <- function(removed_cols, factor_names){

  if (is.null(removed_cols)){
    restricted_levels <- lapply(factor_names, FUN=function(i){NULL})
    names(restricted_levels) <- factor_names
    return(restricted_levels)
  }
  restricted_levels <- strsplit(removed_cols, split='(?<=\\))-', perl = T)
  restricted_levels <- lapply(restricted_levels, FUN=function(i){do.call('rbind', strsplit(i, split='\\(|\\)$'))})

  restricted_levels <- lapply(factor_names, FUN=function(j){
    #Get levels that have *some* restrictions
    r_j <- lapply(restricted_levels, FUN=function(i){
      match_j <- i[,1] == j
      has_j <- any(match_j)
      if (has_j){
        return(i[!match_j,])
      }else{
        return(NULL)
      }
    })
    r_j <- do.call('rbind', r_j)
    if (is.null(r_j)){return(r_j)}
    r_j <- lapply(split(r_j[,2], r_j[,1]), unique)
    # if (length(r_j) != 1){
    #   m <- paste0('Restrictions on multiple factors may behave strangely; this affects factor "', j, '".')
    #   warning(m)
    # }
    return(r_j)
  })
  names(restricted_levels) <- factor_names
  return(restricted_levels)
}

delta_grad_AME <- function(x){
  data_x <- attr(x, 'X')
  grad_x <- apply(x, MARGIN = 2, FUN=function(i){
    colMeans(Diagonal(x = i * (1 - i)) %*% data_x)
  })
  return(grad_x)
}

prune_empirical_dist <- function(baseline_design,
                                 cjoint_names, unique_choice, no_restrictions,
                                 pure_factorial, restricted_fac, verbose){
  
  if (no_restrictions){
    if (pure_factorial){
      return(list('factorial' = list('data' = baseline_design)))
    }else{
      out <- lapply(unique_choice, FUN=function(i){list('data' = baseline_design)})
      names(out) <- unique_choice
      return(out)
    }
  }
  if (pure_factorial){
    
    in_restrictions <- sapply(names(restricted_fac), FUN=function(v){
      baseline_design[[v]] %in% restricted_fac[[v]]
    })
    
    in_restrictions <- (rowSums(in_restrictions) > 0)
    
    # Return NA if 100% of observations are removed...
    if (all(in_restrictions)){
      return(NA)
    }
    
    report_restriction <- round(mean(in_restrictions) * 100, 2)
    
    if (verbose){
      message(paste0(report_restriction, '% of observations ',
                     'removed from the empirical distribution of AMCE\n', 'because of randomization restrictions.'))
    }
    baseline_data <- baseline_design[!in_restrictions,]
    baseline_data <- list('factorial' = list('data' = baseline_data))
    
  }else{
    
    joint_id <- paste(baseline_design[[cjoint_names$name_group]], 
                      baseline_design[[cjoint_names$name_task]], sep = '~!~')
    
    split_restricted <- lapply(names(restricted_fac), FUN=function(v){
      split(baseline_design[[v]], joint_id)
    })
    names(split_restricted) <- names(restricted_fac)
    stopifnot(all(sapply(split_restricted, FUN=function(i){all(lengths(i) == 2)})))
    split_choice <- split(baseline_design[[cjoint_names$name_choice_order]], joint_id)
    stopifnot(all(lengths(split_choice) == 2))
    # Loop over each choice
    # and remove profiles from THAT SIDE that contain and invalid profile.
    # Randomize over all observed ~SIDE profiles.
    baseline_data <- lapply(unique_choice, FUN=function(p){
      in_restriction <- sapply(names(split_restricted), FUN=function(sr){
        unsplit(mapply(split_choice, split_restricted[[sr]], FUN=function(choice_i, res_i){
          res_i[which(choice_i == p)] %in% restricted_fac[[sr]]
        }), joint_id)
      })
      in_restriction <- rowSums(in_restriction) > 0
      
      # Return NA if 100% of observations are removed...
      if (all(in_restriction)){
        return(list('data' = NA))
      }
      
      report_restriction <- round(mean(in_restriction) * 100, 2)
      
      
      if (verbose){
        message(paste0('For choice "', p, '", ', report_restriction, '% of observations ',
                       'removed from the empirical distribution of AMCE\n', 'because of randomization restrictions.'))
      }
      bd <- baseline_design[!in_restriction,,drop=F]
      return(list('data' = bd))
    })
    names(baseline_data) <- unique_choice
    
    if (any(sapply(baseline_data, is.na))){
      return(NA)
    }
  }
  return(baseline_data)
}