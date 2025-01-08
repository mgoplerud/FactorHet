# FactorHet [![R-CMD-check](https://github.com/mgoplerud/FactorHet/workflows/R-CMD-check/badge.svg)](https://github.com/mgoplerud/FactorHet/actions) [![codecov](https://codecov.io/gh/mgoplerud/FactorHet/branch/master/graph/badge.svg?token=R0RDGBF1OX)](https://app.codecov.io/gh/mgoplerud/FactorHet)

This package estimates heterogeneous effects in factorial and conjoint experiments. Its details are fully described in [Goplerud, Imai, and Pashley (2025)](https://arxiv.org/abs/2201.01357): "Estimating Heterogeneous Causal Effects of High-Dimensional Treatments: Application to Conjoint Analysis". 

The core method is a Bayesian regularized finite mixture-of-experts where moderators can affect an individual's probability of cluster membership and a sparsity-inducing prior fuses together levels of each factor in each cluster while respecting ANOVA-style sum-to-zero constraints described in [Egami and Imai (2019)](https://imai.fas.harvard.edu/research/files/int.pdf). The posterior mode is found using an AECM algorithm with a number of techniques to accelerate convergence. Approximate quantification of uncertainty is provided by examining the Hessian of the log-posterior. Additional details are explained in the paper and (briefly) in the package documentation.

It can be installed from CRAN or the most-to-update version can be installed using `devtools`.

```
# CRAN
install.packages("FactorHet")
# Up-to-Date GitHub Version
library(devtools)
devtools::install_github('mgoplerud/FactorHet')
```

There are two key functions for estimating the model: In most cases, one will prefer to use the `FactorHet_mbo` function to jointly (i) estimate the amount of regularization by minimizing a criterion such as the BIC using model-based optimization and (ii) estimate the final model. However, if one has a specific value of `lambda` of interest, one can fit the model for a fixed amount of regularization using `FactorHet`. A simple example is shown below:

```
fit_FH <- FactorHet_mbo(
  formula = y ~ factor_1 + factor_2 + factor_1 : factor_2, 
  design = design,
  moderator = ~ moderator_1 + moderator_2)
```
In the case of repeated observations, the individual is specified via `group` and the task identifier is specified via `task`. In the case of a conjoint experiment, the profile identifier (i.e. "left" or "right") is specified via `choice_order`. An example is shown below:

```
fit_FH <- FactorHet_mbo(
  formula = y ~ factor_1 + factor_2 + factor_1 : factor_2, 
  design = design, moderator = ~ moderator_1 + moderator_2, 
  group = ~ id, task = ~ task, choice_order = ~ choice_left)
```

Finally, after fitting the model, there are functions to calculate the Average Marginal Effect (AME) and related concepts (e.g. ACE, AMIE). A simple example is shown below:

```
AME(fit_FH)
```

The effects of moderators on cluster membership can be analyzed using two key functions; first, `posterior_by_moderators` shows the estimated distribution of (posterior) cluster membership probabilities by covariates. Second, `margeff_moderators` shows the change in the prior cluster membership as one moderator changes, averaging across all other moderators. This is similar to a marginal effect in a multinomial logistic regression. Example code is shown below:

```
posterior_by_moderators(fit_FH)
margeff_moderators(fit_FH)
```