## Initial Submission

This corrects issues with an initial submission of `FactorHet` to CRAN. I have
made the following changes per the comments from CRAN:

1. DESCRIPTION: I have used usethis::use_tidy_description() to remove excess
white spaces from the "Description"

2. Do not use T/F: I have corrected all instances of T/F as noted in the initial
feedback.

3. All of the exported functions now contain an argument for "Value", possibly
referencing the discussion in "Details" if necessary.

4. `print`: I have removed all uses of `print` except for `print.FactorHet` and
`summary.Factorhet`. They have all been swapped for `message`, `warning` or
`stop`. The only exception is for `print(g)` where `g` is an object from
`ggplot2`. This is equivalent to `plot` but allows for `list(g,g)` to be
rendered correctly.

5. Modifying the global environment: The function in `print_functions.R` that
modifies the global environment has been modified to remove this behavior.

## R CMD check results

There were no ERRORs or WARNINGs. The only NOTE was about an initial submission. 

## Downstream dependencies

There are no downstream dependencies.