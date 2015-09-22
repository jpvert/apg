# APG: Accelerated Proximal Gradient in R

APG is an R package to solve convex but potentially non-smooth optimization problems with an accelerated proximal gradient.

It implements basic functions for least-squares and logistic regression with lasso, elastic net, group lasso or monotonicity constraints, and allows the user to easily implement extensions.

To install APG from R, type:

```{r}
install.packages('devtools')
library(devtools)
install_github("jpvert/apg")
```

Check the [vignette](vignettes/apg_vignette.Rmd) for examples.
