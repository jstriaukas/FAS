[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/FAS)](https://cran.r-project.org/package=FAS)
[![Downloads](http://cranlogs.r-pkg.org/badges/grand-total/FAS)](https://cran.rstudio.com/web/packages/FAS/index.html) 
[![Downloads](http://cranlogs.r-pkg.org/badges/FAS)](http://www.r-pkg.org/pkg/FAS)

# FAS

FAS - Inference methods for Factor-Augmented Sparse Regression 

## About

The *FAS* package implements *inference* methods for high-dimensional inference on sparse regression coefficient. In addition,  the package is equipped with the functionality for adaptively selecting the LASSO regression tuning parameter for linear and factor-augmented sparse regression. 

The approach is based on adapetively bootstrapping the effective noise of the LASSO regression, for more details see [^1][^2]. 

## Run to install the package

```{r }
# CRAN version - 1.0.0
install.packages(FAS) 

# Development version - 1.0.0
# install.packages("FAS")
library(devtools)
install_github("jstriaukas/FAS")
```

## Acknowledgements

Jad Beyhum gratefully acknowledges financial support from the Research Fund KU Leuven through the grant STG/23/014. Jonas Striaukas gratefully acknowledges the financial support from the European Commission, MSCA-2022-PF Individual Fellowship, Project 101103508. Project Acronym: MACROML.

## References

[^1]: Lederer, J., & Vogt, M. (2021). Estimating the Lasso's effective noise. The Journal of Machine Learning Research, 22(1), 12658-12689. https://www.jmlr.org/papers/volume22/20-539/20-539.pdf

[^2]: Beyhum, J., & Striaukas, J. (2023). Tuning-free testing of factor regression against factor-augmented sparse alternatives. arXiv preprint arXiv:2307.13364.

