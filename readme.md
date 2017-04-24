# HMSC

The Hierarchical Modelling of Species Communities (HMSC) framework described in Ovaskainen et al. (2017). This framework uses Bayesian hierarchical modelling to account for environment, traits and phylogeny to model species communities. In addition, this framework also include and spatially (or temporally) autocorrelated latent variables to measure association among species. This R package implements all aspects of Ovaskainen et al. (2017) but is a work in progress that aims at going beyond the ideas presented in Ovaskainen et al. (2017).

This version of the R package is in development.

## Installing HMSC

To install this development version of HMSC you can use the ```install_github``` function in the [```devtools```](http://cran.r-project.org/web/packages/devtools/index.html) package, like this:

```{r}
# install devtools if you haven't already
# install.packages('devtools')

# load the package
library(devtools)

# install HMSC from github
install_github('HMSC', 'guiblanchet')

# and load it
library(HMSC)
```

This installing is somewhat long, it take roughly 15 minutes, do not despair. 

Note that many parts of the R package have been implemented using ```Rcpp```, ```RcppArmadillo``` and C++11. This has been known to cause some installation issues especially under Mac OS X. Essentially, these problems arise because of C/C++ and Fortran compiler related issues. So far, these problems have been solved by following the explanation presented in section 2.10 and 2.16 of http://dirk.eddelbuettel.com/code/rcpp/Rcpp-FAQ.pdf.

More information on this issue is given here: http://thecoatlessprofessor.com/programming/rcpp-rcpparmadillo-and-os-x-mavericks-lgfortran-and-lquadmath-error/

## Reporting bugs

If you find a bug in the code or have suggestions for improvements, please let me know via [the issues reporting system](https://github.com/guiblanchet/HMSC/issues)

# Reference

Ovaskainen, O., G. Tikhonov, A. Norberg, F. Guillaume Blanchet, L. Duan, D. Dunson, T. Roslin, and N. Abrego. 2017. How to make more out of community data? A conceptual framework and its implementation as models and software. Ecology Letters 20:561–576.
