# is2 - Importance Sampling Squared <img src="man/figures/hexlogo_small.png" align="right"/> #

Importance Sampling Squared is a package designed to take the output of a `pmwg` Particle Metropolis within Gibbs (a set of model parameter and random effect estimates) and calculate the marginal likelihood of the model using Importance Sampling.

Marginal Likelihoods are the current gold standard for model selection and comparison and are typically used in calculating Bayes factors. This package implements the methods outlined in [Robustly estimating the marginal likelihood for cognitive models via importance sampling (Tran et. al)](https://link.springer.com/article/10.3758/s13428-020-01348-w). The IS2 method is a robust estimation method that accounts for model flexibility and provides unbiased estimates of the marginal likelihood. Marginal likelihood estimates allow you to assess the fit of a model and the model's flexibility by integrating the likelihood across the prior predictive space of a model. In hierarchical models, obtaining the marginal likelihood is difficult, as the likelihood function is the density of the data with the random effects integrated out when viewed as a function of the group-level parameters; an integral which is often unavailable (computationally or as it is intractable). Despite this integral being intractable, IS2 allows a method of estimating the marginal likelihood when the likelihood is intractable but can be estimated in unbiasedly. 

The method works by first sampling from the posterior via a sampling scheme such as MCMC (or here, PMwG). These draws are then used to create the importance distribution for the fixed parameters. This importance distribution is constructed by fitting a mixture of normal or Student t distributions to these MCMC samples. We then construct conditional proposal parameters - called particles - for each subject. The marginal likelihood is then estimated unbiasedly which is combined with the importance distribution. From this method, the importance sampling procedure is in itself an importance sampling procedure which can be used to estimate the likelihood. 

## Installation

To install the latest stable version you can use the following command to install from CRAN:

`install.packages("is2")`

If you want the (possibly unstable) development version, you can also install the package using devtools as follows:

`devtools::install_github('newcastlecl/is2', ref="develop")`

This package is tested and should work on all versions of R > 4.0.1, however instructions on installing to an earlier version of R are included below.

## Using the package

The `pmwg` package documentation available at the bookdown site https://newcastlecl.github.io/samplerDoc/ includes some worked examples using the `is2` package. In particular [Chapter 5](https://newcastlecl.github.io/samplerDoc/estimating-the-marginal-likelihood-via-importance-sampling-is2.html) includes a detailed example you can follow along to.

Also available online is the package documentation at https://newcastlecl.github.io/is2/ which consists of this README, a Reference of help documentation for individual functions, a list of changes to the project over time and more.
