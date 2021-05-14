# is2 - Importance Sampling Squared

Importance Sampling Squared is a package designed to take the output of a `pmwg` Particle Metropolis within Gibbs (a set of model parameter and random effect estimates) and calculate the marginal likelihood of the model using Importance Sampling.

## Installation

To install the latest stable version you can use the following command to install from CRAN:

`install.packages("is2")`

If you want the (possibly unstable) development version, you can also install the package using devtools as follows:

`devtools::install_github('newcastlecl/is2', ref="develop")`

This package is tested and should work on all versions of R > 4.0.1, however instructions on installing to an earlier version of R are included below.

## Using the package

The `pmwg` package documentation available at the bookdown site https://newcastlecl.github.io/samplerDoc/ includes some worked examples using the `is2` package. In particular [Chapter 5](https://newcastlecl.github.io/samplerDoc/estimating-the-marginal-likelihood-via-importance-sampling-is2.html) includes a detailed example you can follow along to.

Also available online is the package documentation at https://newcastlecl.github.io/is2/ which consists of this README, a Reference of help documentation for individual functions, a list of changes to the project over time and more.
