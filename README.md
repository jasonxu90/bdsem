bdsem
=====

EM algorithm for inference in discretely observed birth-shift-death processes

## Description
This package includes functions for inference in discretely observed birth-shift-death processes. All routines are based on methods that apply more generally to multi-type Markov branching processes; our representation of the birth-death-shift process is a specific two-type case of such a process. See [1] for details on the model specification and methodology.

Capabilities include functions that compute finite-time transition probabilities and restricted sufficient statistics of the process, based on spectral probability generating function techniques using a Fast Fourier Transform. These routines allow for the evaluation of the discrete data observed likelihood, leading to maximum likelihood inference for discretely observed multi-type branching processes given any black box optimization routine. Furthermore, inference via the Expectation-Maximization (EM) algorithm is fully implemented, which often outperforms direct optimization of the likelihood in practice. The package also includes functions to forward simulate realizations of the birth-death-shift process or other multi-type branching process via the Gillespie algorithm [2].


## Installation
The package can be installed directly from github using the `devtools` package, which can easily be installed using the command `install.packages("devtools")`.
See https://github.com/hadley/devtools for more details.

To install `bdsem`, run the following:
```r
library(devtools)
install_github("jasonxu90/bdsem")
```

## Vignette
By default, the vignette is not compiled during the installation. We provide a vignette, ``BDSEM: inferring MLE rates of discretely observed birth-death-shift processes'', that walks through simulation and inference functions necessary to recreate all simulation studies in [1]. Although examples are smaller-scale, they still require quite some computing time. Users are encouraged to run all functions in the vignette locally by reducing the number of simulations, initial restarts, etc when necessary. 

The source `.Rnw` file, as well as a `.pdf` of the vignette compiled using `knitr`, are included in the repository. To automatically compile the vignette upon package installation, instead install using the line
`install_github("jasonxu90/bdsem", build_vignettes = TRUE)`. Note that this option will be significantly slower than installing the package alone.


## References
1. Xu, J., Guttorp, P., Kato‐Maeda, M., & Minin, V. N. (2015). Likelihood‐based inference for discretely observed birth–death‐shift processes, with applications to evolution of mobile genetic elements. Biometrics, 71(4), 1009-1021.

2.  Gillespie DT (1977) "Exact stochastic simulation of coupled chemical reactions," *Journal of Physical Chemistry,* 81(25):2340-2361.
