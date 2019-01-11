# htspop

## Overview

`htspop` is a R package designed for analysis and computation of population
genetic statistics using high-dimensional biallelic data. It includes the
following statistics:

- *Reich's f2, f3, and f4* statistics, similar to treemix implementation.
- *Patterson's D* statistics (also known *ABBA/BABA* test).
- *jackknife* mean estimator.
- *Reich, Weir & Cockerham, Hudson, and Wright Fst* and their bootstrap
  estimation.
- *Nei's standard* and *Da* genetic distance and a bootstrap estimator.

## Instalation

You may install the development version, using `devtools`

```R
# devtools instalation
devtools::install_github("andremrsantos/htspop")
```

## Usage

```R
## Simulate genotype matrix
geno <- matrix(sample(0:2, 100, replace = TRUE), ncol = 10)
## Convert into Allele Count structure
ac <- ac_matrix_from_genotype(geno) 
## Run F Statistics
jackknife(f4_stat(ac, c(1, 2, 3, 4)))
jackknife(f3_stat(ac, c(1, 2, 3)))
## Run D Statistics
jackknife(d_stat(ac, c(1, 2, 3, 4)))
## Compute Weir & Cockerham Fst and Nei Standard Genetic Distance
fst(ac, "wc")
nei(ac)
```

## Citation

If you use `htspop`, please specify the version and cite:

> Ribeiro-dos-Santos, AM, de Souza, SJ (2018) Htspop: high-troughput sequencing
> population genetic functions.

## Contact

Create an issue to report bugs, propose new functions or ask for help. Please take in consideration this project is under development.

André M. Ribeiro-dos-Santos, andremrsantos@gmail.com
