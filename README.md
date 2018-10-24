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
hello_world()
```

## Citation

If you use `htspop`, please specify the version and cite:

> Ribeiro-dos-Santos, AM, de Souza, SJ (2018) Htspop: high-troughput sequencing
> population genetic functions.

## Contact

André M. Ribeiro-dos-Santos, andremrsantos@gmail.com