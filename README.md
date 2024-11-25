# Principal Subsimplex Analysis

This repository contains the code for replication of the results in the paper "Principal Subsimplex Analysis" by Hyeon Lee, Kassel Liam Hingee, Janice L. Scealy, Andrew T. A. Wood, Eric Grunsky, and J. S. Marron (2024+). This repository also contains a package "PSA".

The sub-directory `PSA` contains the contents for the package "PSA".

The sub-directory `utils` contains additional R functions that are used to produce the figures in the paper.

The sub-directory `Simulation` contains the materials for the simulation study provided in Section 5. The materials include the code for data generation, the generated data, code for analysis of the data, and the figures.

The sub-directory `Diatom` contains the materials for the analysis of the Diatom data set provided in Section 6. The materials include the data set, code for analysis of the data, and the figures.

The sub-directory `Figures` contains the figures in the other sections. Some of the figures are generated by "figures_generator.R".

## Installation

To install the package "PSA", please use the following code.
```{r}
library(devtools)
devtools::install_github("haneone33/Principal-Subsimplex-Analysis/PSA")
```

## Description

The main function of "PSA" is `psa`. It provides an implementation of PSA via Subsimplies (PSA-S) and PSA via Suborthants (PSA-O). See **Example** for the usage. The output includes scores, loading vectors, lower dimensional compositional representation, and the corresponding vertices of the lower dimensional simplices.

The function `compare_analysis` is used to compare the results of PSA with those of three benchmark methods (Principal Component Analysis (PCA), Power Transformed PCA, Log-ratio PCA). `compare_analysis` provides the results of the methods in the same format.

## Example

```{r}
# sample data generation
set.seed(1)
X = to_simplex(matrix(rnorm(0, 1, 30), 6, 5))

# apply PSA-S and PSA-O only
X.psas = psa('s', X)
X.psao = psa('o', X)

# apply PSA and benchmark methods at the same time
X.res = compare_analysis(X)
```
