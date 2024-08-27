# Principal Nested Simplices

This repository contains the code for replication of the results in the paper "Principal Nested Simplices" by anonymous authors. 

The sub-directory `PSA` contains the main functions of Principal Nested Simplices (PSA).

The sub-directory `Functions` contains R functions that are used to produce the figures in the paper.

The sub-directory `Simulation` contains the materials for the simulation study provided in Section 5. The materials include the code for data generation, the generated data, code for analysis of the data, and the figures.

The sub-directory `Diatom` contains the materials for the analysis of the Diatom data set provided in Section 6. The materials include the data set, code for analysis of the data, and the figures.

## Installation

```{r}
source('PSA_init.R')
```

## Description

The function `psa` provides an implementation of PSA via Subsimplies (PSA-S) and PSA via Suborthants (PSA-O). The output includes scores, loading vectors, lower dimensional compositional representation, and the corresponding vertices of the lower dimensional simplices.

The function `compare_analysis` is a wrapper function for applying two PSA methods and three benchmark methods (Principal Component Analysis (PCA), Power Transformed PCA, Log-ratio PCA) to the same data set. `compare_analysis` provides the results of the methods in the same format.

The function `plotdendrogram2` produces a dendrogram representation of PSA-S or PSA-O.

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

# produce dendrogram representations
plotdendrogram2(X.psas)
plotdendrogram2(X.psao)
```
