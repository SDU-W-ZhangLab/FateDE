# FateDE

**The recent advancement of single-cell RNA sequencing technologies is leading to a more complete understanding of lineage specification during cell development. Meanwhile, many mathematical models of gene regulatory networks have been proposed to link molecular regulatory mechanisms with observed cell state. However, it remains challenging that neither a data-only nor a theory-only approach can be considered sufficient for unraveling the underlying genetic programs of cell-fate decision events. Here, we present a framework that attempts to steer the cell fate determinant identification by exploring the dynamic behavior of stochastic differentiation equation model (SDE). We first tracked the time-evolution to a toggle switch circuit system with double-negative feedback which is known to control binary developmental decision, and further quantify the dynamic behavior of the functional circuit. Then, we model the network basing on analyzing gene pair co-expression and co-variation patterns, we developed a computational method to predict putative cell fate determinants, using the pseudo-time of single-cell sequencing (scRNA-seq) during differentiation.** 

插入流程图

## Installation

You can install the development version of FateDE from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Biores410/FateDE")
library(FateDE)
```


## Steps

1.  Data preprocessing (if required) can be implemented using the Seurat or other package.
2.  Infer trajectory  (if required)  can be implemented using the Monocle or other package.
3.  Fitting can be implemented using function *fateDE_fit_one()*.
4.  *Determine_transition_state()* function can be used to identify pre-somatic cell cells and differentiated cells
5.  FateDE can be implemented using function *FateDE_main()*.



## Examples

- Analysis for CNS  simulation data : *[FateDE_CNS_simulation.R](https://github.com/Biores410/FateDE/blob/master/examples/FateDE_CNS_simulation.R)*

- Analysis for crest data : *[FateDE_crest.R](https://github.com/Biores410/FateDE/blob/master/examples/FateDE_crest.R)*

- Analysis for hepatoblast data : *[FateDE_hepatoblast.R](https://github.com/Biores410/FateDE/blob/master/examples/FateDE_hepatoblast.R)*

The data used can be loaded using the following code:

```R
library(FateDE)

data(simulation) # CNS simulation
data(crest) # crest
data(hepatoblast) # hepatoblast
```

