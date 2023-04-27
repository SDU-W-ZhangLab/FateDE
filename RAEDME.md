# FateDI

**The recent advancement of single-cell RNA sequencing technologies is leading to a more complete understanding of lineage specification during cell development. Meanwhile, many mathematical models of gene regulatory networks have been proposed to link molecular regulatory mechanisms with observed cell state. However, it remains challenging that neither a data-only nor a theory-only approach can be considered sufficient for unraveling the underlying genetic programs of cell-fate decision events. Here, we present a framework that attempts to steer the cell fate determinant identification by exploring the dynamic behavior of stochastic differentiation equation model (SDE). We first tracked the time-evolution to a toggle switch circuit system with double-negative feedback which is known to control binary developmental decision, and further quantify the dynamic behavior of the functional circuit. Then, we model the network basing on analyzing gene pair co-expression and co-variation patterns, we developed a computational method to predict putative cell fate determinants, using the pseudo-time of single-cell sequencing (scRNA-seq) during differentiation.** 

插入流程图

## Contents

- [Functions](#Functions)

- [Steps](#Steps)
- [Examples](#Examples)



## Functions

First, you need to load the files in the *functions* folder：

```R
source('utils.r')
library("reticulate")
pd <- import("pandas")
source_python("base.py")
source('plot.r')
```

- **base.py：including relevant codes for Gaussian processes**

  - ```R
    # Searching for gene sets that increase first and then decrease later
    function find_genes(df1, df2, G1, G2, fate1, fate2)
    ```

    | Arguments  |                                                              |
    | ---------- | ------------------------------------------------------------ |
    | **df1**    | a dataframe,  cells by row and features by columns. The first column is pseudotime and the second column represents the state of the cell, the rest are listed as genes. |
    | **df2**    | the same as df1                                              |
    | **G1**     | a vector, destiny specific genes representing df1, i.e. a set of genes that continuously increase in expression over time in df1. |
    | **G2**     | the same as G1                                               |
    | **fate1**  | a character string, the Name of Fate.                        |
    | **fate2**  | the same as fate1                                            |
    | **return** | a dataframe, each gene's score, including *FSV, P, Q*        |

- **utils.r：the functions of Correlation Analysis**

  - ```R
    # Fitting gene expression of one branch 
    function fateDE_fit_one(df_notfit)
    ```

    | Arguments     |                                                |
    | ------------- | ---------------------------------------------- |
    | **df_notfit** | a dataframe, df1_notfit Branch 1 unfitted data |
    | **return**    | a dataframe of fitting                         |

  - ```R
    # Fitting gene expression of two branches
    function fateDE_fit_two(df1_notfit,df2_notfit) 
    ```

    | Arguments                  |                            |
    | -------------------------- | -------------------------- |
    | **df1_notfit, df2_notfit** | a dataframe, unfitted data |
    | **return**                 | a list of fitting          |

  - ```R
    # Calculating Score
    function FateDE_cal_Score(Cartesian_Product,df1,df2,State=c('State 1','State 2','State 3'))
    ```

    | Arguments             |                                      |
    | --------------------- | ------------------------------------ |
    | **Cartesian_Product** | Cartesian product matrix             |
    | **df1, df1**          | fitting data                         |
    | **State**             | The state to which each cell belongs |
    | **return**            | Score for each gene pair             |

- **plot.r：drawing related functions**

  - *plot_Genes_Heatmap* : draw a thermal diagram
  - *plot_SingleGene_trend* : curve of a single gene over pseudo time on two branches
  - *plot_DoubleGenes_trend* : a two-dimensional map of genes

## Steps

1.  Data preprocessing (if required) can be implemented using the Seurat or other package.
2.  Infer trajectory  (if required)  can be implemented using the Monocle or other package.
3.  Fitting can be implemented using function *fateDE_fit_one()* or *fateDE_fit_two()*.
4.  Find destiny specific genes can be implemented using cluster or SSlogis.
5.  Searching for gene sets that increase first and then decrease later can be implemented using function *find_genes()*.
6.  Correlation calculation can be implemented using function *FateDE_cal_Score()*.



## Examples

- Analysis for crest data : *[FateDE_crest.r]()*

- Analysis for hepatoblast data : *[FateDE_hepatoblast.r]()*



FateDE supports Mac OS X and Windows. 