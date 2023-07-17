# tingverse

<img width="385" alt="Screenshot 2023-07-17 at 2 37 42 PM" src="https://github.com/ccnawrocki/tingverse/assets/68296470/96c23198-2781-42e1-b0ef-f1e7a81c3541">


The tingverse is a set of basic functions for Bioinformatics at Ting Lab of the Massachusetts General Hospital Cancer Center. It was originally developed by a Cole Nawrocki in 2022, a summer intern aiming to speed up analyses. Much of the principles and code that are integral to the tingverse were taught to Cole by Mike Raabe. 
In 2023, the tingverse was updated extensively by Cole, resulting in version 2.0.

## How To Install
You can install and load the tingverse package by running the following lines in R:

```
install.packages('devtools')
library(devtools)
devtools::install_github('ccnawrocki/tingverse')
library(tingverse)
```

### Dependencies
The tingverse has multiple dependencies. R makes updating dependencies very painful, so, for simplicity's sake, install and load each dependency as follows in order to use the tingverse.

### GSEA Dependencies
The GSEA functionality of the tingverse is a wrapper for clusterProfiler functions that depend on certain gene set databases available for download through Bioconductor. These two databases are about 80 MB each, which is not trivial. Thus, if you do not plan to use the GSEA functionality of the tingverse, then do not download these databases. Otherwise, you can download these databases as follows: 

```
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("org.Mm.eg.db")
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db
```
Once you have downloaded the clusterProfiler package and these databases, you will have access to the tingverse's full capabilities.
