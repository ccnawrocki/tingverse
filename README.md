# tingverse

<img width="385" alt="Screenshot 2023-07-17 at 2 37 42 PM" src="https://github.com/ccnawrocki/tingverse/assets/68296470/96c23198-2781-42e1-b0ef-f1e7a81c3541">

## Disclaimer
`tingverse` is a set of __basic__ functions for Bioinformatics at the Ting Lab of the Massachusetts General Hospital Cancer Center. It was one of my projects when I was a __summer intern__ in the Ting Lab. I worked on it during the summer of 2022 and the summer of 2023. It showcases what I learned over the course of those two summers, and it is __no longer updated__, nor is it representative of current best practices that I employ in my work. 

### Acknowledgements
Mike Raabe and Peter Richieri are the two people who mentored me during these two summers. Much of the code in this repository was influenced by their previous work.

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

```
install.packages("tidyverse")
install.packages("ggrepel")
install.packages("parallel")
library(tidyverse)
library(ggrepel)
library(parallel)
```

### GSEA Dependencies
The GSEA functionality of the tingverse is a wrapper for clusterProfiler functions that depend on certain gene set databases available for download through Bioconductor. These two databases are about 80 MB each, which is not trivial. Thus, if you do not plan to use the GSEA functionality of the tingverse, then do not download these databases. Otherwise, you can download these databases as follows: 

```
install.packages("BiocManager")
library(BiocManager)
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("org.Mm.eg.db")
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db
```
Once you have downloaded the clusterProfiler package and these databases, you will have access to the tingverse's full capabilities.
