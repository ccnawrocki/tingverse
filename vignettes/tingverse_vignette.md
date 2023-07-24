# tingverse vignette

Required packages:
```
library(tingverse)
library(tidyverse)
library(org.Hs.eg.db)
library(clusterProfiler)
library(parallel)
library(ggrepel)
library(knitr)
```

## Read in Example Data 
This is old data from a Nanostring GeoMx experiment. `geom` is a normalized counts table. `metadata` is the meta data for each AOI in the experiment. scRNA-seq data could be used as well. You'll need to use your own data.

```
geom <- read.delim("/Users/cnawrocki/Library/CloudStorage/OneDrive-UniversityofVirginia/Ting Lab/Summer 2022/COLE/HCC/GeoMx/Data/qNormCounts.tsv")
metadata <- read.csv("/Users/cnawrocki/Library/CloudStorage/OneDrive-UniversityofVirginia/Ting Lab/Summer 2022/COLE/HCC/GeoMx/Data/meta_small.csv")

# Some basic manipulations to make the data files compatible. 
colnames(geom) <- substr(colnames(geom),14,23)
metadata$Sample_ID_short<-gsub('-','.',metadata$Sample_ID_short)
geom <- geom[,colnames(geom) %in% metadata$Sample_ID_short]
rownames(metadata)<-metadata$Sample_ID_short
metadata <- metadata[metadata$Sample_ID_short %in% colnames(geom),]
```

## Differential Expression
Use the `differential.expression()` function first. Then, use the `volcano.plot()` function to visualize the results. FUse `?differential.expression` and `?volcano.plot` to view documentation.

```
DEA <- differential.expression(counts.table = geom,meta.data = metadata,meta.feature = "Tissue_Type",value.1 = "Tumor",value.2 = "Vessel",test = "wilcox")
volcano.plot(DEA,color.vals = c("green3","grey","red3"))
```
![volcano](https://github.com/ccnawrocki/tingverse/assets/68296470/7d78fe5d-a36a-404c-8193-fc0ee96fb1a5)


## GSEA
`?easy.gsea` will show documentation.

```
GSEA_obj <- easy.gsea(DEA,order.stat = "log2.FC",ontology = "BP",organism = "human")

# Viewing the results
kable(head(as.data.frame(GSEA_obj[,1:10])),format = 'markdown')
```
<img width="807" alt="table" src="https://github.com/ccnawrocki/tingverse/assets/68296470/dee148e5-d799-4973-8694-540c38b6b572">


#### Looking at a Specific Gene Set
We can subset for only gene sets with 50 or more genes in them. Then, we can output the GSEA plot for the gene set with the highest enrichment score. `gseaplot` is a `ClusterProfiler` function.

```
results_oi <- GSEA_obj@result[order(GSEA_obj@result$enrichmentScore,decreasing = T) & GSEA_obj@result$setSize > 50,]
gseaplot(GSEA_obj,geneSetID = results_oi$ID[1])
```
![GSEA](https://github.com/ccnawrocki/tingverse/assets/68296470/66ece6df-7043-4685-965e-1a59028eff21)


#### Dotplot Based on Our Results
`dotplot` is a `ClusterProfiler` function. Here, we will look at the top 40 gene sets and separate based on if they are activated or deactivated.

```
dotplot(GSEA_obj,showCategory=40,split=".sign") + facet_grid(.~.sign) + theme(axis.text.y = element_text(size = 2))
```
![dotplot](https://github.com/ccnawrocki/tingverse/assets/68296470/fe9d5ab9-3d2d-4536-b0ca-a212adb1b2e6)


#### Heatplot Based on Our Results
`heatplot` is a `ClusterProfiler` function. Here, we will look at the top 50 gene sets and visualize the LFC for each gene.

```
heatplot(GSEA_obj,showCategory = 50,foldChange=GSEA_obj@geneList) + theme(axis.text.y = element_text(size = 4),axis.text.x = element_text(size=2,angle=60))
```
![heatplot](https://github.com/ccnawrocki/tingverse/assets/68296470/675bc94f-ad4d-48c4-9f73-878a125983cb)


## Enrichment Analysis and Barplot
If we do not want to look at GSEA results, we can utilize this tool. We can look at all genes with a p.adj that is significant, run the analysis on them, then visualize the top 30 gene sets.

```
genes_oi <- rownames(DEA[DEA$p.adj < 0.05,])
EA_obj <- enrichGO(gene = genes_oi,OrgDb = "org.Hs.eg.db",keyType = "SYMBOL",ont = "BP")
barplot(EA_obj, showCategory = 30) + theme(axis.text.y = element_text(size = 4))
```
![barplot](https://github.com/ccnawrocki/tingverse/assets/68296470/9f3d93cc-8b81-43b7-bdcc-6112e167c55f)


