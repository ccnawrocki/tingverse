#' Easily perform GSEA
#'
#' This function performs GSEA without the headache of getting files in specific 
#' formats that the original GSEA function requires. It also eliminates the many 
#' parameters the original GSEA function requires by simply setting them as the 
#' values that will be most often used by the Ting Lab. One should revert to the 
#' original GSEA function, if trying to perform GSEA with atypical parameters.
#'
#' @param normdata A normalized and quality-controlled counts table with the 
#' first column being Gene IDs. Class: data frame.
#' @param metadata A meta data table that has an ID column with identical IDs to
#' the column names of normdata. Class: data frame.
#' @param directory The working directory to send output files to. Class: string.
#' @param feature The name of the column in metadata that contains the two 
#' groups for which differential expression should be performed. Class: string.
#' @param lab.1 One label in the feature column that is to be compared to lab.2.
#' Class: string.
#' @param lab.2 One label in the feature column that is to be compared to lab.1.
#' Class: string.
#' @param ID The name of the ID column in metadata. Class: string.
#' @param gene_set The file path to a gene set in .gmt format. This file is 
#' downloadable from the GSEA website. Class: string.
#' @param num_perm The number of permutations used to obtain the p value. 1000 
#' is recommended. Class = integer.
#' @export
EasyGSEA <- function(normdata,metadata,directory,feature,lab.1,lab.2,ID,gene_set,num_perm){
  
  normdata <- cbind('NAME'=normdata[,1],'DESCRIPTION'=NA,normdata[,2:ncol(normdata)])
  
  ID_col <- which(colnames(metadata) == ID)
  feature_col <- which(colnames(metadata) == feature)
  subsetted_meta <- metadata[,c(ID_col,feature_col)]
  colnames(subsetted_meta) <- c('Sample_ID',feature)
  
  phenotypes_list <- data.frame('Sample_ID'=colnames(normdata[,3:length(normdata)])) %>% 
    inner_join(subsetted_meta, by='Sample_ID')
  phenotypes_list$phenotype <- 0
  for (group in unique(phenotypes_list[,feature])) {
    phenotypes_list$phenotype[phenotypes_list[[feature]]==group]<-group
  }
  
  phenotypes_list <- phenotypes_list[order(phenotypes_list[[feature]]),]
  phenotypes_list <- phenotypes_list[phenotypes_list[[feature]]==lab.1|phenotypes_list[[feature]]==lab.2,]
  
  normdata <- normdata[,c('NAME','DESCRIPTION',phenotypes_list[,'Sample_ID'])]
  print(normdata)
  
  phens <- paste(phenotypes_list$phenotype[1])
  for (i in 2:length(phenotypes_list$phenotype)) {
    phens <- paste(phens,phenotypes_list$phenotype[i],sep=' ')
  }
  
  setwd(directory)
  fileConn<-file('phenos.cls')
  writeLines(c(
    paste(nrow(phenotypes_list),n_distinct(phenotypes_list$phenotype),1,sep=' '),
    paste('#',lab.1,lab.2,sep=' '),
    phens
  ),fileConn)
  close(fileConn)
  
  
  GSEA::GSEA(input.ds = normdata,
             input.cls = paste(directory,'phenos.cls',sep='/'),
             gs.db = gene_set,
             output.directory = directory, 
             doc.string = paste('GSEA',lab.1,'vs',lab.2,sep='-'),
             reshuffling.type = "sample.labels",
             nperm = num_perm, weighted.score.type =  1,
             nom.p.val.threshold = -1,
             fwer.p.val.threshold = -1,
             fdr.q.val.threshold   = 0.25,
             topgs = 20,
             adjust.FDR.q.val = F,gs.size.threshold.min = 15,
             gs.size.threshold.max = 500,
             reverse.sign = F, 
             preproc.type = 0,
             random.seed = 09132001,
             perm.type = 0,
             fraction = 1.0,
             replace = F,
             save.intermediate.results = F,
             use.fast.enrichment.routine = T
  )
  
}