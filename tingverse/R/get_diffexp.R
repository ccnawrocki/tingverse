#' Get differential expression
#' 
#' This function returns a table with differential expression data for genes in 
#' a counts table.
#' 
#' @param df A normalized and quality controlled counts table. Class: data frame.
#' @param de A data frame of sample names and groups to perform DE on.
#' @export
get_diffexp<-function(df,de){
  # df is a counts table
  # de is a data frame of sample names and groups to perform DE on
  typ<- as.factor(de[,2])
  diffexp <- foreach(cl = levels(typ), .combine=cbind) %do% {
    idx<-typ==cl
    idx[is.na(idx)] <- FALSE
    rowMeans(df[,idx],na.rm=T)
  } # Calculate group means.
  diffexp <- as.data.frame(diffexp)
  colnames(diffexp) <- levels(typ)
  test<-combn(1:ncol(diffexp),2)
  for (n in 1:ncol(test)){
    cn<-colnames(diffexp)
    g1<-test[1,n]
    g2<-test[2,n]
    t1<-levels(typ)[g1]
    t2<-levels(typ)[g2]
    
    diffexp$LFC<-log2(diffexp[,g1] / diffexp[,g2])
    diffexp$p <- foreach(i = 1:nrow(df), .combine=c) %do% {
      wilcox.test(as.numeric(df[i, typ==t1]), as.numeric(df[i, typ==t2]), 
                  exact=FALSE)$p.value
    } # Calculate Wilcoxon Rank Sum p-values
    diffexp$padj<-(p.adjust(diffexp$p,method='fdr'))
    colnames(diffexp)<-c(cn,paste0('LFC_',t1,'_',t2),paste0('p_',t1,'_',t2),
                         paste0('padj_',t1,'_',t2))
  }  
  return(diffexp)
}
