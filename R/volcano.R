#' Make a volcano plot to show differential expression
#' 
#' This function can be used to make a volcano plot showing which genes are 
#' differentially expressed across two groups of AOIs or cells. Differential
#' expression is calculated using a wilcox rank sum test for difference in mean 
#' expression across the two groups.
#' 
#' @param normdata A normalized and quality-controlled counts table with the 
#' column names being Sample IDs or Cell IDs and the row names being gene IDs. 
#' Class: data frame.
#' @param metadata A meta data table that has an ID column with identical IDs to
#' the column names of normdata. Class: data frame.
#' @param feature The name of the column in metadata that contains the two 
#' groups for which differential expression should be performed. Class: string.
#' @param lab.1 A label in the feature column that is to be compared to lab.2.
#' Class: string. Note: needs to be before lab.2 alphabetically.
#' @param lab.2 A label in the feature column that is to be compared to lab.1.
#' Class: string. Note: needs to be after lab.1 alphabetically.
#' @param ID The name of the ID column in metadata. Class: string.
#' @param myCols A vector of three colors to use in the plot. Default: blue, 
#' black, and red.
#' @export
volcano <- function(normdata,metadata,feature,lab.1,lab.2,ID,
                    myCols=c('blue','black','red')) {
  
  idx<-c(grep(lab.1,metadata[[feature]]),grep(lab.2,metadata[[feature]]))
  df<-normdata[,idx]
  de<-data.frame(metadata[,ID],metadata[,feature])[idx,]
  
  dea <- get_diffexp(df,de)
  colnames(dea)[3:5]<-c('LFC','p','padj')
  
  myCol<-myCols
  names(myCol)<-c(lab.1,'Neither',lab.2)
  dea$delabel<-'Neither'
  dea$delabel[dea[,5]<.05 & dea[,3]>1]<-lab.1
  dea$delabel[dea[,5]<.05 & dea[,3]< -1]<-lab.2
  dea$lbl<-NA
  dea$lbl[dea$delabel != 'Neither'] <- rownames(dea)[dea$delabel != 'Neither']
  tmp<-ggplot(data=dea,aes(x=LFC, y=-log10(padj),col=delabel, label=lbl)) +
    geom_point() + theme_minimal() + geom_text_repel() +
    labs(title=paste(lab.1,'vs.',lab.2,sep=' ')) + 
    scale_color_manual(values=myCol) +
    geom_vline(xintercept=c(-1,1),col="gray") + 
    geom_hline(yintercept=-log10(.05),col="gray")
  print(tmp)
  
}