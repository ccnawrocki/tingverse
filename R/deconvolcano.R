#' Make a deconvolcano plot
#' 
#' This function prints a volcano plot that depicts differential containment of 
#' immune cell types across two groups of AOIs in Nanostring GeoMx data.
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
#' @param n_cell_proxy The name of the column in metadata that describes a 
#' variable that can be used as a proxy for the number of cells in each AOI. 
#' Often this is 'AOINucleiCount' but not always. Class: string.
#' @param neg_probe_names The names of the negative probes as they appear in the 
#' row names of normdata. Class: vector.
#' @export
deconvolcano <- function(normdata,metadata,feature,lab.1,lab.2,ID,n_cell_proxy,neg_probe_names) {
  
  idx<-c(grep(lab.1,metadata[[feature]]),grep(lab.2,metadata[[feature]]))
  exp<-normdata[,idx] %>% as.matrix()
  bg<-derive_GeoMx_background(norm=exp,probepool=rep(1,nrow(exp)),
                              negnames=neg_probe_names)
  dec<-spatialdecon(norm=exp,bg=bg)
  
  metadata <- metadata[idx,c(ID,n_cell_proxy,feature)]
  colnames(metadata) <- c('Sample_ID',n_cell_proxy,feature)
  
  immune_props <- as.data.frame(dec$prop_of_all) %>% t() %>% 
    cbind( 'Sample_ID' = colnames(dec$prop_of_all)) %>% as.data.frame()
  
  metadata <- inner_join(metadata,immune_props,by='Sample_ID')
  
  for (celltype in colnames(metadata[,4:21])) {
    new_col <- as.numeric(unlist(metadata[,n_cell_proxy])) * as.numeric(unlist(metadata[,celltype]))
    metadata[,celltype] <- new_col
  }
  
  immune_counts <- group_by(metadata, .data[[feature]]) %>% 
    summarise(N=sum(.data[[n_cell_proxy]])) %>% as.data.frame() 
  
  for (celltype in colnames(metadata[,4:21])) {
    counts <- tapply(metadata[[celltype]],metadata[,feature],sum)
    immune_counts[,celltype] <- floor(counts)
  }
  
  p <- c()
  LFC <- c()
  for (celltype in colnames(immune_counts[,3:20])) {
    TEST<-prop.test(immune_counts[,celltype],immune_counts$N,
                    alternative='two.sided',correct=F)
    p <- c(p,TEST$p.value)
    LFC <- c(LFC, log2((immune_counts[1,celltype]/immune_counts[1,'N']) / (immune_counts[2,celltype]/immune_counts[2,'N'])))
  }
  
  decon_dea <- column_to_rownames(immune_counts,var=feature)[,2:19]
  decon_dea <- data.frame(lab.1=as.numeric(decon_dea[1,])/immune_counts[1,'N'],
                          lab.2=as.numeric(decon_dea[2,]/immune_counts[2,'N']))
  rownames(decon_dea) <- colnames(immune_counts[,3:20])
  decon_dea <- cbind(decon_dea,'LFC'=LFC,'p.value'=p)
  
  myCol<-c('blue','black','red')
  names(myCol)<-c(lab.1,'Neither',lab.2)
  decon_dea$delabel<-'Neither'
  decon_dea$delabel[decon_dea[,4]<.05 & decon_dea[,3]>1]<-lab.1
  decon_dea$delabel[decon_dea[,4]<.05 & decon_dea[,3]< -1]<-lab.2
  decon_dea$lbl<-NA
  decon_dea$lbl[decon_dea$delabel != 'Neither'] <- rownames(decon_dea)[decon_dea$delabel != 'Neither']
  
  deconvolcano_plot <- ggplot(data=decon_dea,aes(x=LFC, y=-log10(p),col=delabel, label=lbl)) +
    geom_point() + theme_minimal() + geom_text_repel() +
    labs(title=paste(lab.1,'vs',lab.2,sep='_')) + scale_color_manual(values=myCol) +
    geom_vline(xintercept=c(-1,1),col="gray") + 
    geom_hline(yintercept=-log10(.05),col="gray")
  print(deconvolcano_plot)
  
}
