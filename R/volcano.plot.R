#' Make a volcano plot to show differential expression
#' 
#' This function can be used to make a volcano plot showing which genes are 
#' differentially expressed across two groups of AOIs or cells. Differential
#' expression is calculated using a wilcox rank sum test or a student's t test 
#' for difference in mean expression across the two groups.
#' 
#' @param diff.exp.results The data frame that is outputted by the diff.exp 
#' function.
#' @param color.vals A vector of three color values. The first color will 
#' correspond to genes up-regulated in the first group being compared in one's 
#' differential expression analysis. The second color will correspond to genes 
#' up-regulated in neither group. The third color will correspond to genes 
#' up-regulated in the third group. 
#' @param log2.FC.cutoff The LFC cutoff for genes to be labeled on the plot.
#' @param p.adj.cutoff The p.adj cutoff for genes to be labeled on the plot.
#' @param label.res Higher values mean more overlaps are allowed for gene 
#' labeling purposes on the plot.
#' @export
volcano.plot <- function(diff.exp.results,color.vals=c('blue','black','red'),
                    log2.FC.cutoff=1,p.adj.cutoff=0.05,label.res=50) {
  
  names(color.vals)<-c(colnames(diff.exp.results)[2],'Neither',colnames(diff.exp.results)[3])
  diff.exp.results$delabel<-'Neither'
  diff.exp.results$delabel[diff.exp.results[,5]<p.adj.cutoff & diff.exp.results[,4]>log2.FC.cutoff]<-colnames(diff.exp.results)[2]
  diff.exp.results$delabel[diff.exp.results[,5]<p.adj.cutoff & diff.exp.results[,4]< -log2.FC.cutoff]<-colnames(diff.exp.results)[3]
  diff.exp.results$lbl<-NA
  diff.exp.results$lbl[diff.exp.results$delabel != 'Neither'] <- rownames(diff.exp.results)[diff.exp.results$delabel != 'Neither']
  v.plot<-ggplot(data=diff.exp.results,aes(x=log2.FC, y=-log10(p.adj),col=delabel, label=lbl)) +
    geom_point() + theme_minimal() + geom_text_repel(max.overlaps = label.res) +
    labs(title=paste(colnames(diff.exp.results)[2],'vs.',colnames(diff.exp.results)[3],sep=' ')) + 
    scale_color_manual(values=color.vals) +
    geom_vline(xintercept=c(-log2.FC.cutoff,log2.FC.cutoff),col="gray") + 
    geom_hline(yintercept=-log10(p.adj.cutoff),col="gray")
  print(v.plot)
  
}
