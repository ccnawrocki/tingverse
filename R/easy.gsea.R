#' Do easy GSEA
#' 
#' This function is a wrapper for the gseGO function of the clusterProfiler
#' package. It implements the function in a simpler way.
#' 
#' @param de.table The output of the diff.exp function. A table containing 
#' differential expression results.
#' @param order.stat The statistic in de.table that will be used to order the 
#' differentially expressed genes. Typically, one will use "log2.FC" as this 
#' parameter's value.
#' @param ontology Either "BP", "MP", or "CC". BP is for biological 
#' processes, MP is for molecular function, and CC is for cellular component.
#' @param organism Either "human" or "mouse". Based on the value, the function 
#' will reference the corresponding database of genes.
#' @export
easy.gsea <- function(de.table,order.stat,ontology,organism) {
  if (organism == "human") {
    organism.DB <- "org.Hs.eg.db"
  }
  if (organism == "mouse") {
    organism.DB <- "org.Mm.eg.db"
  }
  ordered_de_table <- de.table[order(de.table[[order.stat]],decreasing = T),]
  ordered_de_stat <- ordered_de_table[[order.stat]]
  names(ordered_de_stat) <- rownames(ordered_de_table)
  gsea <- gseGO(ordered_de_stat,ontology,OrgDb = organism.DB,keyType = "SYMBOL",eps = 1e-300)
  return(gsea)
}
