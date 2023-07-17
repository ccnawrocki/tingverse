#' Get differential expression
#' 
#' This function returns a table with differential expression data for genes in 
#' a counts table.
#' 
#' @param counts.table A normalized and quality controlled counts table. Class: 
#' data frame.
#' @param meta.data A data frame of sample names and groups to perform 
#' differential expression on. The row names must exist as column names of 
#' counts.table.
#' @param meta.feature The name of the column in meta.data for which the groups 
#' being compared exist.
#' @param value.1 One of the two labeled groups being compared. This will be a 
#' level of the meta.feature column.
#' @param value.2 The second of the labeled groups being compared. This will be 
#' a level of the meta.feature column.
#' @param test The statistical test that will be run to compare the two groups. 
#' This test must be "wilcox" or "t" and nothing else.
#' @export
diff.exp <- function(counts.table, meta.data, meta.feature, value.1, value.2, test = "wilcox") {
  indexed_meta_1 <- meta.data[meta.data[[meta.feature]] == value.1,]
  indexed_meta_2 <- meta.data[meta.data[[meta.feature]] == value.2,]
  indexed_counts_1 <- counts.table[,rownames(indexed_meta_1)] 
  indexed_counts_2 <- counts.table[,rownames(indexed_meta_2)]
  means_1 <- rowMeans(indexed_counts_1,na.rm = T) 
  means_2 <- rowMeans(indexed_counts_2,na.rm = T) 
  diff_exp <- data.frame(rownames(counts.table),means_1,means_2)
  colnames(diff_exp) <- c("Gene.ID",value.1,value.2)
  diff_exp$log2.FC <- log2(diff_exp[[value.1]] / diff_exp[[value.2]])
  diff_exp_ls <- list()
  for (gene in rownames(diff_exp)) {
    diff_exp_ls[[gene]] <- list()
    diff_exp_ls[[gene]][[value.1]] <- as.numeric(indexed_counts_1[gene,])
    diff_exp_ls[[gene]][[value.2]] <- as.numeric(indexed_counts_2[gene,])
  }
  do_test <- function(ls_of_rows,test_name) {
    if (test_name == "wilcox") {
      p <- wilcox.test(ls_of_rows[[1]],ls_of_rows[[2]],exact=F)$p.value
      return(p)
    }
    if (test_name == "t") {
      p <- t.test(ls_of_rows[[1]],ls_of_rows[[2]])$p.value
      return(p)
    }
    else {
      print("NOT A SUPPORTED TEST.")
      p <- NA
      return(p)
    }
  }
  p_vals <- unlist(mclapply(diff_exp_ls,do_test,test))
  p_adjs <- p.adjust(p_vals,method = "fdr")
  diff_exp$p.val <- p_vals
  diff_exp$p.adj <- p_adjs
  return(diff_exp)
}
