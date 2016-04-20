#' Perform a rank-based single sample gene set analysis.
#'
#' This function performs a rank-based single sample gene set analysis for
#' a set of samples (\code{data}) and a list of gene sets (\code{geneset_list}).
#'
#' @param data Either a set of samples in \code{data.frame} format
#'             or a sample in a named array format.
#' @param geneset_list a named list of gene sets, each of list element is named
#'                     and a list ofgenes.
#' @param alternative see 'a_sample_gene_set_rank_test' for detail.
#' @param test.method \code{c('ks.test', 'wilcox.test')}.  Defaults to 'ks'.
#' @export
#' @examples
#' gene_set_rank_test(data, a_geneset_list)

gene_set_rank_test <- function(data, geneset_list, row_names = NULL, alternative = "two.sided", test.method = "ks") {
  if (is.data.frame(data)) {
    n_samples <- ncol(data)
    gsas_list <- vector("list", length(n_samples))
    row_names_x <- row.names(data)
    for (ii in 1:n_samples) {
      x <- data[, ii]
      a_gsa <- a_sample_gene_set_rank_test(x, geneset_list,
                                           row_names = row_names_x,
                                           alternative = alternative,
                                           test.method = test.method)
      gsas_list[[ii]] <- a_gsa
      print(sprintf("%s (%d / %d) processed", colnames(data)[ii], ii, n_samples))
    }
    names(gsas_list) <- colnames(data)
    gsas_list
  } else {
    a_sample_gene_set_rank_test(data, geneset_list,
                                row.names = row_names,
                                alternative = alternative,
                                test.method = test.method)
  }
}

