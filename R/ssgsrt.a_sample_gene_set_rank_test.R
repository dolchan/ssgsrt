#' Perform a rank-based single sample gene set analysis.
#'
#' This function performs a rank-based single sample gene set analysis for
#' a sample (\code{x}) and a list of gene sets (\code{geneset_list}).
#'
#' @param x a sample in a named array format.
#' @param geneset_list a named list of gene sets, each of list element is named
#'                     and a list ofgenes.
#' @param alternative \code{c('less','greater','two.sided')}.
#'                    'less' ('greater') to check if the ranks of given gene set
#'                    is higher (lower) than the rest.
#' @param test.method \code{c('ks.test', 'wilcox.test')}.  Defaults to 'ks'.
#' @export
#' @examples
#' a_sample_gene_set_rank_test(x, a_geneset_list)
a_sample_gene_set_rank_test <- function(x, geneset_list, alternative = "two.sided", test.method = "ks") {
  xr <- rank(x)
  names(xr) <- names(x)

  gsa_score <- vector("list", length(geneset_list))

  for (ii in 1:length(geneset_list)) {
    c1 <- xr[geneset_list[[ii]]]
    c2 <- setdiff(xr, c1)

    test.methods <- c("ks.test", "wilcox.test")
    if (pmatch(test.method, test.methods) == 1) {
      test.stat <- ks.test(c1, c2, alternative = alternative)
    }
    else {
      test.stat <- wilcox.test(c2, c1, alternative = alternative)
    }

    # to prevent assigning Inf to z.stat
    # if (test.stat$p.value == 1)
    #   test.stat$p.value <- 1 - .Machine$double.eps
    # if (test.stat$p.value == 0)
    #   test.stat$p.value <- .Machine$double.eps

    gsa_score[[ii]] <-
      c(z.stat = qnorm(test.stat$p.value),
        p.value = test.stat$p.value,
        test.stat$statistic)
  }

  names(gsa_score) <- names(geneset_list)
  gsa_score
}
