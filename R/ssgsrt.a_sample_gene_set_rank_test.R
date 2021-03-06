#' Perform a rank-based single sample gene set analysis.
#'
#' This function performs a rank-based single sample gene set analysis for
#' a sample (\code{x}) and a list of gene sets (\code{geneset_list}).
#'
#' @param x a sample in a named array format.
#' @param geneset_list a named list of gene sets, each of list element is named
#'                     and a list ofgenes.
#' @param row_names a list of (gene) names of x. Useful when names(x) is NULL.
#' @param alternative \code{c('less','greater','two.sided')}.
#'                    If 'less'/'greater', test if the overall gene expression of given gene set
#'                    is lower/higher than the rest.  If 'two.sided', test either lower or higher.
#'                    Defaults to 'two.sided'.
#' @param test.method \code{c('ks.test', 'wilcox.test')}.  Defaults to 'ks'.
#' @return a list of z.stat, D (ks.test) or W (wilcox.test) stat, p.value, and adj.p.value (FDR) for each gene set.
#'         z.stat indicates if the gene set expression is "lower" (negative) or "higher" (positive).
#'         We use p-value estimated from the actual statistical test (ks.test or wilcox.test) and
#'         use normal distribution to compute "pseudo" z statistics.
#' @export
#' @examples
#' a_sample_gene_set_rank_test(x, a_geneset_list)
a_sample_gene_set_rank_test <- function(x, geneset_list,
                                        row_names = NULL,
                                        alternative = 'two.sided', test.method = "ks") {

  gsa_score <- vector("list", length(geneset_list))

  if (!is.null(row_names)) {
    if (length(row_names) < length(x))
      print("length(row_names) should be equal to length(x)")
    else
      names(x) <- row_names
  }

  for (ii in 1:length(geneset_list)) {
    c1_idx <- match(geneset_list[[ii]], names(x))
    c1_idx <- c1_idx[!is.na(c1_idx)]
    c1 <- x[c1_idx]
    c2 <- x[-c1_idx]

    test.methods <- c("ks.test", "wilcox.test")
    method.choice <- pmatch(test.method, test.methods)
    if (method.choice == 1) {
      test.stat <- ks.test(c2, c1, alternative = alternative)
    }
    else {
      if (method.choice == 2) {
        test.stat <- wilcox.test(c1, c2, alternative = alternative)
      }
      else {
        stop(sprintf('test.method should be either "ks.test" or "wilcox.test"'))
      }
    }

    # to prevent assigning Inf to z.stat
    # if (test.stat$p.value == 1)
    #   test.stat$p.value <- 1 - .Machine$double.eps
    # if (test.stat$p.value == 0)
    #   test.stat$p.value <- .Machine$double.eps

    # now to compute "z.stat" which is to indicate if the gene set expression is
    #   "lower" (negative) or "higher" (positive)
    # we use p-value estimated from the actual statistical test and use normal distribution
    #   to compute "pseudo" z statistics.
    alternatives <- c("less", "greater", "two.sided")
    alternative.choice = pmatch(alternative, alternatives)

    if (alternative.choice == 1) {
      z.stat <- qnorm(test.stat$p.value)
    }
    else {
      if (alternative.choice == 2) {
        z.stat <- qnorm(1-test.stat$p.value)
      }
      else {
        if (alternative.choice == 3) {
          z.stat <- qnorm(1-test.stat$p.value) * sign(mean(c1)-mean(c2))
        }
        else {
          stop(sprintf('alternative should be one of "less", "greater", or "two.sided"'))
        }
      }
    }

    gsa_score[[ii]] <-
      c(z.stat = z.stat,
        test.stat$statistic,
        p.value = test.stat$p.value)
  }

  # p-value correction for multiple testing across gene sets
  pvals <- get_p_value(list(gsa_score), adjusted = FALSE)
  pvals.adj <- p.adjust(pvals, method="BH")

  for (ii in 1:length(gsa_score)) {
    gsa_score[[ii]] <- c(gsa_score[[ii]], adj.p.value = pvals.adj[ii])
  }

  names(gsa_score) <- names(geneset_list)
  gsa_score
}
