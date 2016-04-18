#' Extract a specific component from gene set analysis (\code{gsa_list}) and
#' create a matrix.
#'
#' @param gsa_list a list of gene set analysis from \code{gene_set_rank_test}.
#' @param which a name of a list element to extract. Defaults to 'z.stat'.
#' @return a matrix (M gene sets x N samples) of selected element.
#' @export

get_a_value <- function(gsa_list, which = 'z.stat') {
  n_row <- length(gsa_list[[1]])
  n_col <- length(gsa_list)

  gsa_z_stat <- matrix(nrow = n_row, ncol = n_col)

  for (ii in 1:n_col) {
    gsa_z_stat[, ii] <- sapply(gsa_list[[ii]], '[', which)
  }

  colnames(gsa_z_stat) <- names(gsa_list)
  rownames(gsa_z_stat) <- names(gsa_list[[1]])

  gsa_z_stat
}

#' Extract \code{'z.stat'} element
#' @export

get_z_stat <- function(gsa_list) {
  get_a_value(gsa_list, 'z.stat')
}

#' Extract \code{'p.value'} element
#' @export

get_p_value <- function(gsa_list) {
  get_a_value(gsa_list, 'p.value')
}

#' Extract \code{'D'} element
#' @export

get_D_stat <- function(gsa_list) {
  get_a_value(gsa_list, 'D')
}
