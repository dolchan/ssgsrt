#' Read a file in MSigDB/gmt format into a list
#'
#' This function reads a GMT file and creates a list of gene sets.
#'
#' @param gmtfile a text file in 'gmt' format.
#'        See 'parse_a_line_in_gmt' function for detail.
#' @param min.size The minimum gene set size to read.  Defaults to 0.
#' @param max.size The maximum gene set size to read.  Defaults to no limit.
#' @param simplify Defaults to TRUE.  See 'output' for detail.
#' @return a list of gene set.
#'        If TRUE, the return structure is simplified to a list of genesets
#'        each of which is a list of genes and the name of each geneset is
#'        the name of gene set from MSigDB file, the first element of each line.
#'        If FALSE, it returns a list of three element list each element of which
#'        is \code{name}, \code{url}, and \code{geneset}.
#' @export
#' @examples
#' read_MSigDB('c2.cp.biocarta.v4.0.symbols.gmt.txt')

read_MSigDB <- function(gmtfile, min.size = 0, max.size = -1, simplify = TRUE) {
  gstemp <- read.table(gmtfile, sep="\n")
  genesets <- apply(gstemp, 1, parse_a_line_in_gmt)

  list.name = array("", length(genesets))

  # convert to named list -- this will make the access to the list much easier.
  for (ii in 1:length(genesets)) {
    list.name[ii] = genesets[[ii]]$name
  }
  names(genesets) = list.name

  # Prune out ones that are too long/short
  for (ii in length(genesets):1) {
    if (length(genesets[[ii]]$geneset) < min.size ||
        (max.size > 0 && length(genesets[[ii]]$geneset) > max.size))
      genesets[[ii]] <- NULL
  }

  if (simplify) {
    genesets <- sapply(genesets, "[", "geneset")
    names(genesets) <- sub(".geneset$", "", names(genesets))
  }

  return(genesets)
}
