#' Parse a line in 'gmt' format
#'
#' This function parses a line of string in 'gmt' format,
#'   for example a line in a MSigDB file.
#'
#' @param aLine a string in 'gmt' format. There are assumed three main elements.
#'        Each element is 'tab' delitmeted.  The first item is a name, the second
#'        item is 'url' or something like that, and the last item is a list of genes,
#'        separated by 'space'.
#' @return a list of three elements: \code{name}, \code{url}, \code{geneset}
#' @keywords gmt parse
#' @export
#' @examples
#' parse_a_line_in_gmt('a_name\tan_url\tgene1 gene2 gene3')
#' @author Seungchan Kim (dolchan@gmail.com)
#'

parse_a_line_in_gmt <- function(aLine) {
  aLine <- strsplit(as.character(aLine), "\\s")[[1]]

  return(list(name = aLine[1], url = aLine[2], geneset = aLine[-c(1:2)]))
}
