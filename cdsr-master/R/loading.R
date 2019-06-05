#' Apply arxspan to ccle conversion to matrix or list of matrices (assuming cell line names are in the rownames)
#'
#' @param dat: matrix or list of matrices
#'
#' @return matrix or list of matrices with rownames converted to CCLE ID
#' @export
#'
#' @examples
map_arxspan_to_ccle <- function(dat) {
  stopifnot(is.data.frame(dat) | is.matrix(dat) | is.list(dat))
  if (is.data.frame(dat) | is.matrix(dat)) {
    if (all(str_sub(rownames(dat), 1, 4) == 'ACH-')) {
      rownames(dat) <- celllinemapr::arxspan.to.ccle(rownames(dat))
    }
  }
  if (is.list(dat)) {
    dat <- llply(dat, function(df) {
      if (all(str_sub(rownames(df), 1, 4) == 'ACH-')) {
        rownames(df) <- celllinemapr::arxspan.to.ccle(rownames(df))
      }
      return(df)
    })
  }
  return(dat)
}


#' Extract hugo symbol column names from data matrices using "HGNC (ENTREZ)" format
#'
#' @param dat: matrix or list of matrices
#'
#' @return
#' @export
#'
#' @examples
extract_hugo_symbol_colnames <- function(dat) {
  stopifnot(is.data.frame(dat) | is.matrix(dat) | is.list(dat))
  if (is.data.frame(dat) | is.matrix(dat)) {
    if (all(!is.na(stringr::str_match(colnames(dat), ' \\([0-9\\&]+\\)$')[,1]))) {
      colnames(dat) <- str_match(colnames(dat), '^(.+) \\([0-9\\&]+\\)$')[,2]
    }
    dat <- dat[, !is.na(colnames(dat)) & colnames(dat) != ''] #any columns that didn't have a Hugo symbol are dropped here
  }
  if (is.list(dat)) {
    dat <- llply(dat, function(df) {
      if (all(!is.na(stringr::str_match(colnames(df), ' \\([0-9\\&]+\\)$')[,1]))) {
        colnames(df) <- str_match(colnames(df), '^(.+) \\([0-9\\&]+\\)$')[,2]
      }
      df <- df[, !is.na(colnames(df)) & colnames(df) != '', drop = FALSE] #any columns that didn't have a Hugo symbol are dropped here
      return(df)
    })
  }
  return(dat)
}
