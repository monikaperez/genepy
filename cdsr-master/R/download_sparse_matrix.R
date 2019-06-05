#' download raw sparse matrix data files from taiga to R
#'
#' @param data_name the name of the taiga dataset.
#' @param version_number the version number of the taiga dataset 
#' @param file_name the name of the taiga file
#'
#' @return a sparse matrix of type S4
#'
#' @export download_sparse_matrix

download_sparse_matrix <- function(data_name, version_number, file_name) {
  require(taigr)
  require(Matrix)

  data <- Matrix::readMM(taigr::download.raw.from.taiga(data.name=data_name, data.version=version_number, data.file=file_name))
  
  return(data)
}
