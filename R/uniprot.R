#' Using r code to connect uniprot web api
#'
#' You can choose one id as well as another id as output id
#' Read uniprot documentation to know more information
#'
#' @param query vector of protein ids
#' @param inputid type of input id, character
#' @param outputid type of output id, character
#' @param fmt output format
#'
#' @return a dataframe
#' @export
#'
#' @examples
#' idMapping(query=proid, inputid="ACC", outputid="P_ENTREZGENEID", mft="fmt")
#'
idMapping <- function(query, inputid, outputid, fmt){
  query = paste(query, collapse=",")
  r=httr::POST('http://www.uniprot.org/uploadlists/',
               body=list(from=inputid, to=outputid, format=fmt,  query=query), encode="form")
  cont=httr::content(r, type="text")
  result=readr::read_tsv(cont)
}
