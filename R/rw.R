##' Read and write cns file
##'
##' \itemize{
##'   \item \code{read.cns()}: Read in cns file generated from \code{CNVkit} as a \code{data.frame}.
##'   \item \code{write.cns()}: Write standard cns format file as table.
##' }
##'
##' @title Standard read in cns file
##' @param cnspath The path of cns file.
##' @return
##' \itemize{
##'   \item \code{read.cns()}: A \code{data.frame}.
##'   \item \code{write.cns()}: A table separated \code{txt} format file.
##' }
##'
##' @examples
##' require('magrittr')
##'
##' cnsFile <- system.file('extdata', 'example.filter.cns', package = 'CNVanno')
##' cnsMat <- cnsFile %>% read.cns
##' \dontrun{
##' ## write cns file
##' write.cns(cnsMat, 'cnsMat.txt')
##' }
##'
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @importFrom utils read.table
##' @rdname rwcns
##' @export
##'
read.cns <- function(cnspath) {
  cnsin <- read.table(cnspath,
                      header = TRUE,
                      sep = '\t',
                      stringsAsFactors = FALSE)

  return(cnsin)
}

##' @rdname rwcns
##' @param cns The standard cns format.
##' @inheritParams read.cns
##' @importFrom utils write.table
##' @export
##'
write.cns <- function(cns, cnspath) {
  write.table(cns,
              cnspath,
              sep = '\t',
              row.names = FALSE)
}

## cnsFile <- '/home/Yulong/RESEARCH/ciphergene/CNVscript/CNVanno/inst/extdata/example.filter.cns'
