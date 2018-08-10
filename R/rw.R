##' Read and write cnvkit cns file
##'
##' \itemize{
##'   \item \code{read.cnvkit()}: Read in cns file generated from \code{CNVkit} as a \code{data.frame}.
##'   \item \code{write.cnvkit()}: Write standard cnv format file as table.
##' }
##'
##' @title Standard read and write multiple cnv files
##' @param cnvpath The path of cnv file.
##' @return
##' \itemize{
##'   \item \code{read.cnvkit()}: A \code{data.frame}.
##'   \item \code{write.cnvkit()}: A table separated \code{txt} format file.
##' }
##'
##' @examples
##' require('magrittr')
##'
##' cnsFile <- system.file('extdata', 'example.cnvkit', package = 'CNVanno')
##' cnsMat <- cnsFile %>% read.cnvkit
##' \dontrun{
##' ## write cns file
##' write.cnvkit(cnsMat, 'cnsMat.txt')
##' }
##'
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @importFrom utils read.table
##' @rdname rwcnv
##' @export
##'
read.cnvkit <- function(cnvpath) {
  cnsin <- read.table(cnvpath,
                      header = TRUE,
                      sep = '\t',
                      stringsAsFactors = FALSE)

  return(cnsin)
}

##' @rdname rwcnv
##' @param cns The standard cns format.
##' @inheritParams read.cnvkit
##' @importFrom utils write.table
##' @export
##'
write.cnv <- function(cnv, cnspath) {
  write.table(cnv,
              cnspath,
              sep = '\t',
              row.names = FALSE)
}

## cnsFile <- '/home/Yulong/RESEARCH/ciphergene/CNVscript/CNVanno/inst/extdata/example.cnvkit'
