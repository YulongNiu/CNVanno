##' Read and write cnvkit CNV file
##'
##' \itemize{
##'   \item \code{read.cnvkit()}: Read in CNV generated from \code{CNVkit} as a \code{data.frame}.
##'   \item \code{write.cnvkit()}: Write standard CNV format file as table.
##'   \item \code{read.cnvnator()}: Read in CNV generated from \code{CNVnator} (version 0.3.3) as a \code{data.frame}.
##' }
##'
##' @title Standard read and write multiple CNV files
##' @param cnvpath The path of cnv file.
##' @return
##' \itemize{
##'   \item \code{read.cnvkit()}: A \code{data.frame}.
##'   \item \code{read.cnvkit()}: A \code{data.frame}.
##'   \item \code{write.cnvkit()}: A table separated \code{txt} format file.
##' }
##'
##' @examples
##' require('magrittr')
##'
##' cns <- system.file('extdata', 'example.cnvkit', package = 'CNVanno')
##' cnsMat <- cns %>% read.cnvkit
##'
##' nator <- system.file('extdata', 'example.cnvnator', package = 'CNVanno')
##' natorMat <- nator %>% read.cnvnator
##' \dontrun{
##' ## write CNV file
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
##' @inheritParams read.cnvkit
##' @importFrom magrittr  %>%
##' @importFrom stringr str_extract str_sub
##' @export
##'
read.cnvnator <- function(cnvpath) {
  cnsin <- read.table(cnvpath,
                      header = FALSE,
                      sep = '\t',
                      stringsAsFactors = FALSE)

  ## separate cnv
  chr <- cnsin[, 2] %>%
    str_extract('^.*:') %>%
    str_sub(1, -2)

  start <- cnsin[, 2] %>%
    str_extract(':\\d+') %>%
    str_sub(2, -1) %>%
    as.numeric

  end <- cnsin[, 2] %>%
    str_extract('-\\d+') %>%
    str_sub(2, -1) %>%
    as.numeric

  natorRes <- data.frame(chromosome = chr,
                         start = start,
                         end = end,
                         stringsAsFactors = FALSE)

  natorRes %<>% cbind.data.frame(cnsin[, c(1, 4:9)],
                                 stringsAsFactors = FALSE)
  colnames(natorRes)[-1:-3] <- c('type', 'normalized_RD', 'e-val1', 'e-val2', 'e-val3', 'e-val4', 'q0')

  return(natorRes)
}


##' @rdname rwcnv
##' @param cnv The standard CNV format.
##' @param savepath Save path of output CNV files.
##' @inheritParams read.cnvkit
##' @importFrom utils write.table
##' @export
##'
write.cnvkit <- function(cnv, savepath) {
  write.table(cnv,
              savepath,
              sep = '\t',
              row.names = FALSE)
}

## cnsFile <- '/home/Yulong/RESEARCH/ciphergene/CNVscript/CNVanno/inst/extdata/example.cnvkit'
