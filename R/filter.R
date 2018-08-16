##' @include AllClasses.R AllGenerics.R
NULL


##' Filter the CNV file generated from \code{CNVkit} and \code{CNVnator}
##'
##' For the \code{CNVkit}
##' \itemize{
##'   \item 1. filter the \code{cns} with cn number.
##'   \item 2. filter low sequence depth.
##'   \item 3. filter X/Y chromosomes
##' }
##'
##' For the \code{CNVnator}
##' \itemize{
##'   \item 1. filter X/Y chromosomes
##' }
##'
##' @title CNV filter
##' @param cngain A \code{numeric} value as the gain copy number threshold, and the default value is 2.
##' @param cnloss A \code{numeric} value as the loss copy number threshold, and the default value is 2.
##' @param dep A \code{numeric} value as the sequence depth, and the default value is 0.01
##' @param sexchrom A \code{logic} value whether filter sex chromosomes (X and Y), the default value is \code{TRUE}.
##' @param rawkit The standard \code{CNVkit} output as an \code{RawCNV} object.
##' @return A filtered \code{RawCNV} object.
##' @examples
##' require('magrittr')
##'
##' kitf <- system.file('extdata', 'example.cnvkit', package = 'CNVanno') %>% read_cnvkit %>% filter_cnvkit
##'
##' natorf <- system.file('extdata', 'example.cnvnator', package = 'CNVanno') %>% read_cnvnator %>% filter_cnvnator
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @importFrom dplyr transmute filter
##' @importFrom magrittr %>% %<>%
##' @rdname filterraw
##' @export
##'
filter_cnvkit <- function(rawkit, cngain = 2, cnloss = 2, dep = 0.01, sexchrom = TRUE) {

  ## step1: filter gain and loss cnv
  ## step2: filter low depth

  l <- rawkit@params %>%
    transmute((cn > cngain | cn < cnloss) & depth > dep) %>%
    unlist

  ## step3: filter sex chromosomes
  if (sexchrom) {
    lsex <- rawkit@rawCNV %>%
      transmute(chromosome != 'chrX' & chromosome != 'chrY') %>%
      unlist
    l %<>% `&`(lsex)
  } else {}

  cnvkitf <- new('RawCNV',
                 rawCNV = filter(rawkit@rawCNV, l),
                 params = filter(rawkit@params, l),
                 method = rawkit@method)

  return(cnvkitf)

}


##' @param rawnator The standard \code{CNVnator} output as an \code{RawCNV} object.
##' @inheritParams filter_cnvkit
##' @importFrom magrittr  %>%
##' @importFrom dplyr select mutate
##' @rdname filterraw
##' @export
##'
filter_cnvnator <- function(rawnator, sexchrom = TRUE) {

  ## step1: filter sex chromosomes
  if (sexchrom) {
    l <- rawnator@rawCNV %>%
      transmute(chromosome != 'chrX' & chromosome != 'chrY') %>%
      unlist
  } else {}

  cnvnatorf <- new('RawCNV',
                   rawCNV = filter(rawnator@rawCNV, l),
                   params = filter(rawnator@params, l),
                   method = rawnator@method)

  return(cnvnatorf)
}



##' Filtering \code{CoreCNV}.
##'
##' Filer the \code{CoreCNV} according black lists.
##'
##' @title Filter \code{CoreCNV}
##' @inheritParams Filter
##' @return A \code{CoreCNV} object.
##' @examples
##' library('magrittr')
##'
##' nator <- system.file('extdata', 'exampleseg.cnvnator', package = 'CNVanno') %>%
##'   read_natorkit %>%
##'   filter_natorkit %>%
##'   Segment(natorf, interlen = 10L)
##'
##' natorf <- Filter(kitf, **, overlaprate = 0.5)
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @importFrom magrittr %>%
##' @importFrom dplyr mutate
##' @rdname Filter-methods
##' @exportMethod Filter
##'
## setMethod(f = 'Filter',
##           signature = c(core = 'CoreCNV', blacklist = 'tbl_df', overlaprate = 'double'),
##           definition = function(core, blacklist, overlaprate = 0.5, ...) {

##           })

