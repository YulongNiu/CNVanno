##' @include AllClasses.R
NULL


##' Filter the CNV file generated from \code{CNVkit}
##'
##' For the \code{CNVkit}
##' \itemize{
##'   \item 1. filter the \code{cns} with cn number.
##'   \item 2. filter low sequence depth.
##'   \item 3. filter X/Y chromosomes
##' }
##'
##'
##' @title CNV filter
##' @param cngain A \code{numeric} value as the gain copy number threshold, and the default value is 2.
##' @param cnloss A \code{numeric} value as the loss copy number threshold, and the default value is 2.
##' @param dep A \code{numeric} value as the sequence depth, and the default value is 0.01
##' @param sexchrom A \code{logic} value whether filter sex chromosomes (X and Y), the default value is \code{TRUE}.
##' @param RawCNVkit The standard \code{CNVkit} output as an \code{RawCNV} object.
##' @return A filtered \code{RawCNV} object.
##' @examples
##' require('magrittr')
##'
##' kitf <- system.file('extdata', 'example.cnvkit', package = 'CNVanno') %>% read_cnvkit %>% filter_cnvkit
##'
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @importFrom dplyr transmute filter
##' @importFrom magrittr %>% %<>%
##' @export
##'
filter_cnvkit <- function(RawCNVkit, cngain = 2, cnloss = 2, dep = 0.01, sexchrom = TRUE) {

  ## step1: filter gain and loss cnv
  ## step2: filter low depth

  l <- RawCNVkit@params %>%
    transmute((cn > cngain | cn < cnloss) & depth > dep) %>%
    unlist

  ## step3: filter sex chromosomes
  if (sexchrom) {
    lsex <- RawCNVkit@rawCNV %>%
      transmute(chromosome != 'chrX' & chromosome != 'chrY') %>%
      unlist
    l %<>% `|`(lsex)
  } else {}

  cnvkitf <- new('RawCNV',
                 rawCNV = filter(RawCNVkit@rawCNV, l),
                 params = filter(RawCNVkit@params, l),
                 method = RawCNVkit@method)

  return(cnvkitf)

}



## FilterCore <- function(cnvcore, filter) {

## }
