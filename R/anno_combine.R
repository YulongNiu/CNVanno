x##' @include AllClasses.R AllGenerics.R
NULL


##' Combine parameters from the raw CNV outputs
##'
##' Merge the parameters from the raw CNV outpus.
##' @title Combine parameters from \code{RawCNV}
##' @inheritParams CombinePara
##' @return A \code{tbl_df} object.
##' @examples
##' require('magrittr')
##' require('dplyr')
##' data(nator)
##' data(kit)
##'
##' mergecnv <- Merge(list(nator, kit), reciprate = 0.5, n = 2)
##' kitraw <- system.file('extdata', 'example.cnvkit', package = 'CNVanno') %>% read_cnvkit
##' kitraw@params %<>% select(log2, cn)
##'
##' tmp1 <- compressPara_(kitraw)
##'
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @importFrom magrittr %>%
##' @importFrom dplyr do group_by ungroup select everything
##' @importFrom tibble tibble
##' @importFrom methods new
##' @rdname CombinePara-methods
##' @exportMethod CombinePara
##'
SetMethod(f = 'CombinePara',
          signature = c(core = 'CoreCNV', raw = 'RawCNV', reciprate = 'numeric'),
          definition = function(core, raw, reciprate, ...) {

          })


##' CombinePara internal functions
##'
##' \itemize{
##'   \item \code{compressPara_()}: Compress the parameters of \code{RawCNV} into one column and combine it with the `rawCNV` slot.
##'   \item \code{combineRow_()}: Combine each row of the \code{RawCNV}.
##' }
##' @title Internal functions for CombinePara
##' @param raw A \code{RawCNV}.
##' @return
##' \itemize{
##'   \item \code{compressPara_()}: A \code{tbl_df}.
##'   \item \code{combineRow_()}: A \code{string}.
##' }
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @importFrom magrittr %>%
##' @importFrom dplyr bind_cols expr funs mutate_if transmute_all
##' @importFrom tidyr unite
##' @rdname combineparautility
##' @keywords internal
##'
compressPara_ <- function(raw) {

  res <- raw@params %>%
    mutate_if(is.numeric, round, 2) %>% ## round numeric columns
    transmute_all(.funs = funs(paste(expr(.), ., sep = '='))) %>% ## paste column names
    unite(params, 1 : ncol(.), sep = ',') %>% ## compress all columns
    bind_cols(raw@rawCNV, .)

  return(res)
}


##' @param rawpara A \code{tbl_df}
##' @inheritParams Merge
##' @inheritParams filterRow_
##' @importFrom dplyr transmute filter
##' @importFrom magrittr %<>% %>%
##' @rdname combineparautility
##' @keywords internal
##'
combineRow_ <- function(corerow, rawpara, reciprate) {

  rawpara %<>% filter(chromosome == corerow$chromosome)
  if (nrow(rawpara) == 0) {
    return('')
  } else {}

  res <- corerow %>%
    OverlapRegionRate(rawpara) %>%
    transmute(l = fRate > reciprate & tRate > reciprate) %>% ## select by reciprate
    unlist %>%
    filter(rawpara, .) %>%
    .$params %>% ## extract params
    paste(collapse = ';') %>%
    ifelse(nchar(.) > 0, ., '') ## empty params to ''

  return(res)
}
