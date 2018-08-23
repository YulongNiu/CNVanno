##' @include AllClasses.R AllGenerics.R region.R
NULL


##' Segmentation of \code{RawCNV}.
##'
##' Merge the concatenated CNV at the same chromosome with the same type (gain and loss).
##'
##' @title Segmentation
##' @inheritParams Segment
##' @return A \code{CoreCNV} object.
##' @examples
##' library('magrittr')
##'
##' kit <- system.file('extdata', 'exampleseg.cnvkit', package = 'CNVanno') %>%
##'   read_cnvkit %>%
##'   filter_cnvkit
##' kitseg <- Segment(kit, gap = 10L)
##' kitsegmore <- Segment(kitseg, gap = 2000000L)
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @importFrom magrittr %>%
##' @importFrom methods new
##' @importFrom dplyr mutate
##' @rdname Segment-methods
##' @exportMethod Segment
##'
setMethod(f = 'Segment',
          signature = c(raw = 'RawCNV', gap = 'integer'),
          definition = function(raw, gap, ...) {

            res <- raw@rawCNV %>%
              mutate(method = raw@method) %>%
              segmentcore_(gap = gap) %>%
              new('CoreCNV', coreCNV = .)

            return(res)
          })



##' @inheritParams Segment
##' @importFrom magrittr %>%
##' @importFrom methods new
##' @rdname Segment-methods
##' @exportMethod Segment
##'
setMethod(f = 'Segment',
          signature = c(raw = 'CoreCNV', gap = 'integer'),
          definition = function(raw, gap, ...) {

            res <- raw@coreCNV %>%
              segmentcore_(gap = gap) %>%
              new('CoreCNV', coreCNV = .)

            return(res)
          })


##' Segment internal functions
##'
##' @title Internal functions for segment
##' @param corecnv A \code{tbl_df} contains "chromosome", "start", "end", "type", and "method" columns.
##' @inheritParams Segment
##' @return A \code{tbl_df}.
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @importFrom magrittr %>% %<>%
##' @importFrom dplyr group_by mutate do select everything ungroup
##' @importFrom tibble tibble
##' @keywords internal
##'
segmentcore_ <- function(corecnv, gap) {

  corecnv %<>%
    group_by(chromosome, type, method) %>% ## sort reduce
    do(SortRegion(tibble(start = .$start, end = .$end))) %>%
    do(ReduceRegion(tibble(start = .$start, end = .$end), gap = gap)) %>% ## no need to sort again as reduce is done in order
    select(chromosome, start:end, everything()) %>% ## columns in right order
    ungroup

  return(corecnv)
}
