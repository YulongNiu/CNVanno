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
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @importFrom magrittr %>%
##' @importFrom dplyr group_by mutate do select everything ungroup
##' @importFrom methods new
##' @rdname Segment-methods
##' @exportMethod Segment
##'
setMethod(f = 'Segment',
          signature = c(raw = 'RawCNV', gap = 'integer'),
          definition = function(raw, gap, ...) {

            cnvSeg <- raw@rawCNV %>%
              group_by(chromosome, type) %>% ## sort reduce
              do(SortRegion(tibble(start = .$start, end = .$end))) %>%
              do(ReduceRegion(tibble(start = .$start, end = .$end), gap = gap)) %>% ## no need to sort again as reduce is done in order
              mutate(method = raw@method) %>% ## columns in right order
              select(chromosome, start:end, everything()) %>%
              ungroup

            res <- new('CoreCNV', coreCNV = cnvSeg)

            return(res)
          })
