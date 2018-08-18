##' @include AllClasses.R AllGenerics.R
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
##' kit <- system.file('extdata', 'exampleseg.cnvkit', package = 'CNVanno') %>% read_cnvkit %>% filter_cnvkit
##' kitseg <- Segment(kit, gap = 10L)
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @importFrom magrittr %>%
##' @importFrom dplyr mutate
##' @rdname Segment-methods
##' @exportMethod Segment
##'
setMethod(f = 'Segment',
          signature = c(raw = 'RawCNV', gap = 'integer'),
          definition = function(raw, gap, ...) {

            cnvSeg <- raw@rawCNV %>%
              segMerge_(gap = gap) %>%
              mutate(method = raw@method)

            res <- new('CoreCNV', coreCNV = cnvSeg)

            return(res)
          })


##' Segmentation internal functions.
##'
##' \itemize{
##'   \item \code{segPrepare_()}: Separate chromosome and order the star and end.
##'   \item \code{segMergeType_()}: Merege CNVs in the same chromosome in the same type.
##'   \item \code{segMergeChr_()}: Merege CNVs in the same chromosome.
##'   \item \code{segMerge_()}: Merge CNVs.
##' }
##' @title Internal functions for segmentation
##' @param cnv A \code{tbl_df} from the \code{rawCNV} slot of a \code{RawCNV} object.
##' @param chr A \code{character string} indicating the chromosome, like "chr1", "chr2".
##' @return A \code{tbl_df}
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @importFrom dplyr filter arrange select everything bind_cols
##' @importFrom magrittr %<>% %>%
##' @rdname segutility
##' @keywords internal
##'
segPrepare_ <- function(cnv, chr) {

  ## step1: sort rows
  st <- cnv %>%
    filter(chromosome == chr) %>%
    select(start:end) %>%
    SortRegion

  ## step2: sort columns
  cnv %<>%
    select(-(start:end)) %>%
    bind_cols(st, .) %>%
    select(chromosome, everything()) %>%
    arrange(start)

  return(cnv)
}



##' @inheritParams segPrepare_
##' @inheritParams Segment
##' @param type A \code{string} either "gain" or "loss".
##' @importFrom dplyr bind_cols select
##' @importFrom magrittr %>%
##' @rdname segutility
##' @keywords internal
##'
segMergeType_ <- function(cnv, gap, chr, type) {

  region  <- cnv %>%
    select(start, end) %>%
    ReduceRegion(gap = gap)

  cnvSeg <- list(chromosome = chr,
                 start = region$start,
                 end = region$end,
                 type = type) %>%
    bind_cols

  return(cnvSeg)
}


##' @inheritParams segPrepare_
##' @inheritParams Segment
##' @importFrom magrittr %>% %<>%
##' @importFrom dplyr bind_rows arrange
##' @rdname segutility
##' @keywords internal
##'
segMergeChr_ <- function(cnv, chr, gap) {

  cnvList <- split(cnv, cnv$type)
  types <- names(cnvList)

  for (i in seq_along(cnvList)) {
    cnvList[[i]] %<>% segMergeType_(gap = gap,
                                    chr = chr,
                                    type = types[i])
  }

  cnvSeg <- bind_rows(cnvList) %>%
    arrange(start)

  return(cnvSeg)
}


##' @inheritParams Segment
##' @importFrom magrittr %>% %<>%
##' @importFrom dplyr bind_rows arrange
##' @importFrom stringr str_extract
##' @rdname segutility
##' @keywords internal
##'
segMerge_ <- function(cnv, gap) {

  cnvList <- split(cnv, cnv$chromosome)
  chrs <- names(cnvList)

  for (i in seq_along(cnvList)) {
    cnvList[[i]] %<>%
      segPrepare_(chrs[i]) %>%
      segMergeChr_(chr = chrs[i],
                   gap = gap)
  }

  cnvSeg <- bind_rows(cnvList) %>%
    mutate(chrnum = chromosome %>%
             str_extract('\\d+') %>%
             as.numeric) %>%
    arrange(chrnum) %>%
    select(-chrnum)

  return(cnvSeg)
}

