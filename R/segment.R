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
              SortRegionChr %>% ## sort
              segMerge_(gap = gap) %>% ## segment
              mutate(method = raw@method) %>%
              SortRegionChr ## sort again

            res <- new('CoreCNV', coreCNV = cnvSeg)

            return(res)
          })


##' Segmentation internal functions.
##'
##' \itemize{
##'   \item \code{segMergeType_()}: Merege CNVs in the same chromosome in the same type.
##'   \item \code{segMergeChr_()}: Merege CNVs in the same chromosome.
##'   \item \code{segMerge_()}: Merge CNVs.
##' }
##' @title Internal functions for segmentation
##' @param cnv A \code{tbl_df} from the \code{rawCNV} slot of a \code{RawCNV} object.
##' @param chr A \code{string} indicating the chromosome, like "chr1", "chr2".
##' @param type A code{string} "gain", "loss", or "normal"
##' @inheritParams Segment
##' @return A \code{tbl_df}
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @importFrom dplyr filter arrange select everything bind_cols
##' @importFrom magrittr %<>% %>%
##' @rdname segutility
##' @keywords internal
##'
segMergeType_ <- function(cnv, gap, chr, type) {

  cnv %<>%
    select(start, end) %>%
    ReduceRegion(gap = gap) %>%
    mutate(chromosome = chr, type = type) %>%
    select(chromosome, everything())

  return(cnv)
}


##' @inheritParams segMergeType_
##' @inheritParams Segment
##' @importFrom magrittr %>% %<>%
##' @importFrom dplyr bind_rows arrange
##' @rdname segutility
##' @keywords internal
##'
segMergeChr_ <- function(cnv, gap, chr) {

  cnvList <- split(cnv, cnv$type)
  types <- names(cnvList)

  for (i in seq_along(cnvList)) {
    cnvList[[i]] %<>% segMergeType_(gap = gap,
                                    chr = chr,
                                    type = types[i])
  }

  cnvSeg <- bind_rows(cnvList)

  return(cnvSeg)
}

##' @inheritParams segMergeType_
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
      filter(chromosome == chrs[i]) %>%
      segMergeChr_(gap = gap,
                   chr = chrs[i])
  }

  cnvSeg <- bind_rows(cnvList) %>%
    mutate(chrnum = chromosome %>%
             str_extract('\\d+') %>%
             as.numeric) %>%
    arrange(chrnum) %>%
    select(-chrnum)

  return(cnvSeg)
}

