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
##' kitseg <- Segment(kit, interlen = 10L)
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @importFrom magrittr %>%
##' @importFrom dplyr mutate
##' @rdname Segment-methods
##' @exportMethod Segment
##'
setMethod(f = 'Segment',
          signature = c(raw = 'RawCNV', interlen = 'integer'),
          definition = function(raw, interlen = 100L, ...) {

            cnvSeg <- raw@rawCNV %>%
              segMerge_(interlen = interlen) %>%
              mutate(method = raw@method)

            return(cnvSeg)
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
##' @importFrom dplyr filter arrange
##' @rdname segutility
##' @keywords internal
##'
segPrepare_ <- function(cnv, chr) {

  cnv %<>%
    filter(chromosome == chr) %>%
    arrange(start)

  return(cnv)
}



##' @inheritParams segPrepare_
##' @inheritParams Segment
##' @param type A \code{string} either "gain" or "loss".
##' @importFrom tibble tibble
##' @rdname segutility
##' @keywords internal
##'
segMergeType_ <- function(cnv, chr, interlen, type) {

  ## step1: start and end vector
  start <- cnv$start
  end <- cnv$end

  ## step2: select start end
  inter <- c(FALSE, start[-1] - end[-length(end)] < interlen)
  startIdx <- which(!inter)
  endIdx <- c(startIdx[-1] - 1,
              length(inter))

  ## step3: new cnv
  cnvSeg <- tibble(
    chromosome = chr,
    start = start[startIdx],
    end = end[endIdx],
    type = type)

  return(cnvSeg)
}


##' @inheritParams segPrepare_
##' @inheritParams Segment
##' @importFrom magrittr %>% %<>%
##' @importFrom dplyr bind_rows arrange
##' @rdname segutility
##' @keywords internal
##'
segMergeChr_ <- function(cnv, chr, interlen) {

  cnvList <- split(cnv, cnv$type)
  types <- names(cnvList)

  for (i in seq_along(cnvList)) {
    cnvList[[i]] %<>% segMergeType_(chr = chr,
                                    interlen = interlen,
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
segMerge_ <- function(cnv, interlen) {

  cnvList <- split(cnv, cnv$chromosome)
  chrs <- names(cnvList)

  for (i in seq_along(cnvList)) {
    cnvList[[i]] %<>%
      segPrepare_(chrs[i]) %>%
      segMergeChr_(chr = chrs[i],
                  interlen = interlen)
  }

  cnvSeg <- bind_rows(cnvList) %>%
    mutate(chrnum = chromosome %>%
             str_extract('\\d+') %>%
             as.numeric) %>%
    arrange(chrnum) %>%
    select(-chrnum)

  return(cnvSeg)
}

