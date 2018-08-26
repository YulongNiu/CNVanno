##' Internal annotation functions
##'
##' \itemize{
##'   \item \code{AnnoCNVOverlap_()}: Annotation of single CNV with one database with mutual overlap rate.
##'   \item \code{AnnoCNVType_()}: Annotation of single CNV with one database with gain and loss percentage. The input database should be the output of \code{AnnoCNVOverlap_()} function.
##'   \item \code{AnnoGainLossRateCheck_()}: Check the gain/loss annotation.
##' }
##' @title Annotation utilities
##' @param annodb A \code{tbl_df} indicating a single annotation database. In \code{AnnoCNVClinCore()}, it should at least contains "chromosome", "start", "end", "type" (gain or loss), and "clinical_significance" columns.
##' @inheritParams filterRow_
##' @inheritParams Merge
##' @return
##' \itemize{
##'   \item \code{AnnoCNVOverlap_()}: A \code{data.frame} with selected rows. The row number may be zero.
##'   \item \code{AnnoCNVType_()}: Same as the \code{AnnoCNVOverlap_()} function.
##'  \item \code{AnnoGainLossRateCheck_()}: A \code{logic} vector.
##' }
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @importFrom magrittr %>%
##' @importFrom dplyr filter
##' @rdname annointernal
##' @keywords internal
##'
AnnoCNVOverlap_ <- function(corerow,
                            annodb,
                            reciprate = 0.5) {

  ## step1: filter annodb
  annodb %<>% filter(chromosome == corerow$chromosome)
  if (nrow(annodb) == 0) {
    return(filter(annodb, FALSE))
  } else {}

  ## step2: select > reciprate anno
  anno <- OverlapRegionRate(corerow, annodb) %>%
    transmute(fRate > reciprate & tRate > reciprate) %>%
    unlist %>%
    filter(annodb, .)

  return(anno)
}


##' @param typerate A \code{numeric} value indicating the gain/loss rate.
##' @inheritParams AnnoCNVOverlap_
##' @importFrom magrittr %>%
##' @importFrom dplyr filter
##' @rdname annointernal
##' @keywords internal
##'
AnnoCNVType_ <- function(corerow,
                         annodb,
                         typerate = 0.7) {

  ## step1: check `annodb`. No need to check chromosome as AnnoCNVType_ is after AnnoCNVOverlap_
  if (nrow(annodb) == 0) {
    return(filter(annodb, FALSE))
  } else {}

  ## step2: select typerate
  anno <- annodb %>%
    filter(type == corerow$type)

  if (nrow(anno) / nrow(annodb) > typerate) {
    return(anno)
  } else {
    return(filter(annodb, FALSE))
  }
}
