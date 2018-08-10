##' Internal annotation functions
##'
##' \itemize{
##'   \item \code{AnnoCNVOverlap_()}: Annotation of single CNV with one database with mutual overlap rate.
##'   \item \code{AnnoCNVType_()}: Annotation of single CNV with one database with gain and loss percentage. The input database should be the output of \code{AnnoCNVOverlap_()} function.
##'   \item \code{AnnoGainLossRateCheck_()}: Check the gain/loss annotation.
##' }
##' @title Annotation utilities
##' @param cnsSingle A \code{character} vector indicating a single CNV.
##' @param annoSingle A \code{data.frame} indicating a single database.
##' @param mutualRate A \code{numeric} value indicating the mutual overlap rate.
##' @return
##' \itemize{
##'   \item \code{AnnoCNVOverlap_()}: A \code{data.frame} with selected rows. The row number may be zero.
##'   \item \code{AnnoCNVType_()}: Same as the \code{AnnoCNVOverlap_()} function.
##'  \item \code{AnnoGainLossRateCheck_()}: A \code{logic} vector.
##' }
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @importFrom magrittr %>% is_greater_than equals
##' @rdname annointernal
##' @keywords internal
##'
AnnoCNVOverlap_ <- function(cnsSingle,
                            annoSingle,
                            mutualRate = 0.5) {

  regionf <- as.numeric(cnsSingle[2:3])
  anno <- annoSingle[annoSingle[, 1] %in% cnsSingle[1], ]

  anno <- regionf %>%
    OverlapRegion(anno[, 2:3, drop = FALSE]) %>%
    anno[., , drop = FALSE]

  if (nrow(anno) > 0) {
    ## step 2: overlap rate
    anno <- regionf %>%
      OverlapRegionRate(anno[, 2:3, drop = FALSE]) %>%
      is_greater_than(mutualRate) %>%
      rowSums %>%
      equals(2) %>%
      anno[., , drop = FALSE]
  } else {}

  return(anno)
}


##' @param typeColName A \code{character} string indicating the gain/loss column name.
##' @param typeRate A \code{numeric} value indicating the gain/loss rate.
##' @inheritParams AnnoCNVOverlap_
##' @rdname annointernal
##' @keywords internal
##'
AnnoCNVType_ <- function(cnsSingle,
                         annoSingle,
                         typeColName,
                         typeRate = 0.7) {

  if (nrow(annoSingle) > 0) {
    tLogic <- cnsSingle[5] %>%
      is_greater_than(0) %>%
      ifelse('gain', 'loss') %>%
      AnnoGainLossRateCheck_(annoSingle[, typeColName])

    if (sum(tLogic) / nrow(annoSingle) >= typeRate) {
      annoSingle <- annoSingle[tLogic, ,drop = FALSE]
    } else {
      annoSingle <- annoSingle[0, ]
    }

  } else {}

  return(annoSingle)
}

##' @param detectType A \code{character} string, either "gain" or "loss".
##' @param annoType A \code{character} vector containing "gain" and "loss".
##' @importFrom stringr str_detect
##' @rdname annointernal
##' @keywords internal
##'
AnnoGainLossRateCheck_ <- function(detectType, annoType) {

  return(str_detect(annoType, detectType))

}
