##' Annotation core function for the DECIPHER_POPULATION, DGV, and ExAC_POPULATION databases
##'
##' \itemize{
##'   \item \code{AnnoCNVPopuCore()}: Annotation of single CNV to one population database.
##' }
##' @title The CNV population database annotation
##' @inheritParams AnnoCNVOverlap_
##' @inheritParams AnnoCNVType_
##' @return
##' \itemize{
##'   \item \code{AnnoCNVPopuCore()}: A \code{list}
##' }
##' @examples
##' data(CNVdb)
##' data(kit)
##'
##' ## DECIPHER_POPULATION
##' AnnoCNVPopuCore(kit@coreCNV[4, ], CNVdb$DECIPHER_POPULATION)
##'
##' ## DGV
##' AnnoCNVPopuCore(kit@coreCNV[4, ], CNVdb$DGVGRCh37)
##'
##' ## ExAC_POPULATION
##' AnnoCNVPopuCore(kit@coreCNV[3, ], CNVdb$ExAC_POPULATION)
##'
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @importFrom magrittr %>%
##' @importFrom dplyr mutate select everything
##' @rdname popu
##' @export
##'
AnnoCNVPopuCore <- function(corerow,
                            annodb,
                            reciprate = 0.5,
                            typerate = 0.7) {
  res <- vector('list', 2)
  res[2] <- NA
  names(res) <- c('PopuAnno', 'Percentage')

  ## step1: if cns is mapped to the population database.
  ## step2: gain and loss of the 'type' column

  anno <-  annodb %>%
    AnnoCNVOverlap_(corerow, ., reciprate)  %>%
    AnnoCNVType_(corerow, ., typerate) %>%
    mutate(CNV = FormatCorerow_(corerow)) %>%
    select(CNV, everything())
  res[[1]] <- anno

  ## step3: gain/loss percentage
  hasfreq <- c('gain_frequency', 'loss_frequency') %in% colnames(annodb)
  if (sum(hasfreq) == 2 & nrow(anno) > 0) {
    if (corerow$type == 'gain') {
      res[[2]] <- AnnoPopuPercentage_(anno$gain_frequency)
    } else {
      res[[2]] <- AnnoPopuPercentage_(anno$loss_frequency)
    }
  } else {}

  return(res)
}


##' @param annoPopuPer A \code{numeric} vector containing the percentage of CNV in the population.
##' @importFrom magrittr %<>%
##' @rdname popu
##' @return A \code{character} vector summary the CNV percentage.
##' @keywords internal
##'
AnnoPopuPercentage_ <- function(annoPopuPer) {

  popuper <- numeric(3)

  ## threshold <0.1%, 0.1%-1%, >1%
  popuper[1] <- sum(annoPopuPer < 0.001)
  popuper[2] <- sum((annoPopuPer >= 0.001) & (annoPopuPer < 0.01))
  popuper[3] <- length(annoPopuPer) - sum(popuper)

  popuper %<>%
    paste(c('<0.1%', '0.1%-1%', '>1%'), ., sep = ': ') %>%
    paste(collapse = '|')

  return(popuper)
}
