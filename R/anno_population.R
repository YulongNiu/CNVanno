##' Annotation core function for the DECIPHER_POPULATION, DGV, and ExAC_POPULATION databases
##'
##' \itemize{
##'   \item \code{AnnoCNVPopuCore()}: Annotation of single CNV to one population database.
##' }
##' @title The CNV population database annotation
##' @inheritParams AnnoCNVClinCore
##' @param annoPopo A \code{data.frame} of the DECIPHER_POPULATION, DGV, and ExAC_POPULATION databases.
##' @param gainPerColName A \code{character} string indicating the gain column name.
##' @param lossPerColName A \code{character} string indicating the loss column name.
##' @return
##' \itemize{
##'   \item \code{AnnoCNVPopuCore()}: A \code{list}
##' }
##' @examples
##' require('magrittr')
##' data(CNVdb)
##'
##' cnsFile <- system.file('extdata', 'example.cnvkit', package = 'CNVanno')
##' cns <- cnsFile %>% read.cns %>% FilterCNS
##'
##' ## DECIPHER_POPULATION
##' AnnoCNVPopuCore(cns[1, ], CNVdb$DECIPHER_POPULATION)
##'
##' ## DGV
##' AnnoCNVPopuCore(cns[1, ], CNVdb$DGV, 'gain_frequence', 'loss_frequence', 'variantsubtype')
##'
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @rdname popu
##' @export
##'
AnnoCNVPopuCore <- function(cnsSingle,
                            annoPopu,
                            gainPerColName = 'duplication_frequency',
                            lossPerColName = 'deletion_frequency',
                            typeColName = 'cnvtype',
                            mutualRate = 0.5,
                            typeRate = 0.7) {
  res <- vector('list', 2)
  res[2] <- NA
  names(res) <- c('PopuAnno', 'Percentage')

  ## step1: if cns is mapped to the population database.
  ## step2: gain and loss of the 'TYPE' column

  anno <-  annoPopu %>%
    AnnoCNVOverlap_(cnsSingle, ., mutualRate)  %>%
    AnnoCNVType_(cnsSingle, ., typeColName, typeRate)
  res[[1]] <- anno

  ## step3: gain/loss percentage
  if (!(is.na(gainPerColName) |
        is.na(lossPerColName)) &
      nrow(anno) > 0) {
    res[[2]] <- ifelse(cnsSingle[5] > 0,
                       gainPerColName,
                       lossPerColName) %>%
      {AnnoPopuPercentage_(anno[, .])}

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
