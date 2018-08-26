##' CNV ClinGen and ClinVar database annotation
##'
##' \itemize{
##'   \item \code{AnnoCNVClinCore()}: Annotation of single CNV to the ClinGen and ClinVar database.
##' }
##' @title The ClinGen and ClinVar database annotation
##' @inheritParams AnnoCNVOverlap_
##' @inheritParams AnnoCNVType_
##' @return
##' \itemize{
##'    \item \code{AnnoCNVClinCore()}: A \code{list}
##' }
##'
##' @examples
##' data(kit)
##'
##' ## ClinGen
##' AnnoCNVClinCore(kit@coreCNV[3, ], CNVdb$ClinGenGRCh37)
##'
##' ## ClinVar
##' AnnoCNVClinCore(kit@coreCNV[3, ], CNVdb$ClinVarGRCh37)
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @importFrom magrittr %>%
##' @rdname clin
##' @export
##'
AnnoCNVClinCore <- function(corerow,
                            annodb,
                            reciprate = 0.5,
                            typerate = 0.7) {

  res <- vector('list', 3)
  res[2:3] <- NA
  names(res) <- c('ClinAnno', 'conflict', 'summary')

  ## step 1: if the cns is mapped ClinGen/ClinVar database
  ## step 2: gain and loss of the 'TYPE' column
  anno <-  annodb %>%
    AnnoCNVOverlap_(corerow, ., reciprate) %>%
    AnnoCNVType_(corerow, ., typerate)
  res[[1]] <- anno

  if (nrow(anno) > 0) {
    ## step 3: benign/pathogenic
    bpList <- AnnoBenignPathoCheck_(anno$clinical_significance)
    res[[2]] <- bpList[[1]]
    res[[3]] <- bpList[[2]]
  } else {}

  return(res)
}


##' @param annoSig A \code{character} vector containing "benign" and "pathogentic".
##' @importFrom stringr str_detect
##' @importFrom magrittr %>%
##' @rdname clin
##' @return A \code{list} summary the clinical significance.
##' @keywords internal
##'
AnnoBenignPathoCheck_ <- function(annoSig) {

  bNum <- annoSig %>%
    str_detect('[b|B]enign') %>%
    sum
  pNum <- annoSig %>%
    str_detect('[p|P]athogenic') %>%
    sum

  conflict <- ifelse(bNum == 0 | pNum ==0, 'No', 'Yes')

  sumSig <- annoSig %>%
    table %>% {
      paste(names(.), ., sep = ': ')
    } %>%
    paste(collapse = '|')


  res <- list(conflict = conflict,
              summary = sumSig)

  return(res)
}


