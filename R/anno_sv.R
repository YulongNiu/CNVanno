##' CNV dbVar, ClinGen, and ClinVar databases annotation
##'
##' \itemize{
##'   \item \code{AnnoSVCore()}: Annotation of single CNV to the dbVar,  ClinGen, and ClinVar databases.
##' }
##' @title The dbVar, ClinGen, and ClinVar databases annotation
##' @inheritParams AnnoCNVOverlap_
##' @inheritParams AnnoCNVType_
##' @return
##' \itemize{
##'    \item \code{AnnoSVCore()}: A \code{list}
##' }
##'
##' @examples
##' data(CNVdb)
##' data(kit)
##'
##' ## dbVar
##' AnnoSVCore(kit@coreCNV[3, ], CNVdb$dbVarGRCh37)
##'
##' ## ClinGen
##' AnnoSVCore(kit@coreCNV[3, ], CNVdb$ClinGenGRCh37)
##'
##' ## ClinVar
##' AnnoSVCore(kit@coreCNV[3, ], CNVdb$ClinVarGRCh37)
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @importFrom magrittr %>%
##' @importFrom dplyr mutate select everything
##' @rdname clin
##' @export
##'
AnnoSVCore <- function(corerow,
                       annodb,
                       reciprate = 0.5,
                       typerate = 0.7) {

  res <- vector('list', 3)
  res[2:3] <- NA
  names(res) <- c('SVanno', 'conflict', 'summary')

  ## step1: if the cns is mapped ClinGen/ClinVar database
  ## step2: gain and loss of the 'TYPE' column
  anno <-  annodb %>%
    AnnoCNVOverlap_(corerow, ., reciprate) %>%
    AnnoCNVType_(corerow, ., typerate) %>%
    mutate(CNV = FormatCorerow_(corerow)) %>%
    select(CNV, everything())
  res[[1]] <- anno

  if (nrow(anno) > 0) {
    ## step3: benign/pathogenic
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
    sum(na.rm = TRUE) ## may contained NA
  pNum <- annoSig %>%
    str_detect('[p|P]athogenic') %>%
    sum(na.rm = TRUE) ## may contained NA

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


