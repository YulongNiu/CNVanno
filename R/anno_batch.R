##' Batch annotation of CNVs
##'
##' \code{AnnoCNVBatch()}: Annotation of multiple CNVs to ClinVar/CliGen, population, and gene databases.
##'
##' @title The CNV batch annotation
##' @inheritParams write.cns
##' @inheritParams FindCyto
##' @param FUN \code{AnnoCNVClinCore()}, \code{AnnoCNVPopuCore()}, and \code{ANnoCNVGeneCore()}.
##' @param ... Parameters from the \code{FUN} functions
##' @importFrom doParallel registerDoParallel stopImplicitCluster
##' @importFrom foreach foreach %dopar%
##' @importFrom iterators iter
##' @importFrom magrittr %<>%
##' @rdname annobatch
##' @examples
##' require('magrittr')
##' data(CNVdb)
##'
##' cnsFile <- system.file('extdata', 'example.filter.cns', package = 'CNVanno')
##' cns <- cnsFile %>% read.cns %>% FilterCNS
##'
##' ## ClinGenGRCh38
##' AnnoCNVBatch(cns, AnnoCNVClinCore, CNVdb$ClinGenGRCh38)
##'
##' ## ClinVarGRCh38
##' AnnoCNVBatch(cns, AnnoCNVClinCore, CNVdb$ClinVarGRCh38,
##'              typeColName = 'Type',
##'              sigColName = 'ClinicalSignificance')
##' ## DGV
##' AnnoCNVBatch(cns, AnnoCNVPopuCore, CNVdb$DGV,
##'              'gain_frequence', 'loss_frequence',
##'              'variantsubtype')
##'
##' ## DECIPHER_POPULATION
##' AnnoCNVBatch(cns, AnnoCNVPopuCore, CNVdb$DECIPHER_POPULATION)
##'
##' ## ExAC_POPULATION
##' AnnoCNVBatch(cns, AnnoCNVPopuCore, CNVdb$ExAC_POPULATION, NA, NA, 'type')
##'
##' ## OMIM
##' AnnoCNVBatch(cns, AnnoCNVGeneCore, CNVdb$OMIMGRCh38)
##'
##' ## ClinGen_TriHaplo
##' AnnoCNVBatch(cns, AnnoCNVGeneCore, CNVdb$ClinGen_TriHaploGRCh38)
##'
##' ## DECIPHER_Haplo
##' AnnoCNVBatch(cns, AnnoCNVGeneCore, CNVdb$DECIPHER_Haplo)
##'
##' ## ExAC_pLI
##' AnnoCNVBatch(cns, AnnoCNVGeneCore, CNVdb$ExAC_pLI)
##'
##' @rdname clin
##' @export
##'
AnnoCNVBatch <- function(cns, FUN, ..., n = 2) {

  registerDoParallel(cores = n)

  itx <- iter(cns, by = 'row')
  annoList <- foreach(i = itx) %dopar% {
    return(FUN(cnsSingle = i, ...))
  }

  itx <- iter(cns, by = 'row')
  annoNames <- foreach(i = itx, .combine = c) %do% {
    return(paste0(c(i[1], ':', i[2], '-', i[3]), collapse = ''))
  }

  names(annoList) <- annoNames
  annoList %<>% MergeAnno

  ## stop multiple cores
  stopImplicitCluster()

  return(annoList)
}


##' @param annoList A \code{list} containing ClinGen annotation for each called CNV. The 1st element is a \code{data.frame} or \code{NULL}. The 2nd element and the rest are \code{character} vectors.
##' @importFrom magrittr %<>% %>%
##' @rdname annobatch
##' @return A merged \code{list}.
##' @keywords internal
##'
MergeAnno <- function(annoList) {

  annoLen <- length(annoList[[1]])
  res <- vector('list', annoLen)

  for (i in 2:annoLen) {

    res[[i]] <- annoList %>%
      lapply('[', i) %>%
      unlist %>%
      unname
  }

  annoLogic <- annoList %>%
    sapply(function(x) {return(!is.null(x[[1]]))})

  if (sum(annoLogic) > 0) {
    annoList %<>%
      `[`(annoLogic) %>%
      lapply(`[[`, 1)

    cnvNames <- names(annoList)
    annoMat <- do.call(rbind.data.frame,
                       c(annoList,
                       stringsAsFactors = FALSE))

    ## add rownames
    repNum <- sapply(annoList, nrow)
    annoMat %<>% cbind.data.frame(CNV = rep(cnvNames, repNum), ., stringsAsFactors = FALSE)
    rownames(annoMat) <- NULL

    res[[1]] <- annoMat

  } else {}

  return(res)
}
