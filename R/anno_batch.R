##' @include AllClasses.R AllGenerics.R region.R anno_utilities.R anno_population.R anno_gene.R anno_clin.R
NULL

##' Batch annotation of CNVs
##'
##' Annotation of multiple CNVs to ClinVar/ClinGen, population, and gene databases.
##'
##' @title The CNV batch annotation
##' @inheritParams AnnoCNVBatch
##' @importFrom doParallel registerDoParallel stopImplicitCluster
##' @importFrom foreach foreach %dopar%
##' @importFrom iterators iter
##' @importFrom magrittr %<>%
##' @examples
##' data(CNVdb)
##' data(kit)
##'
##' ## ClinGenGRCh37
##' AnnoCNVBatch(kit, AnnoCNVClinCore, CNVdb$ClinGenGRCh37,
##'              reciprate = 0.5, typerate = 0.7, n = 2)
##'
##' ## ClinVarGRCh37
##' AnnoCNVBatch(kit, AnnoCNVClinCore, CNVdb$ClinVarGRCh38,
##'              reciprate = 0.5, typerate = 0.7, n = 2)
##'
##' ## DGVGRCh37
##' AnnoCNVBatch(kit, AnnoCNVPopuCore, CNVdb$DGVGRCh37,
##'              reciprate = 0.5, typerate = 0.7, n = 2)
##'
##' ## DECIPHER_POPULATION
##' AnnoCNVBatch(kit, AnnoCNVPopuCore, CNVdb$DECIPHER_POPULATION,
##'              reciprate = 0.5, typerate = 0.7, n = 2)
##'
##' ## ExAC_POPULATION
##' AnnoCNVBatch(kit, AnnoCNVPopuCore, CNVdb$ExAC_POPULATION,
##'              reciprate = 0.5, typerate = 0.7, n = 2)
##'
##' ## OMIM
##' AnnoCNVBatch(kit, AnnoCNVGeneCore, CNVdb$OMIMGRCh38, n = 2)
##'
##' ## ClinGen_TriHaplo
##' AnnoCNVBatch(kit, AnnoCNVGeneCore, CNVdb$ClinGen_TriHaploGRCh37, n = 2)
##'
##' ## DECIPHER_Haplo
##' AnnoCNVBatch(kit, AnnoCNVGeneCore, CNVdb$DECIPHER_Haplo, n = 2)
##'
##' ## ExAC_pLI
##' AnnoCNVBatch(kit, AnnoCNVGeneCore, CNVdb$ExAC_pLI, n = 2)
##'
##' ## RefGeneGRCh37
##' AnnoCNVBatch(kit, AnnoCNVGeneCore, CNVdb$RefGeneGRCh37, n = 2)
##' @rdname AnnoCNVBatch-methods
##' @exportMethod AnnoCNVBatch
##'
setMethod(f = 'AnnoCNVBatch',
          signature = c(core = 'CoreCNV', FUN = 'function'),
          definition = function(core, FUN, ..., n) {

            registerDoParallel(cores = n)

            itx <- iter(core@coreCNV, by = 'row')
            annoList <- foreach(i = itx) %dopar% {
              return(FUN(i, ...))
            }

            annoList %<>% mergeAnno_

            ## stop multiple cores
            stopImplicitCluster()

            return(annoList)
          })


##' Internal merge annotation functions
##'
##' @title Merge CNV annotation
##' @param annoList A \code{list} containing CNV annotation for each called CNV. The 1st element is a \code{data.frame} or \code{NULL}. The 2nd element and the rest are \code{character} vectors.
##' @importFrom foreach foreach %do%
##' @importFrom dplyr bind_rows
##' @importFrom stringr str_extract
##' @importFrom magrittr %<>% %>%
##' @return A merged \code{list}.
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @keywords internal
##'
mergeAnno_ <- function(annoList) {

  annoLen <- length(annoList[[1]])
  res <- vector('list', annoLen)

  res[[1]] <- foreach(i = seq_along(annoList), .combine = bind_rows) %do% {
    return(annoList[[i]][[1]])
  }

  for (i in seq_len(annoLen)[-1]) {
    res[[i]] <- sapply(annoList, function(x){return(x[[i]])})
  }

  names(res) <- names(annoList[[1]])

  return(res)
}
