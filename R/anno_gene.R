##' Annotation core function for the OMIM, ClinGen_TriHaplo, ExAC_plI, and DECIPHER_Haplo
##'
##' Annotation of single CNV to one gene database
##' @title The CNV gene database annotation
##' @param annoGene A \code{data.frame} of the OMIM, ClinGen_HaploTriplo, ExAC_plI, and DECIPHER_Haplo databases.
##' @inheritParams AnnoCNVClinCore
##' @return A \code{list}
##' @examples
##' require('magrittr')
##' data(CNVdb)
##'
##' cnsFile <- system.file('extdata', 'example.cnvkit', package = 'CNVanno')
##' cns <- cnsFile %>% read.cns %>% FilterCNS
##'
##' ##OMIM
##' AnnoCNVGeneCore(cns[1, ], CNVdb$OMIMGRCh38)
##'
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @rdname cnvgene
##' @export
##'
AnnoCNVGeneCore <- function(cnsSingle,
                            annoGene) {
  res <- vector('list', 2)
  res[2] <- NA
  names(res) <- c('GeneAnno', 'LeaveBlank')

  ## step1: if cns is mapped to the gene databases.
  ## step2: evaluate the relationship of gene and called CNV.
  anno <- annoGene %>%
    AnnoGeneOverlap_(cnsSingle, .)

  if (nrow(anno) > 0) {
    res[[1]] <- anno
  } else {}

  return(res)
}


##' @inheritParams AnnoCNVClinCore
##' @importFrom stringr str_detect
##' @importFrom magrittr %>%
##' @rdname cnvgene
##' @return A \code{data.frame} summary annotated genes/regions.
##' @keywords internal
##'
AnnoGeneOverlap_ <- function(cnsSingle,
                             annoSingle) {

  regionf <- as.numeric(cnsSingle[2:3])
  anno <- annoSingle[annoSingle[, 1] %in% cnsSingle[1], ]

  anno <- regionf %>%
    OverlapRegion(anno[, 2:3, drop = FALSE]) %>%
    anno[., , drop = FALSE]

  if (nrow(anno) > 0) {
    locrela <- regionf %>%
      OverlapRegionRate(anno[, 2:3, drop = FALSE]) %>%
      OverlapRelation_
    anno <- cbind.data.frame(anno[, 1:3, drop = FALSE],
                             OverlapRelation = locrela,
                             anno[, -1:-3, drop = FALSE],
                             stringsAsFactors = FALSE)
  } else {}

  return(anno)
}


##' @param rateMat A \code{numeric} matrix with two columns. 1st column is the rate of overlapped region on overlap region, and 2nd column is the rate of annotated gene/region.
##' @importFrom magrittr %>%
##' @rdname cnvgene
##' @return A \code{character vector} summary coverage of genes/regions.
##' @keywords internal
##'
OverlapRelation_ <- function(rateMat) {

  res <- rateMat %>%
    nrow %>%
    rep('Overlap', .)

  gincLogic <- (rateMat[, 1] == 1) &
    (rateMat[, 2] < 1)
  res[gincLogic] <- 'CNVinGene'

  cingLogic <- (rateMat[, 1] < 1) &
    (rateMat[, 2] == 1)
  res[cingLogic] <- 'GeneinCNV'

  ceqgLogic <- (rateMat[, 1] == 1) &
    (rateMat[, 2] == 1)
  res[ceqgLogic] <- 'CNVeqGene'

  return(res)
}
