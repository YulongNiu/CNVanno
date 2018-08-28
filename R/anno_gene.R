##' Annotation core function for the OMIM, ClinGen_TriHaplo, ExAC_plI, and DECIPHER_Haplo
##'
##' Annotation of single CNV to one gene database
##' @title The CNV gene database annotation
##' @param annoGene A \code{data.frame} of the OMIM, ClinGen_HaploTriplo, ExAC_plI, and DECIPHER_Haplo databases.
##' @inheritParams AnnoCNVClinCore
##' @return A \code{list}
##' @examples
##' data(CNVdb)
##' data(kit)
##'
##' ## ClinGen_TriHaploGRCh37
##' AnnoCNVGeneCore(kit@coreCNV[3, ], CNVdb$ClinGen_TriHaploGRCh37)
##'
##' ## ExAC_pLI
##' AnnoCNVGeneCore(kit@coreCNV[3, ], CNVdb$ExAC_pLI)
##'
##' ## DECIPHER_Haplo
##' AnnoCNVGeneCore(kit@coreCNV[3, ], CNVdb$DECIPHER_Haplo)
##'
##' ##OMIM
##' AnnoCNVGeneCore(kit@coreCNV[1, ], CNVdb$OMIMGRCh38)
##'
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @rdname cnvgene
##' @export
##'
AnnoCNVGeneCore <- function(corerow,
                            annodb) {
  res <- vector('list', 2)
  names(res) <- c('GeneAnno', 'LeaveBlank')

  res[[1]] <- AnnoGeneOverlap_(corerow, annodb)
  res[[2]] <- NA

  return(res)
}


##' @inheritParams AnnoCNVOverlap_
##' @inheritParams AnnoCNVType_
##' @importFrom stringr str_detect
##' @importFrom magrittr %>%
##' @importFrom dplyr bind_cols transmute mutate filter everything
##' @rdname cnvgene
##' @return A \code{tbl_df} summary annotated genes/regions.
##' @keywords internal
##'
AnnoGeneOverlap_ <- function(corerow,
                             annodb) {

  ## step1: if CNV is mapped to the gene databases.
  annodb %<>% filter(chromosome == corerow$chromosome)
  if (nrow(annodb) == 0) {
    return(filter(annodb, FALSE))
  } else {}

  ## step2: evaluate the relationship of gene and called CNV.
  anno <- corerow %>%
    OverlapRegionRate(annodb) %>%
    transmute(overlap_relation = case_when((fRate < 1 & tRate == 1) ~ 'GeneinCNV',
    (fRate == 1 & tRate < 1) ~ 'CNVinGene',
    (fRate == 1 & tRate == 1) ~ 'CNVeqGene',
    (fRate == 0 & tRate == 0) ~ 'nooverlap',
    TRUE ~ 'Overlap')) %>%
    bind_cols(annodb) %>%
    mutate(CNV = FormatCorerow_(corerow)) %>%
    select(CNV, chromosome:end, overlap_relation, everything()) %>%
    filter(overlap_relation != 'nooverlap')

  return(anno)
}

