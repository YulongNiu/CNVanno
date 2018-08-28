##' The CNV region summary table according to Sun's format.
##'
##' \itemize{
##'   \item CNV
##'   \item size(100kb)
##'   \item type
##'   \item method
##'   \item cytoband
##'   \item sample
##' }
##'
##' @title Sun's summary CNV table
##' @inheritParams Cytoband
##' @param sampleType The sample type, like "proband", "mother", "father".
##' @param ... Parameters passed to the \code{Cytoband()} in this package.
##' @return A \code{tbl_df} contains Sun's CNV region format.
##' @examples
##' require('magrittr')
##' data(hg19cyto)
##' data(kit)
##'
##' SunCNVregion(kit, hg19cyto, n = 2)
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @importFrom magrittr %>% %<>%
##' @importFrom stringr str_extract
##' @importFrom dplyr mutate select slice
##' @export
##'
SunCNVregion <- function(core, cyto, sampleType = 'proband', ...) {

  cnv  <- core@coreCNV %>%
    mutate(CNV = paste(chromosome, paste(start, end, sep = '-'), sep = ':')) %>% ## CNV
    mutate(`size(100kb)` = round((end - start) / 100000, 1)) %>% ## size
    mutate(cytoband = Cytoband(core, cyto, ...)) %>% ## cytoband
    mutate(sample = sampleType) %>%
    select(CNV, `size(100kb)`, type, method, cytoband, sample)

  ## order
  ocnv <- core@coreCNV$chromosome %>%
    str_extract('\\d+') %>%
    as.numeric %>%
    order
  cnv %<>% slice(ocnv)

  return(cnv)
}



##' The CNV gene summary table according to Sun's format.
##'
##' ## column names for OMIM
##' \itemize{
##'   \item CNV
##'   \item overlap_relation
##'   \item gene_symbol
##'   \item Phenotypes
##'   \item Phenotypes CN
##'   \item MIMNumber
##'   \item MIMType
##' }
##' ## column names for ClinGenTriHaplo
##' \itemize{
##'   \item Haploinsufficiency Score
##'   \item Triplosensitivity Score
##' }
##'
##' ## column names for DECIPHERHaplo
##' \itemize{
##'   \item haploinsufficiency
##' }
##'
##' ## column names for pLI
##' \itemize{
##'   \item pLI
##' }
##'
## column names for DDG2P
##' \itemize{
##'   \item disease.name
##'   \item DDD.category
##'   \item allelic.requirement
##'   \item mutation.consequence
##'   \item phenotypes
##' }
##'
##' @title Merge genes database
##' @param gdbList A list containing CNV genes database. 1st element is the OMIM database, 2nd is the ClinGenTriHaplo, 3rd is the DECIPHERHaplo, 4th is the ExACpLI, 5th is the DECPHER_DDG2P
##' @return A \code{data.frame}
##' @examples
##' require('magrittr')
##' data(CNVdb)
##' data(kit)
##'
##' refGene <- AnnoCNVBatch(kit, AnnoCNVGeneCore, CNVdb$RefGeneGRCh37, n = 2)
##' gdbList <- list()
##' gdbList[[1]] <- AnnoCNVGeneRefGene2OMIM(refGene[[1]], CNVdb$OMIMGRCh38)
##' gdbList[[2]] <- AnnoCNVBatch(kit, AnnoCNVGeneCore, CNVdb$ClinGen_TriHaploGRCh37, n = 2)[[1]]
##' gdbList[[3]] <- AnnoCNVBatch(kit, AnnoCNVGeneCore, CNVdb$DECIPHER_Haplo, n = 2)[[1]]
##' gdbList[[4]] <- AnnoCNVBatch(kit, AnnoCNVGeneCore, CNVdb$ExAC_pLI, n = 2)[[1]]
##' gdbList[[5]] <- AnnoCNVGeneRefGene2DDG2P(refGene[[1]], CNVdb$DECIPHER_DDG2P)
##'
##' geneTable <- SunCNVgene(gdbList)
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @importFrom magrittr %<>% %>%
##' @importFrom dplyr rename select group_by distinct ungroup mutate slice
##' @importFrom stringr str_extract
##' @export
##'
SunCNVgene <- function(gdbList) {

  ## step1: process dbs
  gdbList[[1]] %<>%
    rename(gene_symbol = `Approved Symbol`) %>%
    select(CNV:overlap_relation, gene_symbol, Phenotypes, PhenotypesCN, MIMNumber, MIMType) %>%
    group_by(CNV) %>%
    distinct(gene_symbol, .keep_all = TRUE) %>%
    ungroup


  gdbList[[2]] %<>%
    rename(gene_symbol = `Gene Symbol`) %>%
    select(CNV:overlap_relation, gene_symbol, `Haploinsufficiency Score`, `Triplosensitivity Score`) %>%
    group_by(CNV) %>%
    distinct(gene_symbol, .keep_all = TRUE) %>%
    ungroup

  gdbList[[3]] %<>%
    mutate(gene_symbol = name %>% strsplit(split = '|', fixed = TRUE) %>% sapply(`[`, 1) %>% unlist) %>%
    select(CNV:overlap_relation, gene_symbol, haploinsufficiency) %>%
    group_by(CNV) %>%
    distinct(gene_symbol, .keep_all = TRUE) %>%
    ungroup

  gdbList[[4]] %<>%
    rename(gene_symbol = gene) %>%
    select(CNV:overlap_relation, gene_symbol, pLI) %>%
    group_by(CNV) %>%
    distinct(gene_symbol, .keep_all = TRUE) %>%
    ungroup

  gdbList[[5]] %<>%
    rename(gene_symbol = gene.symbol) %>%
    select(CNV:overlap_relation, gene_symbol, disease.name:phenotypes) %>%
    select(-disease.mim) %>%
    group_by(CNV) %>%
    distinct(gene_symbol, .keep_all = TRUE) %>%
    ungroup

  ## step2: merge
  gdb <- Reduce(mergeTwoGenedb_, gdbList)

  ## step3: order
  ocnv <- gdb$chromosome %>%
    str_extract('\\d+') %>%
    as.numeric %>%
    order
  gdb %<>% slice(ocnv)

  return(gdb)
}


##' Internal Sun's template annotation functions.
##'
##' \itemize{
##'   \item \code{mergeTwoGenedb_()}: Merge two gene results.
##' }
##'
##' @title Internal Sun annotation functions
##' @param gdb1, gdb2: Annotated gene databases.
##' @return \itemize{
##'   \item \code{mergeTwoGenedb_()}: A \code{tbl_df}.
##' }
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @importFrom magrittr %>% %<>%
##' @importFrom dplyr rename mutate if_else select full_join
##' @keywords internal
##'
mergeTwoGenedb_ <- function(gdb1, gdb2) {

  ## join
  gdb <- full_join(gdb1, gdb2, by = c('CNV', 'gene_symbol'))

  ## remove NA
  gdb %<>%
    rename(chromosome = chromosome.x, start = start.x, end = end.x, overlap_relation = overlap_relation.x) %>%
    mutate(chromosome = chromosome %>% {if_else(is.na(.), chromosome.y, .)}) %>%
    mutate(start = start %>% {if_else(is.na(.), start.y, .)}) %>%
    mutate(end = end %>% {if_else(is.na(.), end.y, .)}) %>%
    mutate(overlap_relation = overlap_relation %>% {if_else(is.na(.), overlap_relation.y, .)}) %>%
    select(-(chromosome.y:overlap_relation.y))

  gdb[] <- lapply(gdb, function(x){return(ifelse(is.na(x), '', x))})

  return(gdb)
}

##' Combine CNVregion and CNVgene tables according to Sun's format.
##'
##' \itemize{
##'   \item \code{CrossRegionGeneTable()}: Merge region and gene tables.
##' }
##'
##' @title Combine CNVregion and CNVgene tables.
##' @param regionTable A \code{tbl_df} of processed region table, which often comes from the \code{SunCNVregion()} in this package.
##' @param geneTable A \code{tbl_df} of gene table, which often comes from the \code{SunCNVgene()} in this package.
##' @return A list containing two \code{tbl_df} represents the processed \code{regionTable} and \code{geneTable} respectively.
##' @examples
##' require('magrittr')
##' data(CNVdb)
##' data(kit)
##' data(hg19cyto)
##'
##' refGene <- AnnoCNVBatch(kit, AnnoCNVGeneCore, CNVdb$RefGeneGRCh37, n = 2)
##' gdbList <- list()
##' gdbList[[1]] <- AnnoCNVGeneRefGene2OMIM(refGene[[1]], CNVdb$OMIMGRCh38)
##' gdbList[[2]] <- AnnoCNVBatch(kit, AnnoCNVGeneCore, CNVdb$ClinGen_TriHaploGRCh37, n = 2)[[1]]
##' gdbList[[3]] <- AnnoCNVBatch(kit, AnnoCNVGeneCore, CNVdb$DECIPHER_Haplo, n = 2)[[1]]
##' gdbList[[4]] <- AnnoCNVBatch(kit, AnnoCNVGeneCore, CNVdb$ExAC_pLI, n = 2)[[1]]
##' gdbList[[5]] <- AnnoCNVGeneRefGene2DDG2P(refGene[[1]], CNVdb$DECIPHER_DDG2P)
##'
##' geneTable <- SunCNVgene(gdbList)
##' regionTable <- SunCNVregion(kit, hg19cyto, n = 2)
##' grList <- CrossRegionGeneTable(regionTable, geneTable)
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @importFrom magrittr %>%
##' @importFrom stringr str_detect
##' @importFrom dplyr filter rename select group_by summarise_all funs left_join
##' @export
##'
CrossRegionGeneTable <- function(regionTable, geneTable) {

  ## step1: separate `geneTable`
  rT <- geneTable %>%
    filter(!str_detect(gene_symbol, 'ISCA-\\d+'))

  rTisca <- geneTable %>%
    filter(str_detect(gene_symbol, 'ISCA-\\d+')) %>%
    rename(ISCA_ID = gene_symbol) %>%
    select(-(chromosome:overlap_relation)) %>%
    group_by(CNV) %>%
    summarise_all(funs(paste(unique(.), collapse = '|')))

  ## step2: merge ICSA to `regionTable`
  gT <- regionTable %>%
    left_join(rTisca, by = 'CNV')
  gT[] <- lapply(gT, function(x){return(ifelse(is.na(x), '', x))})

  res <- list(regionTable = rT,
              geneTable = gT)

  return(res)
}


##' Annotation CNV with the OMIM (GRCh38) database from the RefGene (GRCh37)
##'
##' Step1: Normal annotation CNV with the RefGene (GRCh37) database.
##'
##' Step2: Use Entrez Gene ID as keys to combine RefGene and OMIM database.
##'
##' @title Annotation OMIM (GRCh38) from the RefGene (GRCh37)
##' @param refgeneTable A \code{data.frame} of annotated RefGene table.
##' @param OMIMdb The OMIM database.
##' @return A \code{list} of annotated OMIM.
##' @importFrom magrittr %>%
##' @importFrom dplyr select group_by distinct ungroup inner_join
##' @examples
##' data(CNVdb)
##' data(kit)
##'
##' refGene <- AnnoCNVBatch(kit, AnnoCNVGeneCore, CNVdb$RefGeneGRCh37, n = 2)
##'
##' ## cross OMIM database
##' crossOMIM <- AnnoCNVGeneRefGene2OMIM(refGene[[1]], CNVdb$OMIMGRCh38)
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @export
##'
AnnoCNVGeneRefGene2OMIM <- function(refgeneTable, OMIMdb) {

  ## step1: remove redundant geneid
  refTable <- refgeneTable %>%
    select(CNV:Entrez) %>%
    group_by(CNV) %>%
    distinct(Entrez, .keep_all = TRUE) %>% ## only keep the first same Entrez column
    ungroup

  ## step2: combine with OMIM db and replce OMIM db with GRCh37 location info
  res <- OMIMdb %>%
    select(-(chromosome:end)) %>%
    inner_join(refTable, ., by = 'Entrez')

  return(res)
}


##' Annotation CNV with the DECIPHER_DDG2Pdb (GRCh37) database from the RefGene (GRCh37)
##'
##' Step1: Normal annotation CNV with the RefGene (GRCh37) database.
##'
##' Step2: Use Entrez Gene ID as keys to combine RefGene and DECIPHER_DDG2Pdb database.
##'
##' @title Annotation DECIPHER_DDG2Pdb (GRCh37) from the RefGene (GRCh37)
##' @param DDG2Pdb The DECIPHER_DDG2Pdb database.
##' @inheritParams AnnoCNVGeneRefGene2OMIM
##' @return A \code{list} of annotated DECIPHER_DDG2Pdb
##' @examples
##' require('magrittr')
##' data(CNVdb)
##' data(kit)
##'
##' refGene <- AnnoCNVBatch(kit, AnnoCNVGeneCore, CNVdb$RefGeneGRCh37, n = 2)
##' crossDDG2P <- AnnoCNVGeneRefGene2DDG2P(refGene[[1]], CNVdb$DECIPHER_DDG2P)
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @export
##'
AnnoCNVGeneRefGene2DDG2P <- function(refgeneTable, DDG2Pdb) {

  ## step1: connect RefGene and OMIMdb with Entrez ID
  refTable <- refgeneTable %>%
    select(CNV:Entrez) %>%
    group_by(CNV) %>%
    distinct(Entrez, .keep_all = TRUE) %>% ## only keep the first same Entrez column
    ungroup

  ## step2: combine with OMIM db and replce OMIM db with GRCh37 location info
  res <- DDG2Pdb %>%
        inner_join(refTable, ., by = 'Entrez')

  return(res)
}





## refGene <- AnnoCNVBatch(kit, AnnoCNVGeneCore, CNVdb$RefGeneGRCh37, n = 2)
## gdbList <- list()
## gdbList[[1]] <- AnnoCNVGeneRefGene2OMIM(refGene[[1]], CNVdb$OMIMGRCh38)
## gdbList[[2]] <- AnnoCNVBatch(kit, AnnoCNVGeneCore, CNVdb$ClinGen_TriHaploGRCh37, n = 2)[[1]]
## gdbList[[3]] <- AnnoCNVBatch(kit, AnnoCNVGeneCore, CNVdb$DECIPHER_Haplo, n = 2)[[1]]
## gdbList[[4]] <- AnnoCNVBatch(kit, AnnoCNVGeneCore, CNVdb$ExAC_pLI, n = 2)[[1]]
## gdbList[[5]] <- AnnoCNVGeneRefGene2DDG2P(refGene[[1]], CNVdb$DECIPHER_DDG2P)

## geneTable <- SunMergeGenedb(gdbList)
