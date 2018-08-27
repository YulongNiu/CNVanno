## ##' The CNV summary table from Sun.
## ##'
## ##' \itemize{
## ##'   \item REGION
## ##'   \item SIZE
## ##'   \item RATIO
## ##'   \item COPY
## ##'   \item TYPE
## ##'   \item CYTOBAND
## ##'   \item SAMPLE
## ##'   \item DERIVED
## ##' }
## ##'
## ##' @title Sun's summary CNV table
## ##' @inheritParams filter.cnvkit
## ##' @inheritParams Cytoband
## ##' @param sampleType The sample type, like "proband", "mother", "father".
## ##' @importFrom foreach foreach %do%
## ##' @importFrom iterators iter
## ##' @return A \code{data.frame} of Sun's cns format.
## ##' @examples
## ##' require('magrittr')
## ##' data(hg38cyto)
## ##'
## ##' cnsFile <- system.file('extdata', 'example.cnvkit', package = 'CNVanno')
## ##' cnsSummarySun <- cnsFile %>% read.cnvkit %>% filter.cnvkit %>% SunCNVTable(hg38cyto)
## ##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
## ##' @export
## ##'
## SunCNVTable <- function(cns, cyto, sampleType = 'proband', n = 2) {

##   itx <- iter(cns, by = 'row')
##   REGION <- foreach(i = itx, .combine = c) %do% {
##     return(paste0(c(i[1], ':', i[2], '-', i[3]), collapse = ''))
##   }

##   itx <- iter(cns, by = 'row')
##   SIZE <- foreach(i = itx, .combine = c) %do% {
##     return(round((i[3] - i[2]) / 100000, 1))
##   } %>% unlist

##   CYTOBAND <- FindCyto(cns, cyto, n = n)

##   cnvTable <- data.frame(CNV = REGION,
##                          SIZE = SIZE,
##                          RATIO = cns[, 5],
##                          COPYNUM = cns[, 6],
##                          TYPE = ifelse(cns[, 6] > 2, 'gain', 'loss'),
##                          CYTOBAND = CYTOBAND,
##                          GENE = cns[, 4],
##                          SAMPLE = rep(sampleType, nrow(cns)),
##                          DEPTH = cns[, 7],
##                          PROBES = cns[, 8],
##                          WEIGHT = cns[, 9],
##                          stringsAsFactors = FALSE)

##   colnames(cnvTable)[2] <- 'SIZE(100kb)'

##   return(cnvTable)
## }



##' The CNV summary genes table.
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
##'   \item gene.mim
##'   \item disease.name
##'   \item disease.mim
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
##' geneTable <- SunMergeGenedb(gdbList)
##'
##' ## cnsTable <- cns %>% SunCNVTable(hg38cyto, sampleType = 'proband', n = CORENUM)
##' ## CrossRegionGeneTable(cnsTable, geneTable)
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @importFrom magrittr %<>% %>%
##' @importFrom dplyr rename select group_by distinct ungroup mutate
##' @rdname mergegene
##' @export
##'
SunMergeGenedb <- function(gdbList) {

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
    select(CNV:overlap_relation, gene_symbol, gene.mim:phenotypes) %>%
    group_by(CNV) %>%
    distinct(gene_symbol, .keep_all = TRUE) %>%
    ungroup

  ## res %<>% GenedbPost
  gdb <- Reduce(mergeTwoGenedb_, gdbList)

  return(gdb)
}


##' Internal Sun's template annotation functions.
##'
##' \itemize{
##'   \item \code{mergeTwoGenedb_()}: Merge two gene results.
##' }
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

  return(gdb)
}


## ##' @param mergeTable A \code{data.frame} of processed gene table.
## ##' @return \code{NULL} or the index
## ##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
## ##' @importFrom stringr str_extract
## ##' @importFrom magrittr %<>% %>%
## ##' @rdname mergegene
## ##' @keywords internal
## ##'
## GenedbPost <- function(mergeTable) {
##   ## step1: sort chromosomes
##   oGenes <- mergeTable[, 'Chromosome'] %>%
##     str_extract('\\d+') %>%
##     as.numeric %>%
##     order
##   mergeTable <- mergeTable[oGenes, ]

##   ## step2: select genes
##   mergeTable[, 'Type'] %<>% ifelse(. == '', 'gene', .)
##   mergeTable <- mergeTable[mergeTable[, 'Type'] == 'gene', , drop = FALSE]

##   ## step3: remove column 'Gene.Name' and 'Type'
##   mergeTable <- mergeTable[, c(-7, -12)]

##   return(mergeTable)
## }


## ##' @param regionTable A \code{data.frame} of processed region table, which often comes from the \code{SunCNVTable} functions.
## ##' @param geneTable A \code{data.frame} of gene table, which often comes from the \code{SunMergeGenesdb} functions.
## ##' @return A list containing two \code{data.frame} represents the processed \code{regionTable} and \code{geneTable} respectively.
## ##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
## ##' @importFrom magrittr %<>%
## ##' @importFrom stats aggregate
## ##' @importFrom stringr str_detect
## ##' @rdname mergegene
## ##' @export
## ##'
## CrossRegionGeneTable <- function(regionTable, geneTable) {

##   ## step 1 initial `regionTable`
##   regionTable <- cbind.data.frame(regionTable,
##                                   matrix('',
##                                          nrow = nrow(regionTable),
##                                          ncol = 24),
##                                   stringsAsFactors = FALSE)
##   colnames(regionTable)[-1:-11] <- colnames(geneTable)[-1:-5]
##   colnames(regionTable)[12] <- 'ISCA-ID'

##   ## step 2 select 'ISCA-37488' and merge gene symbols
##   iscaLog <- str_detect(geneTable[, 'Gene.Symbols'], 'ISCA-')
##   iscaMat <- geneTable[iscaLog, -2:-5, drop = FALSE]
##   colnames(iscaMat)[2] <- 'ISCA-ID'
##   geneTable <- geneTable[!iscaLog, , drop = FALSE]

##   unCNV <- aggregate(geneTable[, 'Gene.Symbols'],
##                      by = list(geneTable[, 'CNV']),
##                      function(x) {
##                        return(paste(unique(x), collapse = '|'))
##                      },
##                      drop = FALSE)
##   ## process gene
##   regionTable[, 'GENE'] <- unCNV[match(regionTable[, 'CNV'], unCNV[, 1]), 2]
##   regionTable[, 'GENE'] %<>% ifelse(is.na(.), '', .)

##   ## step3 process isca mat
##   if (sum(iscaLog) > 0) {
##     iscaMat <- aggregate(iscaMat[, -1],
##                          by = list(iscaMat[, 1]),
##                          paste, collapse = '|',
##                          drop = FALSE)

##     iscaIdx <- match(iscaMat[, 1], regionTable[, 1])
##     regionTable[iscaIdx, -1:-11] <- iscaMat[, -1]
##   } else {}

##   res <- list(regionTable = regionTable,
##               geneTable = geneTable)

##   return(res)

## }

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

## geneTable <- SunMergeGenesdb(gdbList)
