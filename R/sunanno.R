##' The CNV summary table from Sun.
##'
##' \itemize{
##'   \item REGION
##'   \item SIZE
##'   \item RATIO
##'   \item COPY
##'   \item TYPE
##'   \item CYTOBAND
##'   \item SAMPLE
##'   \item DERIVED
##' }
##'
##' @title Sun's summary CNV table
##' @inheritParams FilterCNS
##' @inheritParams FindCyto
##' @param sampleType The sample type, like "proband", "mother", "father".
##' @importFrom foreach foreach %do%
##' @importFrom iterators iter
##' @return A \code{data.frame} of Sun's cns format.
##' @examples
##' require('magrittr')
##' data(hg38cyto)
##'
##' cnsFile <- system.file('extdata', 'example.cnvkit', package = 'CNVanno')
##' cnsSummarySun <- cnsFile %>% read.cnvkit %>% FilterCNS %>% SunCNVTable(hg38cyto)
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @export
##'
SunCNVTable <- function(cns, cyto, sampleType = 'proband', n = 2) {

  itx <- iter(cns, by = 'row')
  REGION <- foreach(i = itx, .combine = c) %do% {
    return(paste0(c(i[1], ':', i[2], '-', i[3]), collapse = ''))
  }

  itx <- iter(cns, by = 'row')
  SIZE <- foreach(i = itx, .combine = c) %do% {
    return(round((i[3] - i[2]) / 100000, 1))
  } %>% unlist

  CYTOBAND <- FindCyto(cns, cyto, n = n)

  cnvTable <- data.frame(CNV = REGION,
                         SIZE = SIZE,
                         RATIO = cns[, 5],
                         COPYNUM = cns[, 6],
                         TYPE = ifelse(cns[, 6] > 2, 'gain', 'loss'),
                         CYTOBAND = CYTOBAND,
                         GENE = cns[, 4],
                         SAMPLE = rep(sampleType, nrow(cns)),
                         DEPTH = cns[, 7],
                         PROBES = cns[, 8],
                         WEIGHT = cns[, 9],
                         stringsAsFactors = FALSE)

  colnames(cnvTable)[2] <- 'SIZE(100kb)'

  return(cnvTable)
}

##' The CNV summary genes table for Sun.
##'
##' ## column names for OMIM
##' \itemize{
##'   \item CNV
##'   \item OverlapRelation
##'   \item Gene.Symbols
##'   \item Phenotypes
##'   \item Phenotypes.CN
##'   \item genetic???
##'   \item Mim.Number
##' }
##' ## column names for ClinGenTriHaplo
##' \itemize{
##'   \item Haploinsufficiency.Score
##'   \item Haploinsufficiency.PMID1
##'   \item Haploinsufficiency.PMID2
##'   \item Haploinsufficiency.PMID3
##'   \item Triplosensitivity.Score
##'   \item Triplosensitivity.PMID1
##'   \item Triplosensitivity.PMID2
##'   \item Triplosensitivity.PMID3
##'   \item Date.Last.Evaluated
##'   \item Loss.phenotype.OMIM.ID
##'   \item Triplosensitive.phenotype.OMIM.ID
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
##'   \item organ.specificity.list
##' }
##'
##' ## example column (has blanks)
##' ## chr3:162819884-162903642	CNVinGene	HDLCQ5	[High density lipoprotein cholesterol level QTL 5], 610761 (2)			610761	40				40				2014-03-12			0.020830666	0.003227266288348
##' @title Merge genes database
##' @param gdbList A list containing CNV genes database. 1st element is the OMIM database, 2nd is the ClinGenTriHaplo, 3rd is the DECIPHERHaplo, 4th is the ExACpLI.
##' @return A \code{data.frame}
##' @examples
##' require('magrittr')
##' data(CNVdb)
##' data(hg38cyto)
##'
##' cnsFile <- system.file('extdata', 'example.cnvkit', package = 'CNVanno')
##' cns <- cnsFile %>% read.cnvkit %>% FilterCNS
##'
##' gdbList <- list()
##' gdbList[[1]] <- AnnoCNVBatch(cns, AnnoCNVGeneCore, CNVdb$OMIMGRCh38)[[1]]
##' gdbList[[2]] <- AnnoCNVBatch(cns, AnnoCNVGeneCore, CNVdb$ClinGen_TriHaploGRCh38)[[1]]
##' gdbList[[3]] <- AnnoCNVBatch(cns, AnnoCNVGeneCore, CNVdb$DECIPHER_Haplo)[[1]]
##' gdbList[[4]] <- AnnoCNVBatch(cns, AnnoCNVGeneCore, CNVdb$ExAC_pLI)[[1]]
##' cnsRefGene <- AnnoCNVBatch(cns, AnnoCNVGeneCore, CNVdb$RefGeneGRCh37)
##' gdbList[[5]] <- AnnoCNVGeneRefGene2DDG2P(cnsRefGene[[1]], CNVdb$DECIPHER_DDG2P)
##'
##' geneTable <- SunMergeGenesdb(gdbList)
##'
##' cnsTable <- cns %>% SunCNVTable(hg38cyto, sampleType = 'proband', n = CORENUM)
##' CrossRegionGeneTable(cnsTable, geneTable)
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @importFrom magrittr %<>% %>%
##' @importFrom stringr str_trim
##' @rdname mergegene
##' @export
##'
SunMergeGenesdb <- function(gdbList) {

  ## process dbs
  gdb <- gdbList
  gdb[[1]] %<>% `[`(, c(1:5, 12:14, 6, 7, 10, 11), drop = FALSE)
  gdb[[2]] %<>% `[`(, c(1:5, 10:20), drop = FALSE)
  gdb[[3]] %<>% `[`(, c(1:5, 7), drop = FALSE)
  gdb[[4]] %<>% `[`(, c(1:5, 22), drop = FALSE)
  gdb[[5]] %<>% `[`(, c(1:5, 8:13), drop = FALSE)

  res <- gdb[[1]] %>% cbind.data.frame(.,
                                       matrix('',
                                              nrow = nrow(gdbList[[1]]),
                                              ncol = 19),
                                       stringsAsFactors = FALSE)
  cnames <- c(colnames(gdb[[1]]),
              colnames(gdb[[2]])[-1:-5],
              colnames(gdb[[3]])[-1:-5],
              colnames(gdb[[4]])[-1:-5],
              colnames(gdb[[5]])[-1:-5])

  colnames(res) <- cnames

  ## OMIMkey
  k <- gdbList[[1]][, 'Gene.Symbols']

  ## search ClinGenTriHaplo
  m <- gdb[[2]]
  tk <- gdbList[[2]][, 'Gene.Symbol']
  for(i in seq_along(tk)) {
    midx <- MatchString(tk[i], k)
    if (is.na(midx)) {
      ## add loc
      addEmpty <- c(m[i, 1:5],
                    tk[i],
                    rep('', 6),
                    m[i, -1:-5],
                    rep('', 8))

      names(addEmpty) <- cnames

      res %<>% rbind.data.frame(addEmpty, stringsAsFactors = FALSE)
    } else {
      res[midx, 13:23] <- m[i, -1:-5]
    }
  }

  ## search DECIPHERHaplo
  k <- res[, 'Gene.Symbols']
  m <- gdb[[3]]
  tk <- gdbList[[3]][, 'name'] %>%
    strsplit(split = '|', fixed = TRUE) %>%
    sapply(`[`, 1)
  for(i in seq_along(tk)) {
    midx <- MatchString(tk[i], k)
    if (is.na(midx)) {
      ## add loc
      addEmpty <- c(m[i, 1:5],
                    tk[i],
                    rep('', 17),
                    m[i, -1:-5],
                    rep('', 7))

      names(addEmpty) <- cnames

      res %<>% rbind.data.frame(addEmpty, stringsAsFactors = FALSE)
    } else {
      res[midx, 24] <- m[i, -1:-5]
    }
  }

  ## search ExACpLI
  k <- res[, 'Gene.Symbols']
  m <- gdb[[4]]
  tk <- gdbList[[4]][, 'gene']
  for(i in seq_along(tk)) {
    midx <- MatchString(tk[i], k)
    if (is.na(midx)) {
      ## add loc
      addEmpty <- c(m[i, 1:5],
                    tk[i],
                    rep('', 18),
                    m[i, -1:-5],
                    rep('', 6))

      names(addEmpty) <- cnames

      res %<>% rbind.data.frame(addEmpty, stringsAsFactors = FALSE)
    } else {
      res[midx, 25] <- m[i, -1:-5]
    }
  }

  ## search DDG2P
  k <- res[, 'Gene.Symbols']
  m <- gdb[[5]]
  tk <- gdbList[[5]][, 'gene.symbol']
  for(i in seq_along(tk)) {
    midx <- MatchString(tk[i], k)
    if (is.na(midx)) {
      ## add loc
      addEmpty <- c(m[i, 1:5],
                    tk[i],
                    rep('', 19),
                    m[i, -1:-5])

      names(addEmpty) <- cnames

      res %<>% rbind.data.frame(addEmpty, stringsAsFactors = FALSE)
    } else {
      res[midx, 26:31] <- m[i, -1:-5]
    }
  }

  res %<>% GenedbPost

  return(res)
}


##' @param inputStr A \code{string}.
##' @param searchStr A \code{string vector}.
##' @return A processed \code{data.frame}.
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @rdname mergegene
##' @importFrom magrittr %>%
##' @importFrom stringr str_trim
##' @keywords internal
##'
MatchString <- function(inputStr, searchStr) {
  sRes <- sapply(searchStr, function(x) {
    sStr <- x %>%
      strsplit(split = ',', fixed = TRUE) %>%
      unlist %>%
      str_trim

    return(inputStr %in% sStr)
  })

  return(which(sRes)[1])
  ## return(which(str_detect(searchStr, inputStr))[1])
}


##' @param mergeTable A \code{data.frame} of processed gene table.
##' @return \code{NULL} or the index
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @importFrom stringr str_extract
##' @importFrom magrittr %<>% %>%
##' @rdname mergegene
##' @keywords internal
##'
GenedbPost <- function(mergeTable) {
  ## step1: sort chromosomes
  oGenes <- mergeTable[, 'Chromosome'] %>%
    str_extract('\\d+') %>%
    as.numeric %>%
    order
  mergeTable <- mergeTable[oGenes, ]

  ## step2: select genes
  mergeTable[, 'Type'] %<>% ifelse(. == '', 'gene', .)
  mergeTable <- mergeTable[mergeTable[, 'Type'] == 'gene', , drop = FALSE]

  ## step3: remove column 'Gene.Name' and 'Type'
  mergeTable <- mergeTable[, c(-7, -12)]

  return(mergeTable)
}


##' @param regionTable A \code{data.frame} of processed region table, which often comes from the \code{SunCNVTable} functions.
##' @param geneTable A \code{data.frame} of gene table, which often comes from the \code{SunMergeGenesdb} functions.
##' @return A list containing two \code{data.frame} represents the processed \code{regionTable} and \code{geneTable} respectively.
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @importFrom magrittr %<>%
##' @importFrom stats aggregate
##' @importFrom stringr str_detect
##' @rdname mergegene
##' @export
##'
CrossRegionGeneTable <- function(regionTable, geneTable) {

  ## step 1 initial `regionTable`
  regionTable <- cbind.data.frame(regionTable,
                                  matrix('',
                                         nrow = nrow(regionTable),
                                         ncol = 24),
                                  stringsAsFactors = FALSE)
  colnames(regionTable)[-1:-11] <- colnames(geneTable)[-1:-5]
  colnames(regionTable)[12] <- 'ISCA-ID'

  ## step 2 select 'ISCA-37488' and merge gene symbols
  iscaLog <- str_detect(geneTable[, 'Gene.Symbols'], 'ISCA-')
  iscaMat <- geneTable[iscaLog, -2:-5, drop = FALSE]
  colnames(iscaMat)[2] <- 'ISCA-ID'
  geneTable <- geneTable[!iscaLog, , drop = FALSE]

  unCNV <- aggregate(geneTable[, 'Gene.Symbols'],
                     by = list(geneTable[, 'CNV']),
                     function(x) {
                       return(paste(unique(x), collapse = '|'))
                     },
                     drop = FALSE)
  ## process gene
  regionTable[, 'GENE'] <- unCNV[match(regionTable[, 'CNV'], unCNV[, 1]), 2]
  regionTable[, 'GENE'] %<>% ifelse(is.na(.), '', .)

  ## step3 process isca mat
  if (sum(iscaLog) > 0) {
    iscaMat <- aggregate(iscaMat[, -1],
                         by = list(iscaMat[, 1]),
                         paste, collapse = '|',
                         drop = FALSE)

    iscaIdx <- match(iscaMat[, 1], regionTable[, 1])
    regionTable[iscaIdx, -1:-11] <- iscaMat[, -1]
  } else {}

  res <- list(regionTable = regionTable,
              geneTable = geneTable)

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
##' @examples
##' require('magrittr')
##' data(CNVdb)
##' data(hg38cyto)
##'
##' cnsFile <- system.file('extdata', 'example.cnvkit', package = 'CNVanno')
##' cns <- cnsFile %>% read.cnvkit %>% FilterCNS
##'
##' cnsRefGene <- AnnoCNVBatch(cns, AnnoCNVGeneCore, CNVdb$RefGeneGRCh37)
##' cnsOMIM <- AnnoCNVGeneRefGene2OMIM(cnsRefGene[[1]], CNVdb$OMIMGRCh38)
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @export
##'
AnnoCNVGeneRefGene2OMIM <- function(refgeneTable, OMIMdb) {

  ## step 1: connect RefGene and OMIMdb with Entrez ID
  refTable <- refgeneTable[, c(1:5, 8), drop = FALSE]
  refTable <- aggregate(refTable[, -6], by = list(refTable[, 6]), `[[`, 1)
  colnames(refTable)[1:5] <- c('RefGene',
                               'CNV',
                               'Chromosome',
                               'Genomic.Position.Start',
                               'Genomic.Position.End')
  OMIMdb <- OMIMdb[, -1:-3]

  res <- merge.data.frame(refTable, OMIMdb, by.x = 'RefGene', by.y = 'Entrez.Gene.ID')
  res <- res[, c(2:15, 1, 16:18), drop = FALSE]
  colnames(res)[15] <- 'Entrez.Gene.ID'

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
##' data(hg38cyto)
##'
##' cnsFile <- system.file('extdata', 'example.cnvkit', package = 'CNVanno')
##' cns <- cnsFile %>% read.cnvkit %>% FilterCNS
##'
##' cnsRefGene <- AnnoCNVBatch(cns, AnnoCNVGeneCore, CNVdb$RefGeneGRCh37)
##' cnsDDG2P <- AnnoCNVGeneRefGene2DDG2P(cnsRefGene[[1]], CNVdb$DECIPHER_DDG2P)
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @export
##'
AnnoCNVGeneRefGene2DDG2P <- function(refgeneTable, DDG2Pdb) {

  ## step 1: connect RefGene and OMIMdb with Entrez ID
  refTable <- refgeneTable[, c(1:5, 8), drop = FALSE]
  refTable <- aggregate(refTable[, -6], by = list(refTable[, 6]), `[[`, 1)
  colnames(refTable)[1:5] <- c('RefGene',
                               'CNV',
                               'Chromosome',
                               'Genomic.Position.Start',
                               'Genomic.Position.End')

  DDG2Pdb <- DDG2Pdb[, c(1, 2, 4, 6:10, 15)]

  res <- merge.data.frame(refTable, DDG2Pdb, by.x = 'RefGene', by.y = 'entrezgene')
  res <- res[, -1, drop = FALSE]

  return(res)
}



