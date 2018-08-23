##' Cytobands of hg38 genome from the UCSC goldenpath
##'
##' The default hg38 cytobands.
##'
##' @docType data
##' @name hg38cyto
##' @format A \code{tbl_df}
##' @references \href{http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/}{UCSC goldenpath}
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##'
NULL


##' Cytobands of hg19 genome from the UCSC goldenpath
##'
##' The default hg19 cytobands.
##'
##' @docType data
##' @name hg19cyto
##' @format A \code{tbl_df}
##' @references \href{http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/}{UCSC goldenpath}
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##'
NULL


##' Blacklist of hg19 genome.
##'
##' Combination of low mapbility regions and blacklist from the 10X genome.
##'
##' @docType data
##' @name hg19bl
##' @format A \code{tbl_df}
##' @references \href{ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeDacMapabilityConsensusExcludable.bed.gz}{mapbility_concensus}
##' @references \href{ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeDukeMapabilityRegionsExcludable.bed.gz}{mapbility_region}
##' @references \href{http://cf.10xgenomics.com/supp/genome/hg19/sv_blacklist.bed}{10X genome blacklist}
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##'
NULL


##' Blacklist of hg38 genome.
##'
##' Combination of low mapbility regions and blacklist from the 10X genome.
##'
##' @docType data
##' @name hg38bl
##' @format A \code{tbl_df}
##' @references \href{http://cf.10xgenomics.com/supp/genome/hg38/sv_blacklist.bed}{10X genome blacklist}
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##'
NULL


##' CNVnator filtered example file
##'
##' ## CNVanno version 0.2.0
##' nator <- system.file('extdata', 'example.cnvnator', package = 'CNVanno') %>%
##'   read_cnvnator %>%
##'   filter_cnvnator %>%
##'   Segment(gap = 10L) %>%
##'   FilterBlacklist(bl_cytoband(hg19cyto), overlaprate = 0.5, shortlen = 1000L, gap = 0L, n = 2) %>%
##'   FilterBlacklist(hg19bl, overlaprate = 0.5, shortlen = 1000L, gap = 100000L, n = 2)
##'
##'
##' @docType data
##' @name nator
##' @format A \code{coreCNV} object
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##'
NULL


##' CNVkit filtered example file
##'
##' ## CNVanno version 0.2.0
##' kit <- system.file('extdata', 'example.cnvkit', package = 'CNVanno') %>%
##'   read_cnvkit %>%
##'   filter_cnvkit %>%
##'   Segment(gap = 10L)%>%
##'   FilterBlacklist(bl_cytoband(hg19cyto), overlaprate = 0.5, shortlen = 1000L, gap = 0L, n = 2) %>%
##'   FilterBlacklist(hg19bl, overlaprate = 0.5, shortlen = 1000L, gap = 1000000L, n = 2)
##'
##' @docType data
##' @name kit
##' @format A \code{coreCNV} object
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##'
NULL


##' Annotation databases of CNV
##'
##' \itemize{
##'   \item ClinGenGRCh38 2017-04-01
##'   \item ClinGenGRCh37 2017-04-01
##'   \item ClinGenNCBI36 2017-04-01
##'   \item ClinVarGRCh38
##'   \item ClinVarGRCh37
##'   \item ClinVarNCBI36
##'   \item ClinGen_TriHaploGRCh38 2018-07-14
##'   \item ClinGen_TriHaploGRCh37 2018-07-14
##'   \item DECIPHER_Haplo
##'   \item DECIPHER_POPULATION
##'   \item DECIPHER_DDG2P 2018-07-14
##'   \item DGV
##'   \item ExAC_POPULATION
##'   \item ExAC_pLI
##'   \item OMIMGRCh38 2018-07-13
##'   \item RefGeneGRCh37
##' }
##'
##' @docType data
##' @name CNVdb
##' @format A \code{list}
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##'
NULL
