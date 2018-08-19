##' @include AllClasses.R AllGenerics.R region.R
NULL


##' Filter the CNV file generated from \code{CNVkit} and \code{CNVnator}
##'
##' For the \code{CNVkit}
##' \itemize{
##'   \item 1. filter the \code{cns} with cn number.
##'   \item 2. filter low sequence depth.
##'   \item 3. filter X/Y chromosomes
##' }
##'
##' For the \code{CNVnator}
##' \itemize{
##'   \item 1. filter X/Y chromosomes
##' }
##'
##' @title CNV filter
##' @param cngain A \code{numeric} value as the gain copy number threshold, and the default value is 2.
##' @param cnloss A \code{numeric} value as the loss copy number threshold, and the default value is 2.
##' @param dep A \code{numeric} value as the sequence depth, and the default value is 0.01
##' @param sexchrom A \code{logic} value whether filter sex chromosomes (X and Y), the default value is \code{TRUE}.
##' @param rawkit The standard \code{CNVkit} output as an \code{RawCNV} object.
##' @return A filtered \code{RawCNV} object.
##' @examples
##' require('magrittr')
##'
##' kitf <- system.file('extdata', 'example.cnvkit', package = 'CNVanno') %>%
##'   read_cnvkit %>%
##'   filter_cnvkit
##'
##' natorf <- system.file('extdata', 'example.cnvnator', package = 'CNVanno') %>%
##'   read_cnvnator %>%
##'   filter_cnvnator
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @importFrom dplyr transmute filter
##' @importFrom magrittr %>% %<>%
##' @rdname filterraw
##' @export
##'
filter_cnvkit <- function(rawkit, cngain = 2, cnloss = 2, dep = 0.01, sexchrom = TRUE) {

  ## step1: filter gain and loss cnv
  ## step2: filter low depth

  l <- rawkit@params %>%
    transmute((cn > cngain | cn < cnloss) & depth > dep) %>%
    unlist

  ## step3: filter sex chromosomes
  if (sexchrom) {
    lsex <- rawkit@rawCNV %>%
      transmute(chromosome != 'chrX' & chromosome != 'chrY') %>%
      unlist
    l %<>% `&`(lsex)
  } else {}

  cnvkitf <- new('RawCNV',
                 rawCNV = filter(rawkit@rawCNV, l),
                 params = filter(rawkit@params, l),
                 method = rawkit@method)

  return(cnvkitf)

}


##' @param rawnator The standard \code{CNVnator} output as an \code{RawCNV} object.
##' @inheritParams filter_cnvkit
##' @importFrom magrittr  %>%
##' @importFrom dplyr select mutate
##' @rdname filterraw
##' @export
##'
filter_cnvnator <- function(rawnator, sexchrom = TRUE) {

  ## step1: filter sex chromosomes
  if (sexchrom) {
    l <- rawnator@rawCNV %>%
      transmute(chromosome != 'chrX' & chromosome != 'chrY') %>%
      unlist
  } else {}

  cnvnatorf <- new('RawCNV',
                   rawCNV = filter(rawnator@rawCNV, l),
                   params = filter(rawnator@params, l),
                   method = rawnator@method)

  return(cnvnatorf)
}



##' Filtering \code{CoreCNV}.
##'
##' Filer the \code{CoreCNV} according black lists.
##'
##' @title Filter \code{CoreCNV} according to a give blacklist
##' @inheritParams FilterBlacklist
##' @return A \code{CoreCNV} object.
##' @examples
##' library('magrittr')
##' data(hg19cyto)
##'
##' nator <- system.file('extdata', 'example.cnvnator', package = 'CNVanno') %>%
##'   read_cnvnator %>%
##'   filter_cnvnator %>%
##'   Segment(gap = 10L)
##'
##' ## natorf <- FilterBlacklist(nator, bl_cytoband(hg19cyto), overlaprate = 0.5, n = 1)
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @importFrom doParallel registerDoParallel stopImplicitCluster
##' @importFrom foreach foreach %dopar%
##' @importFrom iterators iter
##' @importFrom magrittr %>%
##' @importFrom dplyr bind_rows
##' @importFrom methods new
##' @references \href{http://penncnv.openbioinformatics.org/en/latest/user-guide/annotation/#filtering-cnv-calls-by-user-specified-criteria}{cytoband extend}
##' @rdname FilterBlacklist-methods
##' @exportMethod FilterBlacklist
##'
setMethod(f = 'FilterBlacklist',
          signature = c(core = 'CoreCNV', blacklist = 'tbl_df', overlaprate = 'numeric'),
          definition = function(core, blacklist, overlaprate, n, ...) {
            core <- core@coreCNV

            registerDoParallel(cores = n)
            itx <- iter(core, by = 'row')

            coreFilter <- foreach(i = itx) %dopar% {
              eachCore <- filterRow_(i,
                                     blacklist = blacklist,
                                     overlaprate = overlaprate)
              return(eachCore)
            }

            ## stop multiple cores
            stopImplicitCluster()

            coreFilter %<>%
              bind_rows %>%
              new('CoreCNV', coreCNV = .)

            return(coreFilter)
          })



##' Filtering internal functions.
##'
##' \itemize{
##'   \item \code{filterSeg_()}: For the blacklist (<= `overlaprate`) in the same chromosome, segment the regions of input `regionf`. If the `regionf` wrap up the `bl` (====~~~~===), `regionf` is cut into two segments. The minimum left and right regions are kept.
##'   \item \code{filterRow_()}: Filter input one cnv (one row).
##' }
##'
##' @title Internal functions for filtering
##' @inheritParams OverlapRegionRate
##' @inheritParams segMergeType_
##' @param rateMat The output of \code{OverlapRegionRate()} in this package.
##' @return A \code{tbl_df} of filtered fragments.
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @importFrom magrittr %<>% %>%
##' @importFrom dplyr filter mutate distinct bind_rows rename everything
##' @rdname filterutility
##' @keywords internal
##'
filterSeg_ <- function(regionf, rateMat, chr, type) {
  # Example:
  ## rM <- tibble(start = c(9L, 2L, 5L, 8L, 3L, 1L, 15L, 10L), end = c(15L, 6L, 6L, 9L, 20L, 4L, 32L, 9L)) %>% SortRegion %>% ReduceRegion(gap = 0L)
  ## rf <- c(4L, 10L)
  ## rMf <- OverlapRegionRate(rf, rM) %>% filter(fRate > 0 & fRate <= 0.5)
  ## filterSeg_(rf, rMf, chr = 'chr6', type = 'gain')

  ## step1: split regions
  minf <- min(regionf)
  maxf <- max(regionf)
  ## intersect =====~~~~~ or ~~~~=====
  rateMatInter <- rateMat %>%
    filter(tRate < 1) %>% ## select intersect regions
    mutate(segstart = if_else(minf < maxstart, minf, minend)) %>%
    mutate(segend = if_else(maxf > minend, maxf, maxstart))

  ## in ====~~~~===
  rateMatCoverTemp  <- rateMat %>%
    filter(tRate == 1)

  rateMatCoverLeft <- rateMatCoverTemp %>%
    mutate(segstart = minf) %>%
    mutate(segend = maxstart)

  rateMatCoverRight <- rateMatCoverTemp %>%
    mutate(segstart = minend) %>%
    mutate(segend = maxf)

  seg <- bind_rows(list(rateMatInter,
                        rateMatCoverLeft,
                        rateMatCoverRight)) %>%
    select(segstart, segend)

  ## step2: select min regions on left and right, respectively
  leftSeg <- seg %>%
    filter(segstart == minf) %>%
    distinct() %>%
    filter(segend == min(segend))

  rightSeg <- seg %>%
    filter(segend == maxf) %>%
    distinct() %>%
    filter(segstart == max(segstart))

  seg <- bind_rows(list(leftSeg, rightSeg)) %>%
    mutate(chromosome = chr, type = type) %>%
    rename(start = segstart, end = segend) %>%
    select(chromosome, everything())

  return(seg)

}


##' @inheritParams filterSeg_
##' @inheritParams OverlapRegionRate
##' @inheritParams FilterBlacklist
##' @importFrom magrittr %<>% %>%
##' @importFrom dplyr filter bind_cols
##' @rdname filterutility
##' @keywords internal
##'
filterSegCover_ <- function(regionf, rateMat, overlaprate, shortlen) {
  ## step1: check the sum rate of intersection regions is larger than `overlaprate`
  rateMat %<>% filter(tRate == 1)

  sumRate <- rateMat %>%
    select(fRate) %>%
    sum

  if(sumRate >= overlaprate) {
    return(NULL)
  } else {
    ## step2: get gaps
    start <- (rateMat$minend + 1L) %>%
      c(min(regionf), .)
    end <- (rateMat$maxstart - 1L) %>%
      c(max(regionf))

    coverRegion <- list(start = start, end = end) %>%
      bind_cols %>%
      filter(start - end + 1 > shortlen) ## filter too short regions
  }

  return(coverRegion)
}


##' @inheritParams filterSeg_
##' @inheritParams OverlapRegionRate
##' @inheritParams FilterBlacklist
##' @importFrom magrittr %<>% %>%
##' @importFrom dplyr filter bind_cols
##' @rdname filterutility
##' @keywords internal
##'
filterSegInter_ <- function(regionf, rateMat, shortlen) {

  minf <- min(regionf)
  maxf <- max(regionf)

  ## intersect =====~~~~~ or ~~~~==== (only these tow cases)
  ## no need to reduce
  interRegion <- rateMat %>%
    filter(tRate < 1) %>% ## select intersect regions
    mutate(start = if_else(minf < maxstart, minf, minend)) %>%
    mutate(end = if_else(maxf > minend, maxf, maxstart)) %>% ## start/end
    select(star:end) %>%
    filter(start - end + 1 > shortlen)

  return(interRegion)

}


##' @inheritParams FilterBlacklist
##' @param corerow A row of the CNV in a \code{tbl_df} form.
##' @importFrom magrittr %<>% %>%
##' @importFrom dplyr select filter
##' @rdname filterutility
##' @keywords internal
##'
filterRow_ <- function(corerow, blacklist, overlaprate) {

  blacklist %<>%
    filter(chromosome == corerow$chromosome)

  if (nrow(blacklist) == 0) {
    ## case 1: no same chromosomes in blacklist
    return(corerow)
  } else {}

  rateMat <- blacklist %<>%
    select(start, end) %>%
    OverlapRegionRate(select(corerow, start, end), .)

  frate <- rateMat$fRate
  if (sum(frate >= overlaprate) > 0) {
    ## case 2: has overlap region with > overlaprate
    return(NULL)
  }
  else if (all(frate == 0)) {
    ## case 3: 0 overlap regions
    return(corerow)
  }
  else {
    ## case 4: has > 0 & < overlaprate overlap regions
    seg <- rateMat %>%
      filter(fRate > 0) %>%
      filterSeg_(select(corerow, start, end),
                 .,
                 corerow$chromosome,
                 corerow$method)
    return(seg)
  }

}



##' Filter the cytoband.
##'
##' Construct the blacklist of cytoband. The \code{extend} can be set as 100,000 500,000 or 1000,000. The "acen", "gvar", and "stalk" will be kept in the black list.
##'
##' @title Blacklist of cytoband
##' @inheritParams Cytoband
##' @inheritParams OverlapRegion
##' @return A \code{tbl_df} represents the black list of cytoband.
##' @examples
##' data(hg19cyto)
##'
##' bl_cytoband(hg19cyto, extend = 5e5L)
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @importFrom magrittr %>%
##' @importFrom dplyr filter mutate
##' @export
##'
bl_cytoband <- function(cyto, extend = 5e5L) {

  ## step1: extend bl regions
  cytobl <- cyto %>%
    filter(color %in% c('acen', 'gvar', 'stalk')) %>% ## filter
    mutate(start = if_else(start > extend, start - extend, 0L)) %>%
    mutate(end = end + extend) %>% ## extend
    SortRegionChr %>% ## sort cyto
    ReduceRegionChr(gap = 0L) %>% ## reduce regions
    SortRegionChr ## sort again

  return(cytobl)
}
