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
##' Filer the \code{CoreCNV} according black lists. Multiple black lists can be applied as the example show below.
##'
##' @title Filter \code{CoreCNV} according to a give blacklist
##' @inheritParams FilterBlacklist
##' @return A \code{CoreCNV} object.
##' @examples
##' library('magrittr')
##' data(hg19cyto)
##' data(hg19bl)
##'
##' nator <- system.file('extdata', 'example.cnvnator', package = 'CNVanno') %>%
##'   read_cnvnator %>%
##'   filter_cnvnator %>%
##'   Segment(gap = 10L)
##'
##' ## filter based on cytoband blacklist
##' hg19cytobl <- bl_cytoband(hg19cyto)
##' natorf <- FilterBlacklist(nator, hg19cytobl, overlaprate = 0.5, shortlen = 10L, gap = 0L, n = 2)
##'
##' ## more filter based on pre-built blacklist
##' natorf <- FilterBlacklist(natorf, hg19bl, overlaprate = 0.5, shortlen = 10L, gap = 1000L, n = 2)
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @importFrom doParallel registerDoParallel stopImplicitCluster
##' @importFrom foreach foreach %dopar%
##' @importFrom iterators iter
##' @importFrom magrittr %>%
##' @importFrom dplyr bind_rows do group_by ungroup select everything
##' @importFrom tibble tibble
##' @importFrom methods new
##' @references \href{http://penncnv.openbioinformatics.org/en/latest/user-guide/annotation/#filtering-cnv-calls-by-user-specified-criteria}{cytoband extend}
##' @rdname FilterBlacklist-methods
##' @exportMethod FilterBlacklist
##'
setMethod(f = 'FilterBlacklist',
          signature = c(core = 'CoreCNV', blacklist = 'tbl_df', overlaprate = 'numeric', shortlen = 'integer', gap = 'integer'),
          definition = function(core, blacklist, overlaprate, shortlen, gap, n, ...) {
            core <- core@coreCNV

            registerDoParallel(cores = n)
            itx <- iter(core, by = 'row')

            coreFilter <- foreach(i = itx) %dopar% {
              eachCore <- filterRow_(i,
                                     blacklist = blacklist,
                                     overlaprate = overlaprate,
                                     shortlen = shortlen)
              return(eachCore)
            }

            ## stop multiple cores
            stopImplicitCluster()

            coreFilter %<>%
              bind_rows %>%
              group_by(chromosome, type, method) %>% ## sort reduce
              do(SortRegion(tibble(start = .$start, end = .$end))) %>%
              do(ReduceRegion(tibble(start = .$start, end = .$end), gap = gap)) %>%
              select(chromosome, start:end, everything()) %>% ## colmuns in right order
              ungroup %>%
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
##' @inheritParams FilterBlacklist
##' @param rateMat The output of \code{OverlapRegionRate()} in this package.
##' @param chr A \code{string} indicating the chromosome, like "chr1", "chr2".
##' @param type A code{string} "gain", "loss", or "normal"
##' @param method A \code{string} indicating the CNV calling method.
##' @return A \code{tbl_df} of filtered fragments.
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @importFrom magrittr %<>% %>%
##' @importFrom dplyr filter mutate distinct bind_rows everything
##' @rdname filterutility
##' @keywords internal
##'
filterSeg_ <- function(regionf, rateMat, overlaprate, shortlen, chr, type, method) {
  # Example:
  ## rM <- tibble(start = c(1L, 4L, 8L, 14L, 18L, 25L), end = c(2L, 6L, 10L, 16L, 27L, 30L))
  ## rf <- c(5L, 20L)
  ## rMf <- OverlapRegionRate(rf, rM) %>% filter(fRate > 0 & fRate <= 0.5)
  ## filterSegCover_(rf, rMf, 0.6, 1L)
  ## filterSegInter_(rf, rMf, 0.6)
  ## filterSeg_(rf, rMf, 0.6, 1L, chr = 'chr6', type = 'gain')

  ## step 1: split regions
  ## case 1: intersect =====~~~~~ or ~~~~=====
  interRegion <- filterSegInter_(regionf = regionf,
                                 rateMat = rateMat,
                                 shortlen = shortlen)

  ## case 2: in ====~~~~===
  coverRegion <- filterSegCover_(regionf = regionf,
                                 rateMat = rateMat,
                                 overlaprate = overlaprate,
                                 shortlen = shortlen)

  keepRegion <- list(interRegion, coverRegion) %>%
    bind_rows

  ## step2: select min regions on left and right, respectively
  leftSeg <- keepRegion %>%
    filter(start == min(regionf)) %>%
    distinct %>%
    filter(end == min(end))

  rightSeg <- keepRegion %>%
    filter(end == max(regionf)) %>%
    distinct %>%
    filter(start == max(start))

  midSeg <- keepRegion %>%
    filter(start != min(regionf) & end != max(regionf))

  seg <- bind_rows(list(leftSeg, rightSeg, midSeg)) %>%
    SortRegion %>% ## sort
    mutate(chromosome = chr, type = type, method = method) %>% ## other columns
    select(chromosome, everything())

  return(seg)

}


##' @inheritParams filterSeg_
##' @importFrom magrittr %<>% %>%
##' @importFrom dplyr filter select bind_cols
##' @rdname filterutility
##' @keywords internal
##'
filterSegCover_ <- function(regionf, rateMat, overlaprate, shortlen) {
  ## step1: check the sum rate of intersection regions is larger than `overlaprate`

  rateMat %<>% filter(tRate == 1)

  sumRate <- rateMat %>%
    select(fRate) %>%
    unlist %>%
    sum

  ## sum cover region rate > overlaprate is filtered
  ## no tRate == 1 is filtered
  if (sumRate > 0 & sumRate < overlaprate) {
    ## step2: get gaps
    start <- (rateMat$minend + 1L) %>%
      c(min(regionf), .)
    end <- (rateMat$maxstart - 1L) %>%
      c(max(regionf))

    coverRegion <- list(start = start, end = end) %>%
      bind_cols %>%
      filter(end - start + 1 > shortlen) ## filter too short regions

    return(coverRegion)
  } else {
    return(filter(regionf, FALSE))
  }
}


##' @inheritParams filterSeg_
##' @importFrom magrittr %<>% %>%
##' @importFrom dplyr filter bind_cols select mutate
##' @rdname filterutility
##' @keywords internal
##'
filterSegInter_ <- function(regionf, rateMat, shortlen) {

  minf <- min(regionf)
  maxf <- max(regionf)

  ## intersect =====~~~~~ or ~~~~==== (only these tow cases)
  ## no need to reduce
  rateMat %<>% filter(tRate < 1 & tRate > 0)
  ## no tRate < 1 is filtered
  if (nrow(rateMat) == 0) {
    return(filter(regionf, FALSE))
  } else {
    interRegion <- rateMat %>% ## select intersect regions
      mutate(start = if_else(minf < maxstart, minf, minend + 1L)) %>%
      mutate(end = if_else(maxf > minend, maxf, maxstart - 1L)) %>% ## start/end
      select(start:end) %>%
      filter(end - start + 1 > shortlen)
  }

  return(interRegion)

}


##' @inheritParams FilterBlacklist
##' @param corerow A row of the CNV in a \code{tbl_df} form.
##' @importFrom magrittr %<>% %>%
##' @importFrom dplyr select filter
##' @rdname filterutility
##' @keywords internal
##'
filterRow_ <- function(corerow, blacklist, overlaprate, shortlen) {

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
    return(filter(corerow, FALSE))
  }
  else if (all(frate == 0)) {
    ## case 3: 0 overlap regions
    return(corerow)
  }
  else {
    ## case 4: has > 0 & < overlaprate overlap regions
    seg <- filterSeg_(select(corerow, start, end),
                      rateMat,
                      overlaprate = overlaprate,
                      shortlen = shortlen,
                      corerow$chromosome,
                      corerow$type,
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
##' @importFrom dplyr filter mutate do group_by ungroup
##' @importFrom tibble tibble
##' @references \href{http://penncnv.openbioinformatics.org/en/latest/user-guide/annotation/#filtering-cnv-calls-by-user-specified-criteria}{cytoband extend}
##' @export
##'
bl_cytoband <- function(cyto, extend = 5e5L) {

  cytobl <- cyto %>%
    filter(color %in% c('acen', 'gvar', 'stalk')) %>% ## filter
    mutate(start = if_else(start > extend, start - extend, 0L)) %>%
    mutate(end = end + extend) %>% ## extend
    group_by(chromosome) %>% ## sort reduce sort by chromosome
    do(SortRegion(tibble(start = .$start, end = .$end))) %>%
    do(ReduceRegion(tibble(start = .$start, end = .$end), gap = 0L)) %>%
    do(SortRegion(tibble(start = .$start, end = .$end))) %>%
    ungroup

  return(cytobl)
}


## nator <- system.file('extdata', 'example.cnvnator', package = 'CNVanno') %>%
##   read_cnvnator %>%
##   filter_cnvnator %>%
##   Segment(gap = 10L)


## FilterBlacklist(nator, bl_cytoband(hg19cyto), overlaprate = 0.5, shortlen = 1000L, gap = 0L, n = 2) %>%
##   FilterBlacklist(hg19bl, overlaprate = 0.5, shortlen = 1000L, gap = 100000L, n = 2) %>%
##   slot('coreCNV') %>%
##   write.csv('tmp1.csv')


## FilterBlacklist(nator, hg19bl, overlaprate = 0.5, shortlen = 1000L, gap = 0L, n = 2) %>%
##   FilterBlacklist(bl_cytoband(hg19cyto), overlaprate = 0.5, shortlen = 1000L, gap = 100000L, n = 2) %>%
##   slot('coreCNV') %>%
##   write.csv('tmp2.csv')


## CNVanno:::filterRow_(nator@coreCNV[92, ], hg19bl, 0.5, 1000L)[2, ] %>% CNVanno:::filterRow_(bl_cytoband(hg19cyto), 0.5, 1000L)
## CNVanno:::filterRow_(nator@coreCNV[92, ], bl_cytoband(hg19cyto), 0.5, 1000L)

## kit <- system.file('extdata', 'example.cnvkit', package = 'CNVanno') %>%
##   read_cnvkit %>%
##   filter_cnvkit %>%
##   Segment(gap = 10L)


## FilterBlacklist(kit, bl_cytoband(hg19cyto), overlaprate = 0.5, shortlen = 1000L, gap = 0L, n = 2) %>%
##   FilterBlacklist(hg19bl, overlaprate = 0.5, shortlen = 1000L, gap = 1000000L, n = 2)

## FilterBlacklist(kit, hg19bl, overlaprate = 0.5, shortlen = 1000L, gap = 0L, n = 2) %>%
##   FilterBlacklist(bl_cytoband(hg19cyto), overlaprate = 0.5, shortlen = 1000L, gap = 1000000L, n = 2)
