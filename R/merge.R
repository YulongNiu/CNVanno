##' @include AllClasses.R AllGenerics.R region.R
NULL



##' Merge two \code{CoreCNV}.
##'
##' Merge CNVs with reciprocal overlap regions.
##' @title Merge \code{CoreCNV}
##' @inheritParams Merge
##' @return A \code{CoreCNV} object.
##' @examples
##' require('magrittr')
##' data(nator)
##' data(kit)
##'
##' Merge(list(nator, kit), reciprate = 0.5, n = 2)
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @importFrom magrittr %>%
##' @importFrom dplyr do group_by ungroup select everything
##' @importFrom tibble tibble
##' @importFrom methods new
##' @rdname Merge-methods
##' @exportMethod Merge
##'
setMethod(f = 'Merge',
          signature = c(corelist = 'list', reciprate = 'numeric'),
          definition = function(corelist, reciprate, n, ...) {

            res <- corelist %>%
              mergePrepare_ %>%
              group_by(chromosome, type) %>%
              do(mergeSearch_(tibble(start = .$start, end = .$end, method = .$method), reciprate = reciprate, n = n)) %>%
              ungroup %>%
              select(chromosome, start:end, everything()) %>%
              new('CoreCNV', coreCNV = .)

            return(res)
            })



##' Merge internal functions.
##'
##' \itemize{
##'   \item \code{mergePrepare_()}: Prepare the input \code{CoreCNV} objects.
##'   \item \code{mergeHasOverlap_()}: Has passed overlap regions.
##'   \item \code{mergeReduce_()}: Reduce selected regions.
##' }
##'
##' @title Internal functions for merge
##' @inheritParams Merge
##' @return
##' \itemize{
##'   \item \code{mergePrepare_()}: A \code{tbl_df} merged CNV.
##'   \item \code{mergeHasOverlap_()}: A \code{logical vector}.
##'   \item \code{mergeReduce_()}: A code{tbl_df} with one row.
##' }
##' A \code{tbl_df} of input CNVs.
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @importFrom magrittr %>%
##' @importFrom dplyr group_by do select bind_rows everything ungroup
##' @rdname mergeutility
##' @keywords internal
##'
mergePrepare_ <- function(corelist) {

  core <- corelist %>%
    lapply(function(x) {return(x@coreCNV)}) %>%
    bind_rows ## row bind

  return(core)
}



##' @inheritParams Merge
##' @param regf, regt A \code{tbl_df} including "start", "end", and "method" columns.
##' @importFrom doParallel registerDoParallel stopImplicitCluster
##' @importFrom foreach foreach %dopar%
##' @importFrom iterators iter
##' @rdname mergeutility
##' @keywords internal
##'
mergeHasOverlap_ <- function(regf, regt, reciprate, n) {

  ## Example
  ## regf <- tibble(start = c(1L, 5L, 7L), end = c(8L, 10L, 9L), method = c('nator', 'kit', 'nator'))
  ## regt <- tibble(start = c(2L, 6L, 8L, 17L, 19L), end = c(12L, 9L, 13L, 20L, 31L), method = sample(c('kit', 'nator'), 5, replace = TRUE))
  ## mergeSeg_(regf, regt, 0.5, n = 2)

  registerDoParallel(cores = n)
  itx <- iter(regf, by = 'row')

  resLogic <- foreach(i = itx, .combine = `|`) %dopar% {
    eachRate <- OverlapRegionRate(i, regt)
    return(eachRate$fRate > reciprate & eachRate$tRate > reciprate)
  }

  ## stop multiple cores
  stopImplicitCluster()

  return(resLogic)
}


##' @param reglist A \code{list} contains \code{tbl_df}.
##' @importFrom dplyr mutate bind_rows
##' @importFrom magrittr %<>% %>%
##' @rdname mergeutility
##' @keywords internal
##'
mergeReduce_ <- function(reglist) {

  reg <- lapply(reglist, function(x) {

    m <- x$method %>%
      unique %>%
      paste(collapse = ';')

    x %<>% ## sort reduce
      SortRegion %>%
      ReduceRegion(gap = 0L) %>%
      mutate(method = m)

    return(x)
  }) %>%
  bind_rows

  return(reg)
}


##' @inheritParams Merge
##' @inheritParams mergeHasOverlap_
##' @importFrom dplyr slice bind_rows
##' @importFrom magrittr %<>%
##' @rdname mergeutility
##' @keywords internal
##'
mergeSearch_ <- function(regt, reciprate, n) {

  ## Example
  ## regf <- tibble(start = c(1L, 5L, 7L), end = c(8L, 10L, 9L))
  ## regt <- tibble(start = c(2L, 6L, 8L, 17L, 19L), end = c(12L, 9L, 13L, 20L, 31L))
  ## reg <- bind_rows(regf, regt) %>% SortRegion %>% mutate(method = sample(c('kit', 'nator'), 8, replace = TRUE))
  ## mergeSearch_(reg, 0.5, n = 2)
  res <- list()
  regf <- slice(regt, 1)
  regt <- slice(regt, -1)
  i <- 1
  res[[i]] <- regf

  while(TRUE) {

    ## check last row
    if (nrow(regt) == 0) {
      break
    } else {}

    eachLogic <- mergeHasOverlap_(regf = regf, regt = regt, reciprate = reciprate, n = n)
    if (sum(eachLogic) == 0) {
      ## no passed region
      regf <- slice(regt, 1)
      regt <- slice(regt, -1)
      i <- i + 1
      res[[i]] <- regf
    } else {
      regf <- slice(regt, which(eachLogic))
      regt <- slice(regt, which(!eachLogic))
      res[[i]] %<>% bind_rows(regf)
    }
  }

  ## combine and sort
  res %<>%
    mergeReduce_ %>%
    SortRegion

  return(res)
}


## ## test instansv

## library('CNVanno')
## library('magrittr')
## library('dplyr')

## tmp1 <- kit@coreCNV %>% mutate(size = end - start + 1L) %>% mutate(type = case_when(type == 'gain' ~ 'dup', type == 'loss' ~ 'del')) %>% rename(pos1 = start, pos2 = end, methods = method) %>% select(chromosome:pos2, size, everything())

## tmp2 <- kit@coreCNV %>% mutate(size = end - start + 1L) %>% mutate(type = case_when(type == 'gain' ~ 'dup', type == 'loss' ~ 'del')) %>% rename(pos1 = start, pos2 = end, methods = method) %>% select(chromosome:pos2, size, everything())

## library('intansv')
## methodsMerge(others = bind_rows(tmp1, tmp2), overLapPerDup = 0.5, overLapPerDel=0.5)

