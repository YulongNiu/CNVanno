##' Find regions have intersections
##'
##' \itemize{
##'   \item \code{OverlapRegionRate()}: Check if extended regions have intersections with overlap rates.
##'   \item \code{OverlapRegion()}: Check regions have interactions
##'   \item \code{SortRegion()}: Sort the "start" and "end" columns resulting in the "start" column is less than or equal to that of the "end" column. At last, the "start" column is also sorted.
##'   \item \code{SortRegionChr()}: Sort the input which is split by the "chromosome" columns.
##'   \item \code{ReduceRegion()}: Merge concatenated regions.
##' }
##' @title Overlapped regions
##' @param regionf A \code{numeric} vector with length two, and the value of 1st position is smaller/equal to that of 2nd position.
##' @param regionMat A \code{tbl_df} with at least 2 columns named "start" and "end". In the \code{OverlapRegionRate()} and \code{ReduceRegion()} function, use the \code{SortRegion()} to sort the \code{regionMat} at first.
##' @return
##' \itemize{
##'   \item \code{OverlapRegionRate()}: A \code{logic} value.
##'   \item \code{OverlapRegion()}: A \code{numeric} matrix with four columns. 1st column is the overlap rate of `regionf`, 2nd column is the overlap rate of `regionMat`, 3rd and 4th columns are intersection start and end regions.
##'   \item \code{SortRegion()}: The same object as the input.
##'   \item \code{SortRegionChr()}: The same object as the input.
##'   \item \code{ReduceRegion()}: The same object as the input.
##' }
##' @examples
##' require('magrittr')
##'
##' tMat <- tibble(start = c(1L, 103L, 111L, 49L, 10L), end = c(101L, 112L, 1000L, 8L, 86L)) %>% SortRegion
##' tReg <- c(100L, 110L)
##'
##' ## overlapped regions
##' OverlapRegion(tReg, tMat, extend = 0L)
##'
##' ## overlap rate
##' OverlapRegionRate(tReg, tMat)
##'
##' ## reduce region
##' rMat <- tibble(start = c(1, 8, 14, 15, 19, 34, 40), end = c(12, 13, 19, 29, 24, 35, 46))
##' rMat %<>% SortRegion %>% arrange(start)
##' ReduceRegion(rMat, gap = 0L)
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @importFrom magrittr %>% %<>%
##' @importFrom dplyr mutate if_else select
##' @rdname overlapregion
##' @export
##'
OverlapRegionRate <- function(regionf, regionMat) {

  ## step1: calculate inter length
  regionMat %<>%
    mutate(fLen = max(regionf) - min(regionf) + 1) %>%
    mutate(tLen = end - start + 1) %>%
    mutate(maxstart = if_else(start > min(regionf), start, min(regionf))) %>%
    mutate(minend = if_else(end < max(regionf), end, max(regionf))) %>%
    mutate(interLen = minend - maxstart + 1) %>%
    mutate(interLen = if_else(interLen < 0, 0, interLen))

  ## step2: calculate rate
  regionMat %<>%
    mutate(fRate = interLen / fLen) %>%
    mutate(tRate = interLen / tLen) %>%
    select(fRate, tRate, maxstart:minend)

  return(regionMat)
}


##' @inheritParams OverlapRegionRate
##' @param extend A \code{integer} value indicates the extended length in both direction.
##' @importFrom dplyr mutate transmute
##' @importFrom magrittr %<>% %>%
##' @rdname overlapregion
##' @export
##'
OverlapRegion <- function(regionf, regionMat, extend = 100L) {

  regionMat %<>%
    mutate(start = if_else(start > extend, start - extend, 0)) %>%
    mutate(end = end + extend)

  olLogic <- regionMat %>%
    transmute(start > max(regionf) | end < min(regionf)) %>%
    unlist %>%
    unname

  return(!olLogic)
}



##' @inheritParams OverlapRegionRate
##' @rdname overlapregion
##' @importFrom magrittr %<>% %>%
##' @importFrom dplyr transmute if_else arrange
##' @export
##'
SortRegion <- function(regionMat) {

  startr <- regionMat %>%
    transmute(start = if_else(start < end, start, end)) %>%
    unlist

  endr <- regionMat %>%
    transmute(end = if_else(end > start, end, start)) %>%
    unlist

  ## replace start and end
  regionMat$start <- startr
  regionMat$end <- endr

  regionMat %<>% arrange(start)

  return(regionMat)
}


##' @param regionMatChr A \code{tbl_df} with at least columns named "chromosome", "start", and "end".
##' @importFrom magrittr %<>% %>%
##' @importFrom dplyr bind_rows
##' @export
##'
SortRegionChr <- function(regionMatChr) {

  regionMatChr %<>%
    split(., .$chromosome) %>%
    lapply(SortRegion) %>%
    bind_rows

  return(regionMatChr)
}


##' @inheritParams OverlapRegionRate
##' @inheritParams Segment
##' @rdname overlapregion
##' @importFrom magrittr %>%
##' @importFrom dplyr bind_cols
##' @export
##'
ReduceRegion <- function(regionMat, gap) {

  ## `regionMat` must be sorted by row and by column
  start <- regionMat$start
  end <- regionMat$end

  inter <- c(FALSE, start[-1] - end[-length(end)] - 1 <= gap)
  startIdx <- which(!inter)
  endIdx <- c(startIdx[-1] - 1,
              length(inter))

  res <- list(start = start[startIdx],
              end = end[endIdx]) %>%
    bind_cols

  return(res)
}
