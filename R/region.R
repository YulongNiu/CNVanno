##' Find regions have intersections
##'
##' \itemize{
##'   \item \code{OverlapRegionRate()}: Check if extended regions have intersections with overlap rates.
##'   \item \code{OverlapRegion()}: Check regions have interactions
##'   \item \code{SortRegion()}: Sort the 1st and 2nd columns resulting in 1st column is less than or equal to that of 2nd column.
##' }
##' @title Overlapped regions
##' @param regionf A \code{numeric} vector with length two, and the value of 1st position is smaller/equal to that of 2nd position.
##' @param regionMat A \code{tbl_df} with 2 columns. In the \code{OverlapRegionRate()} function. The 1st column should be equal or smaller than the 2nd column, otherwise use the \code{SortRegion()} sort the \code{regionMat} at first.
##' @return
##' \itemize{
##'   \item \code{OverlapRegionRate()}: A \code{logic} value.
##'   \item \code{OverlapRegion()}: A \code{numeric} matrix with four columns. 1st column is the overlap rate of `regionf`, 2nd column is the overlap rate of `regionMat`, 3rd and 4th columns are intersection start and end regions.
##'   \item \code{SortRegion()}: A same object as the input \code{regionMat}
##' }
##' @examples
##' require('magrittr')
##'
##' tMat <- tibble(start = c(1, 103, 111, 49, 10), end = c(101, 112, 1000, 8, 86)) %>% SortRegion
##' tReg <- c(100, 110)
##'
##' ## overlapped regions
##' OverlapRegion(tReg, tMat, extend = 0)
##'
##' ## overlap rate
##' OverlapRegionRate(tReg, tMat)
##'
##' ## sort matrix at first
##'
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @importFrom magrittr %>% %<>%
##' @importFrom dplyr mutate if_else select
##' @rdname overlapregion
##' @export
##'
OverlapRegionRate <- function(regionf, regionMat) {

  ## step1: calculate inter length
  regionMat %<>%
    mutate(fLen = max(regionf) - min(regionf)) %>%
    mutate(tLen = end - start) %>%
    mutate(maxstart = if_else(start > min(regionf), start, min(regionf))) %>%
    mutate(minend = if_else(end < max(regionf), end, max(regionf))) %>%
    mutate(interLen = minend - maxstart) %>%
    mutate(interLen = if_else(interLen < 0, 0, interLen))

  ## step2: calculate rate
  regionMat %<>%
    mutate(fRate = interLen / fLen) %>%
    mutate(tRate = interLen / tLen) %>%
    select(fRate, tRate, maxstart:minend)

  return(regionMat)
}


##' @inheritParams OverlapRegionRate
##' @param extend A \code{numeric} value indicates the extended length in both direction.
##' @importFrom dplyr mutate transmute
##' @importFrom magrittr %<>% %>%
##' @rdname overlapregion
##' @export
##'
OverlapRegion <- function(regionf, regionMat, extend = 100) {

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
##' @importFrom dplyr transmute if_else bind_cols
##' @export
##'
SortRegion <- function(regionMat) {

  start <- regionMat %>%
    transmute(start = if_else(start < end, start, end))

  end <- regionMat %>%
    transmute(end = if_else(end > start, end, start))

  reg <- bind_cols(start, end)

  return(reg)
}
