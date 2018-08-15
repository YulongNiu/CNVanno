##' Find regions have intersections
##'
##' \itemize{
##'   \item \code{OverlapRegionRate()}: Check if extended regions have intersections with overlap rates.
##'   \item \code{OverlapRegion()}: Check regions have interactions
##'   \item \code{SortMat()}: Sort the 1st and 2nd columns resulting in 1st column is less than or equal to that of 2nd column.
##' }
##' @title Overlapped regions
##' @param regionf A \code{numeric} vector with length two, and the value of 1st position is smaller/equal to that of 2nd position.
##' @param regionMat A \code{tbl_df} with 2 columns. In the \code{OverlapRegionRate()} function, the rows should have overlap with the input \code{reginf}. The 1st column should be equal or smaller than the 2nd column, otherwise use the \code{SortMat()} sort the \code{regionMat} at first.
##' @return
##' \itemize{
##'   \item \code{OverlapRegionRate()}: A \code{logic} value.
##'   \item \code{OverlapRegion()}: A \code{numeric} matrix with two columns. 1st column is the overlap rate of `regionf`, 2nd column is the overlap rate of `regionMat`.
##'   \item \code{SortMat()}: A same object as the input \code{regionMat}
##' }
##' @examples
##' require('magrittr')
##'
##' tMat <- tibble(from = c(1, 103, 111, 8, 10), to = c(101, 112, 1000, 49, 86))
##' tReg <- c(100, 110)
##'
##' ## overlapped regions
##' inLogic <- OverlapRegion(tReg, tMat, extend = 0)
##'
##' ## overlap rate
##' OverlapRegionRate(tReg, tMat[inLogic, ])
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @rdname overlapregion
##' @export
##'
OverlapRegionRate <- function(regionf, regionMat) {

  ## overlap rate, 1st extract overlapped regions
  tLen <- regionMat[, 2] - regionMat[, 1]
  maxfLogic <- regionMat[, 1] < regionf[1]
  regionMat[maxfLogic, 1] <- regionf[1]
  mintLogic <- regionMat[, 2] > regionf[2]
  regionMat[mintLogic, 2] <- regionf[2]

  ## 2nd inter region
  interLen <- regionMat[, 2] - regionMat[, 1]
  fLen <- regionf[2] - regionf[1]

  overlapRate <- cbind(fRate = interLen/fLen,
                       tRate = interLen/tLen)

  return(overlapRate)
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
    mutate(from = if_else(from - extend > 0, from - extend, 0)) %>%
    mutate(to = to + extend)

  olLogic <- regionMat %>%
    transmute(from > max(regionf) | to < min(regionf)) %>%
    unlist %>%
    unname

  return(!olLogic)
}

##' @inheritParams OverlapRegionRate
##' @rdname overlapregion
##' @export
##'
SortMat <- function(regionMat) {

  ## sort regionMat
  sortLogic <- regionMat[, 1] > regionMat[, 2]
  if (sum(sortLogic) > 0) {
    regionMat[sortLogic, 1:2, drop = FALSE] <- regionMat[sortLogic, 2:1, drop = FALSE]} else {}

  return(regionMat)
}


