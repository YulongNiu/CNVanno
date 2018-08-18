##' @param raw a \code{RawCNV} object.
##' @param gap An \code{integer} indicate the allowed length between two regions (at same chromosome with same type).
##' @param ... Additional parameters.
##' @rdname Segment-methods
##' @keywords internal
##'
setGeneric(name = 'Segment',
           def = function(raw, gap, ...){standardGeneric('Segment')})



##' @param core a \code{CoreCNV} object.
##' @param blacklist A code{tbl_df} with at least three columns. 1st column is chromosome in the same format with \code{core}. 2nd and 3rd columns are start and end positions. Other columns can be included but will not be used. The blacklist must be reduced (use \code{ReduceRegionChr()}) and sorted (use \code{SortRegionChr()}).
##' @param overlaprate The threshold of overlap rate (overlaplen/CNVlen).
##' @inheritParams Cytoband
##' @param ... Additional parameters.
##' @rdname FilterBlacklist-methods
##' @keywords internal
##'
setGeneric(name = 'FilterBlacklist',
           def = function(core, blacklist, overlaprate, n, ...){standardGeneric('FilterBlacklist')})



##' @param core A \code{CoreCNV} object.
##' @param cyto A \code{tbl_df} represents the cytoband.  A code{tbl_df} with at least three columns. 1st column is chromosome in the same format with \code{core}. 2nd and 3rd columns are start and end positions. 4th column is the cytoband. Other columns can be included but will not be used.
##' @param n The number of threads for parallel computation.
##' @rdname Cytoband-methods
##' @keywords internal
##'
setGeneric(name = 'Cytoband',
           def = function(core, cyto, n, ...){standardGeneric('Cytoband')})


