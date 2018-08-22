##' @param raw A \code{RawCNV} object.
##' @param gap An \code{integer} indicate the allowed length between two regions (at same chromosome with same type).
##' @param ... Additional parameters.
##' @rdname Segment-methods
##' @keywords internal
##'
setGeneric(name = 'Segment',
           def = function(raw, gap, ...){standardGeneric('Segment')})



##' @param core A \code{CoreCNV} object.
##' @param cyto A \code{tbl_df} represents the cytoband.  A code{tbl_df} with at least three columns. 1st column is chromosome in the same format with \code{core}. 2nd and 3rd columns are start and end positions. 4th column is the cytoband. Other columns can be included but will not be used.
##' @param n The number of threads for parallel computation.
##' @rdname Cytoband-methods
##' @keywords internal
##'
setGeneric(name = 'Cytoband',
           def = function(core, cyto, n, ...){standardGeneric('Cytoband')})


##' @param core A \code{CoreCNV} object.
##' @param blacklist A code{tbl_df} with at least three columns. 1st column is chromosome in the same format with \code{core}. 2nd and 3rd columns are start and end positions. Other columns can be included but will not be used. The blacklist must be reduced (use \code{ReduceRegionChr()}) and sorted (use \code{SortRegionChr()}).
##' @param overlaprate The threshold of overlap rate (overlaplen/CNVlen). The CNVs with smaller than or equal to this rate will be filtered.
##' @param shortlen A code{integer}. The CNVs with length shorter than or equal to this value will be filtered.
##' @inheritParams Cytoband
##' @inheritParams Segment
##' @param ... Additional parameters.
##' @rdname FilterBlacklist-methods
##' @keywords internal
##'
setGeneric(name = 'FilterBlacklist',
           def = function(core, blacklist, overlaprate, shortlen, gap, n, ...){standardGeneric('FilterBlacklist')})



##' @param corelist A list of \code{CoreCNV} objects used to merge.
##' @param reciprate The reciprocal overlap rate between two CNVs.
##' @param ... Additional parameters.
##' @inheritParams Cytoband
##' @rdname Merge-methods
##' @keywords internal
##'
setGeneric(name = 'Merge',
           def = function(corelist, reciprate, n, ...){standardGeneric('Merge')})

