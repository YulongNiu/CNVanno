##' @param raw a \code{RawCNV} object.
##' @param interlen An code{integer number} indicate the integer length between CNV (at same chromosome with same type).
##' @param ... Additional parameters.
##' @rdname Segment-methods
##' @keywords internal
##'
setGeneric(name = 'Segment',
           def = function(raw, interlen, ...){standardGeneric('Segment')})
