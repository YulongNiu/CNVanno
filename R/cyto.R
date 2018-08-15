##' Find cytoband of given CNVs.
##'
##' Map the cytoband of given CNVs.
##'
##' @title Find cytoband
##' @importFrom doParallel registerDoParallel stopImplicitCluster
##' @importFrom foreach foreach %dopar%
##' @importFrom iterators iter
##' @return A \code{character} vector.
##' @examples
##' require('magrittr')
##' data(hg19cyto)
##'
##' kit <- system.file('extdata', 'exampleseg.cnvkit', package = 'CNVanno') %>% read_cnvkit %>% filter_cnvkit %>% Segment(interlen = 10L)
##'
##' cnsFilteredCyto <- cnsFile %>% read.cnvkit %>% filter.cnvkit %>% FindCyto(hg19cyto)
##'
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @rdname Cytoband-methods
##' @exportMethod Cytoband
##'
Cytoband <- function(cnvcore, cyto, n = 2) {

  registerDoParallel(cores = n)

  itx <- iter(cnvcore, by = 'row')

  cnvCyto <- foreach(i = itx, .combine = c) %dopar% {
    eachCyto <- cyto[cyto[, 1] %in% i[1], , drop = FALSE]
    eachLogic <- OverlapRegion(i[2:3], eachCyto[, 2:3, drop = FALSE])
    eachCyto <- paste(eachCyto[eachLogic, 4], collapse = ';')
    return(eachCyto)
  }

  ## stop multiple cores
  stopImplicitCluster()

  return(cnvCyto)
}


