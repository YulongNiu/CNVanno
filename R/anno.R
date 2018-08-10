##' Find cytoband of given cns.
##'
##' @title Find cytoband
##' @param cyto A standard cytoband file.
##' @param n The number of threads for parallel computation.
##' @inheritParams write.cns
##' @inheritParams OverlapRegion
##' @importFrom doParallel registerDoParallel stopImplicitCluster
##' @importFrom foreach foreach %dopar%
##' @importFrom iterators iter
##' @return A \code{character} vector.
##' @examples
##' require('magrittr')
##' data(hg19cyto)
##'
##' cnsFile <- system.file('extdata', 'example.cnvkit', package = 'CNVanno')
##' cnsFilteredCyto <- cnsFile %>% read.cnvkit %>% FilterCNS %>% FindCyto(hg38cyto)
##'
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @export
##'
FindCyto <- function(cns, cyto, n = 2) {

  registerDoParallel(cores = n)

  itx <- iter(cns, by = 'row')

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






