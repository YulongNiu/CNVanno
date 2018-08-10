##' Find cytoband of given cns.
##'
##' @title Find cytoband
##' @param cyto A standard cytoband file.
##' @param n The number of threads for parallel computation.
##' @param cnvcore A \code{data.frame} represent the CNV core data structure. 1st column is the chromosome like \code{chr1}, \code{chrX}, and \code{chrM}. 2nd and 3rd columns are start and end positions, respectively.
##' @importFrom doParallel registerDoParallel stopImplicitCluster
##' @importFrom foreach foreach %dopar%
##' @importFrom iterators iter
##' @return A \code{character} vector.
##' @examples
##' require('magrittr')
##' data(hg19cyto)
##'
##' cnsFile <- system.file('extdata', 'example.cnvkit', package = 'CNVanno')
##' cnsFilteredCyto <- cnsFile %>% read.cnvkit %>% filter.cnvkit %>% FindCyto(hg19cyto)
##'
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @export
##'
FindCyto <- function(cnvcore, cyto, n = 2) {

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


