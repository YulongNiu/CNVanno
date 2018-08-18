##' @include AllClasses.R AllGenerics.R region.R
NULL


##' Find cytoband of given CNVs.
##'
##' Map the cytoband of given CNVs.
##'
##' @title Find cytoband
##' @return A \code{character} vector.
##' @examples
##' require('magrittr')
##' data(hg19cyto)
##'
##' kit <- system.file('extdata', 'exampleseg.cnvkit', package = 'CNVanno') %>%
##'   read_cnvkit %>%
##'   filter_cnvkit %>%
##'   Segment(gap = 10L)
##'
##' kitcyto <- Cytoband(kit, hg19cyto, n = 2)
##'
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @importFrom doParallel registerDoParallel stopImplicitCluster
##' @importFrom foreach foreach %dopar%
##' @importFrom iterators iter
##' @importFrom dplyr filter select
##' @importFrom magrittr %<>% %>%
##' @rdname Cytoband-methods
##' @exportMethod Cytoband
##'
setMethod(f = 'Cytoband',
          signature = c(core = 'CoreCNV', cyto = 'tbl_df'),
          definition = function(core, cyto, n, ...) {

            core  <- core@coreCNV

            registerDoParallel(cores = n)

            itx <- iter(core, by = 'row')

            cnvCyto <- foreach(i = itx, .combine = c) %dopar% {

              eachCyto <- filter(cyto, chromosome == i$chromosome)
              eachLogic <- eachCyto %>%
                OverlapRegion(select(i, start:end), ., 0L)
              eachCyto %<>%
                filter(eachLogic) %>%
                select(cytoband) %>%
                unlist %>%
                paste(collapse = ';')

              return(eachCyto)
            }

            ## stop multiple cores
            stopImplicitCluster()

            return(cnvCyto)
          })


