##' @include AllClasses.R
NULL

##' Read in CNV file
##'
##' \itemize{
##'   \item \code{read_cnvkit()}: Read in CNV generated from \code{CNVkit} (version 0.9.3) as a \code{data.frame}.
##'   \item \code{read_cnvnator()}: Read in CNV generated from \code{CNVnator} (version 0.3.3) as a \code{data.frame}.
##' }
##'
##' @title Standard read in CNV files
##' @param cnvpath The path of cnv file.
##' @return
##' \itemize{
##'   \item \code{read_cnvkit()}/\code{read_cnvkit}: A \code{data.frame}.
##' }
##'
##' @examples
##' require('magrittr')
##'
##' kit <- system.file('extdata', 'example.cnvkit', package = 'CNVanno') %>% read_cnvkit
##'
##' nator <- system.file('extdata', 'example.cnvnator', package = 'CNVanno') %>% read.cnvnator
##'
##' \dontrun{
##' ## write CNV file
##' write.cnvkit(cnsMat, 'cnsMat.txt')
##' }
##'
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @importFrom readr read_tsv
##' @importFrom dplyr select mutate case_when
##' @rdname rwcnv
##' @export
##'
read_cnvkit <- function(cnvpath) {
  cnvin <- read_tsv(cnvpath)

  cnv <- cnvin %>%
    mutate(type = case_when(cn > 2 ~ 'duplication',
                            cn < 2 ~ 'deletion',
                            TRUE ~ 'normal')) %>%
    select(chromosome:end, type)

  params <- cnvin %>%
    select(-(chromosome:end))

  cnvkit <- new('RawCNV',
                cnv,
                params = params,
                method = 'CNVkit')

  return(cnvkit)
}


##' @inheritParams read_cnvkit
##' @importFrom magrittr  %>%
##' @importFrom stringr str_extract str_sub
##' @importFrom readr read_tsv
##' @importFrom dplyr
##' @rdname rwcnv
##' @export
##'
read_cnvnator <- function(cnvpath) {
  cnvin <- read_tsv(cnvpath,
                    col_names = c('type', 'location', 'length', 'normalized_RD', 'e_val1', 'e_val2', 'e_val3', 'e_val4', 'q0')) %>% select(-length)


  ## separate cnv
  cnv  <- cnvin %>%
    select(type:location) %>%
    mutate(chromosome = location %>%
             as.character() %>%
             str_extract('^.*:') %>%
             str_sub(1, -2)) %>%
    mutate(start = location %>%
             as.character() %>%
             str_extract(':\\d+') %>%
             str_sub(2, -1) %>%
             as.integer()) %>%
    mutate(end = location %>%
             as.character() %>%
             str_extract('-\\d+') %>%
             str_sub(2, -1) %>%
             as.integer())%>%
    select(chromosome:end, type)

  params  <- cnvin %>%
    select(-(type:location))

  cnvnator <- new('RawCNV',
                cnv,
                params = params,
                method = 'CNVnator')

  return(cnvnator)
}

