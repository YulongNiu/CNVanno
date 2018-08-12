##' This class represents the raw data structure of CNV.
##'
##' @slot .Data An \code{tbl_df} from the \code{tibble} package. It contains four columns: 1st is the chromosome, 2nd is the CNV start position, 3rd is the CNV end position, and 4th is the CNV type ("deletion", "duplication", "normal")
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @importClassesFrom tibble tbl_df
##' @exportClass RawCNV
##'
setClass(Class = 'RawCNV',
         slots = c(params = 'tbl_df', method = 'character'),
         contains = 'tbl_df')

