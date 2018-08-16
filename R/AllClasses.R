##~~~~~~~~~~~~~~~~tbl_df and tbl class~~~~~~~~~~~
setOldClass(c('tbl_df', 'tbl', 'data.frame'))
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##' This class represents the raw data structure of CNV.
##'
##' @slot rawCNV A \code{tbl_df} from the \code{tibble} package. It contains four columns: 1st is the chromosome, 2nd is the CNV start position, 3rd is the CNV end position, and 4th is the CNV type ("gain", "loss", "normal")
##' @slot params A \code{tbl_df} from the \code{tibble} package. Each column is additional parameters.
##' @slot method A \code{character}
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @exportClass RawCNV
##'
setClass(Class = 'RawCNV',
         slots = c(rawCNV = 'tbl_df', params = 'tbl_df', method = 'character'))


##' This class represents the CNV structure used for further analysis.
##'
##' @slot .Data An \code{tbl_df} from the \code{tibble} package. It contains five columns: 1st is the chromosome, 2nd is the CNV start position, 3rd is the CNV end position, 4th is the CNV type ("gain", "loss", "normal"), 5th is the CNV calling method.
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @exportClass RawCNV
##'
setClass(Class = 'CoreCNV',
         slots = c(coreCNV = 'tbl_df'))

