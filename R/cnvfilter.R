##' Filter the cns file generated from \code{CNVkit}
##'
##' \itemize{
##'   \item 1. filter the cns with cn number is 2 (no variance).
##'   \item 2. filter log2 value, and the default threshold is absolute 0.5.
##'   \item 3. filter low sequence depth.
##'   \item 4. filter X/Y chromosomes
##' }
##'
##' @title CNS filter
##' @param cngain A \code{numeric} value as the gain copy number threshold, and the default value is 2.
##' @param cnloss A \code{numeric} value as the loss copy number threshold, and the default value is 2.
##' @param log2value A \code{numeric} value as the absolute log2value threshold, and the default value is absolute 0.5.
##' @param depth A \code{numeric} value as the sequence depth, and the default value is 0.01
##' @param sexchrom A \code{logic} value whether filter sex chromosomes (X and Y), the default value is \code{TRUE}.
##' @inheritParams write.cns
##' @return A filtered cns format
##' @examples
##' require('magrittr')
##'
##' cnsFile <- system.file('extdata', 'example.cnvkit', package = 'CNVanno')
##' cnsFiltered <- cnsFile %>% read.cns %>% FilterCNS
##'
##' @author Yulong Niu \email{yulong.niu@@hotmail.com}
##' @export
##'
FilterCNS <- function(cns, cngain = 2, cnloss = 2, log2value = 0.5, depth = 0.01, sexchrom = TRUE) {

  ## filter gain and loss cnv
  cnsMat <- cns[(cns[, 'cn'] > cngain) | (cns[, 'cn'] < cnloss), , drop = FALSE]

  ## filter log2
  cnsMat <- cnsMat[abs(cnsMat[, 'log2']) >= log2value, , drop = FALSE]

  ## filter low depth
  cnsMat <- cnsMat[cnsMat[, 'depth'] > depth, , drop = FALSE]

  ## filter sex chromosomes
  if (sexchrom) {
    cnsMat <- cnsMat[(cnsMat[, 'chromosome'] != 'chrX') & (cnsMat[, 'chromosome'] != 'chrY'), , drop = FALSE]
  } else {}

  return(cnsMat)

}


