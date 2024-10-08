% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/cyto.R
\docType{methods}
\name{Cytoband}
\alias{Cytoband}
\alias{Cytoband,CoreCNV,tbl_df-method}
\title{Find cytoband}
\usage{
Cytoband(core, cyto, n, ...)

\S4method{Cytoband}{CoreCNV,tbl_df}(core, cyto, n, ...)
}
\arguments{
\item{core}{A \code{CoreCNV} object.}

\item{cyto}{A \code{tbl_df} represents the cytoband.  A code{tbl_df} with at least three columns. 1st column is chromosome in the same format with \code{core}. 2nd and 3rd columns are start and end positions. 4th column is the cytoband. Other columns can be included but will not be used.}

\item{n}{The number of threads for parallel computation.}
}
\value{
A \code{character} vector.
}
\description{
Find cytoband of given CNVs.
}
\details{
Map the cytoband of given CNVs.
}
\examples{
require('magrittr')
data(hg19cyto)

kit <- system.file('extdata', 'exampleseg.cnvkit', package = 'CNVanno') \%>\%
  read_cnvkit \%>\%
  filter_cnvkit \%>\%
  Segment(gap = 10L)

kitcyto <- Cytoband(kit, hg19cyto, n = 2)

}
\author{
Yulong Niu \email{yulong.niu@hotmail.com}
}
\keyword{internal}
