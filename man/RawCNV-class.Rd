% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllClasses.R
\docType{class}
\name{RawCNV-class}
\alias{RawCNV-class}
\title{This class represents the raw data structure of CNV.}
\description{
This class represents the raw data structure of CNV.
}
\section{Slots}{

\describe{
\item{\code{rawCNV}}{A \code{tbl_df} from the \code{tibble} package. It contains four columns: 1st is the chromosome, 2nd is the CNV start position, 3rd is the CNV end position, and 4th is the CNV type ("gain", "loss", "normal")}

\item{\code{params}}{A \code{tbl_df} from the \code{tibble} package. Each column is additional parameters.}

\item{\code{method}}{A \code{character}}
}}

\author{
Yulong Niu \email{yulong.niu@hotmail.com}
}
