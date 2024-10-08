% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/anno_template.R
\name{AnnoCNVGeneRefGene2OMIM}
\alias{AnnoCNVGeneRefGene2OMIM}
\title{Annotation OMIM (GRCh38) from the RefGene (GRCh37)}
\usage{
AnnoCNVGeneRefGene2OMIM(refgeneTable, OMIMdb)
}
\arguments{
\item{refgeneTable}{A \code{data.frame} of annotated RefGene table.}

\item{OMIMdb}{The OMIM database.}
}
\value{
A \code{list} of annotated OMIM.
}
\description{
Annotation CNV with the OMIM (GRCh38) database from the RefGene (GRCh37)
}
\details{
Step1: Normal annotation CNV with the RefGene (GRCh37) database.

Step2: Use Entrez Gene ID as keys to combine RefGene and OMIM database.
}
\examples{
data(CNVdb)
data(kit)

refGene <- AnnoCNVBatch(kit, AnnoCNVGeneCore, CNVdb$RefGeneGRCh37, n = 2)

## cross OMIM database
crossOMIM <- AnnoCNVGeneRefGene2OMIM(refGene[[1]], CNVdb$OMIMGRCh38)
}
\author{
Yulong Niu \email{yulong.niu@hotmail.com}
}
