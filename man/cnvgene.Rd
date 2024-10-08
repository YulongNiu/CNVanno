% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/anno_gene.R
\name{AnnoCNVGeneCore}
\alias{AnnoCNVGeneCore}
\alias{AnnoGeneOverlap_}
\title{The CNV gene database annotation}
\usage{
AnnoCNVGeneCore(corerow, annodb)

AnnoGeneOverlap_(corerow, annodb)
}
\arguments{
\item{corerow}{A row of the CNV in a \code{tbl_df} form.}

\item{annodb}{A \code{tbl_df} indicating a single annotation database. In \code{AnnoSVCore()}, it should at least contains "chromosome", "start", "end", "type" (gain or loss), and "clinical_significance" columns. In \code{AnnoCNVPopuCore()}, it should at least contains "chromosome", "start", "end", "type" (gain or loss), may contain "gain_frequency"and "loss_frequency". In \code{AnnoCNVGeneCore()}, it should at least contains "chromosome", "start" and "end" columns.}

\item{annoGene}{A \code{data.frame} of the OMIM, ClinGen_HaploTriplo, ExAC_plI, and DECIPHER_Haplo databases.}
}
\value{
A \code{list}

A \code{tbl_df} summary annotated genes/regions.
}
\description{
Annotation core function for the OMIM, ClinGen_TriHaplo, ExAC_plI, and DECIPHER_Haplo
}
\details{
Annotation of single CNV to one gene database
}
\examples{
data(CNVdb)
data(kit)

## ClinGen_TriHaploGRCh37
AnnoCNVGeneCore(kit@coreCNV[3, ], CNVdb$ClinGen_TriHaploGRCh37)

## ExAC_pLI
AnnoCNVGeneCore(kit@coreCNV[3, ], CNVdb$ExAC_pLI)

## DECIPHER_Haplo
AnnoCNVGeneCore(kit@coreCNV[3, ], CNVdb$DECIPHER_Haplo)

##OMIM
AnnoCNVGeneCore(kit@coreCNV[1, ], CNVdb$OMIMGRCh38)

}
\author{
Yulong Niu \email{yulong.niu@hotmail.com}
}
\keyword{internal}
