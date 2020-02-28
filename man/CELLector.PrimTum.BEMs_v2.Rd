\name{CELLector.PrimTum.BEMs_v2}
\alias{CELLector.PrimTum.BEMs_v2}
\docType{data}
\title{
Primary Tumours' Binary Genomic/Epigenomic Event Matrices
}
\description{
A list containing 16 data frames (one for cancer type), identified through TCGA labels. Each of these data frames contains primary tumours' \emph{binary event matrices} (BEMs) with the status (presence/absence) of \emph{cancer functional events} (CFEs, somatic mutations in high-confidence cancer driver genes, copy number altered genomic segments and hypermethylated gene promoters) as defined in [1].
}
\usage{data("CELLector.PrimTum.BEMs_v2")}
\format{
 A named list of binary matrices (with TCGA cancer type labels as names). The entries of each of these matrices indicate the status (Present/Absent) of each CFE (one per row) across primary tumors samples (one per column).
}
\details{
BEMs of primary tumours from the Genomics of Drug Sensitivity in Cancer (GDSC1000, [1]) study. Data is available for 16 different TCGA cancer types: \emph{BLCA, BRCA, COREAD, GBM, HNSC, KIRC, LAML, LGG, LUAD, LUSC, OV, PRAD, SKCM, STAD, THCA, UCEC}.

A decoding table for these labels is available at

\href{https://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources/Data_Summary.html}{https://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources/Data_Summary.html}

Each data frame contains primary tumour samples on the columns and CFEs on the rows, with entries indicating the presence/absence of a given CFE in a given primary tumour sample.

Gene symbols as row names indicate high confidence cancer driver genes and the entries in the corresponding rows indicate the presence/absence of somatic mutations. Row names with \emph{cna} as prefix indicates chromosomal segments that are recurrently copy number altered in cancer (RACSs, defined in [1]), \emp{hms} as a prefix indicates hypermethylated gene promoters (HMS, defined in [1]). A list with all the considered CFEs is available in the \code{\link{CELLector.CFEsV2}} data object. A decoding table for the RACSs is available in the \code{\link{CELLector.CFEs.CNAid_decode}} object, with the mapping realised by the values in the \code{CNA_identifier} column. A decoding table for the HMS is available in the \code{\link{CELLector.CFEs.HMSid_decode}} object, with the mapping realised by the values in the \code{hms_id} column.

Please note that the same RACS identifier across multiple cancer types might indicate different chromosomal regions, therefore in order to be decode it should be considered jointly with the TCGA label of the data frame it has been extracted from.
}
\references{
[1] Iorio, F. et al. A Landscape of Pharmacogenomic Interactions in Cancer. Cell 166, 740–754 (2016).
}
\seealso{
\code{\link{CELLector.PrimTum.BEMs}}, \code{\link{CELLector.CellLine.BEMs_v2}}, \code{\link{CELLector.CFEsV2}}, \code{\link{CELLector.CFEs.CNAid_decode}}, \code{\link{CELLector.CFEs.HMSid_decode}}
}
\examples{
data(CELLector.PrimTum.BEMs_v2)
CELLector.PrimTum.BEMs_v2$COREAD[c('BRAF','KRAS','hms67','hms70'),]

}
\keyword{datasets}
