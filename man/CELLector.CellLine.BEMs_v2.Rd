\name{CELLector.CellLine.BEMs_v2}
\alias{CELLector.CellLine.BEMs_v2}
\docType{data}
\title{
Cell Lines' Binary Genomic/Epigenomic Event Matrices
}
\description{
A list containing 16 data frames (one for cancer type), identified through TCGA labels. Each of these data frames contains cell lines' \emph{binary event matrices} (BEMs) with the status (presence/absence) of \emph{cancer functional events} (CFEs, somatic mutations in high-confidence cancer driver genes, copy number altered genomic segments and hypermethylated gene promoters) as defined in [1].
}
\usage{data("CELLector.CellLine.BEMs_v2")}
\format{
A named list of data frames (with TCGA cancer type labels as names). Each of these data frames contains two columns with COSMIC [2] identiefiers and names of cell lines (one per row), respectively, and then binary entries indicating the status of each CFEs (one per column) across cell lines.
}
\details{
BEMs for cell lines from the Genomics of Drug Sensitivity in Cancer (GDSC1000, [1]) panel. Data is available for cell lines matching one among 16 different TCGA cancer types: \emph{BLCA, BRCA, COREAD, GBM, HNSC, KIRC, LAML, LGG, LUAD, LUSC, OV, PRAD, SKCM, STAD, THCA, UCEC}.

A decoding table for these labels is available at

\href{https://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources/Data_Summary.html}{https://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources/Data_Summary.html}.

Each data frame contains cell lines on the rows (with COSMIC identifiers and names, respectively on first and second column) and then a binary matrix with a CFE per column and entries indicating the presence/absence of a given CFE in a given cell line.

Gene symbols as row names indicate high confidence cancer driver genes and the entries in the corresponding rows indicate the presence/absence of somatic mutations. Row names with \emph{cna} as prefix indicates chromosomal segments that are recurrently copy number altered in cancer (RACSs, defined in [1]), \emp{hms} as a prefix indicates hypermethylated gene promoters (HMS, defined in [1]). A list with all the considered CFEs is available in the \code{\link{CELLector.CFEsV2}} data object. A decoding table for the RACSs is available in the \code{\link{CELLector.CFEs.CNAid_decode}} object, with the mapping realised by the values in the \code{CNA_identifier} column. A decoding table for the HMS is available in the \code{\link{CELLector.CFEs.HMSid_decode}} object, with the mapping realised by the values in the \code{hms_id} column.

Please note that the same RACS identifier across multiple cancer types might indicate different chromosomal regions, therefore in order to be decode it should be considered jointly with the TCGA label of the data frame it has been extracted from.

}
\references{
[1] Iorio, F. et al. A Landscape of Pharmacogenomic Interactions in Cancer. Cell 166, 740–754 (2016).

[2] Forbes, S. A. et al. COSMIC: exploring the world’s knowledge of somatic mutations in human cancer. Nucleic Acids Res. 43, D805–11 (2015).
}
\examples{
data(CELLector.CellLine.BEMs)
CELLector.CellLine.BEMs_v2$COREAD[1:10,
                               c("COSMIC_identifier","CellLine",'BRAF','hms70','cna27')]
}
\author{
Francesco Iorio (fi9232@gmail.com)
}
\seealso{
\code{\link{CELLector.CellLine.BEMs}},
\code{\link{CELLector.PrimTum.BEMs_v2}}, \code{\link{CELLector.CFEsV2}}, \code{\link{CELLector.CFEs.CNAid_decode}},
\code{\link{CELLector.CFEs.HMSid_decode}}
}
\keyword{datasets}
