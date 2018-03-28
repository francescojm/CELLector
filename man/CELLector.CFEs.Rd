\name{CELLector.CFEs}
\alias{CELLector.CFEs}
\docType{data}
\title{
Cancer Functional Events
}
\description{
Identifiers of cancer functional events (CFEs, i.e. somatic mutations in high confidence cancer driver genes or chromosomal regions of recurrent copy number amplification/deletion) from [1], which are also present in the binary event matrices of the cell lines and the primary tumours considered in this version of CELLector.
}
\usage{data("CELLector.CFEs")}
A vector of strings with one entry per identifier.
}
\details{
Gene symbols indicate somatic mutations in igh confidence cancer driver genes and entries with \emph{cna} prefix indicate chromosomal segments that are recurrently copy number altered in cancer (RACSs), both defined in [1].

A decoding table for the RACSs is available in the \code{\link{CELLector.CFEs.CNAid_decode}}, with the mapping realised by the values in the \code{CNA_identifier} column.

Please note that the same RACS identifier across multiple cancer types might indicate different chromosomal regions, therefore in order to be decode it should be considered jointly with the TCGA label of the data frame it has been extracted from.
}
\references{
[1] Iorio, F. et al. A Landscape of Pharmacogenomic Interactions in Cancer. Cell 166, 740–754 (2016).
}
\examples{
data(CELLector.CFEs)
head(CELLector.CFEs)
}
\seealso{
\code{\link{CELLector.PrimTum.BEMs}}, \code{\link{CELLector.CellLine.BEMs}}, \code{\link{CELLector.CFEs}}, \code{\link{CELLector.CFEs.CNAid_decode}}
}
\keyword{datasets}
