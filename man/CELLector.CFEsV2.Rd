\name{CELLector.CFEsV2}
\alias{CELLector.CFEsV2}
\docType{data}
\title{
Cancer Functional Events (version 2)
}
\description{
Identifiers of cancer functional events (CFEs, i.e. somatic mutations in high confidence cancer driver genes, chromosomal regions of recurrent copy number amplification/deletion, and informative hypermethylated gene promoters) from [1], which are also present in the binary event matrices (version 2) of the cell lines and the primary tumours considered in this version of CELLector.
}
\usage{data("CELLector.CFEsV2")}
\format{
A vector of strings with one entry per identifier.
}
\details{
Gene symbols indicate somatic mutations in high confidence cancer driver genes, entries with \emph{cna} prefix indicate chromosomal segments that are recurrently copy number altered in cancer (RACSs), entries with \emph(hms) prefix indicate hypermehylated gene promoters, all defined in [1].

A decoding table for the HMS is available in the \code{\link{CELLector.CFEs.HMSid_decode}} object, with the mapping realised by the values in the \code{hms_id} column.

Please note that the same RACS identifier across multiple cancer types might indicate different chromosomal regions, therefore in order to be decode it should be considered jointly with the TCGA label of the data frame it has been extracted from.
}
\references{
[1] Iorio, F. et al. A Landscape of Pharmacogenomic Interactions in Cancer. Cell 166, 740–754 (2016).
}
\examples{
data(CELLector.CFEsV2)
head(CELLector.CFEsV2)
}
\author{
Francesco Iorio (fi9232@gmail.com)
}
\seealso{
\code{\link{CELLector.PrimTum.BEMs_v2}}, \code{\link{CELLector.CellLine.BEMs_v2}},
\code{\link{CELLector.CFEs}}, \code{\link{CELLector.CFEs.CNAid_decode}}, \code{\link{CELLector.CFEs.HMSid_decode}}
}
\keyword{datasets}
