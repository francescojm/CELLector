\name{CELLector.CFEs.HMSid_decode}
\alias{CELLector.CFEs.HMSid_decode}
\docType{data}
\title{
Decoding table for hypermethylated gene promoter cancer functional events
}
\description{
A table with identifiers of cancer functional events (CFEs) involving informative hypermehtylated gene promoters (hms) defined by [1].
}
\usage{data("CELLector.CFEs.HMSid_decode")}
\format{
  A data frame with 731 observations (one for each hms CFE) on the following 7 variables.
  \describe{
    \item{\code{Cancer.Types}}{A TCGA label indicating the cancer type where the hms has been identified as informative (as explained in [1])}
    \item{\code{Genomic.Coordinates}}{The genomic coordinates of the CgP island whose methylation signal was considered in [1] for this hms}
    \item{\code{Discretization.threshold}}{Threshold used to determine the status of this CgP island (hypermethylated or non hypermethylated) defined
    as the best separative ones between distributions of hypermethylated and non hypermethylated samples when looking at the methylation (beta) signal
    of the CgP island under consideration in cancer patients, whole procedure described in [1]}
    \item{\code{dip.statistic}}{Dip statistic resulting from testing the bimodality of the methylation (beta) signal of the CgP island under consideration in
    cancer patients. Whole procedure described in [1]}
    \item{\code{dip.pval}}{Dip p-value from testing the bimodality of the methylation (beta) signal of the CgP island under consideration in
    cancer patients. Whole procedure described in [1]}
    \item{\code{GN}}{Gene whose promoter is considered by the hms}
    \item{\code{hms_id}}{Alphanumerical identifiers of the hms}
  }
}
\details{
This data frame contains the annotation of the CFEs involving the hms appearing in the BEMs of cell lines and primary tumours, contained in the \code{CELLector.CellLine.BEMs_v2} and the \code{CELLector.PrimTum.BEMs_v2} data objects (version 2).

This table is used by the \code{\link{CELLector.hms_look_up}} function to decode the identifier of CFE involving an hms.
}
\references{
[1] Iorio, F. et al. A Landscape of Pharmacogenomic Interactions in Cancer. Cell 166, 740–754 (2016).
}
\examples{
data(CELLector.CFEs.HMSid_decode)
head(CELLector.CFEs.HMSid_decode)

data(CELLector.CellLine.BEMs_v2)
colnames(CELLector.CellLine.BEMs_v2$COREAD)[162]

CELLector.hms_look_up(hms_ID = colnames(CELLector.CellLine.BEMs_v2$COREAD)[162],
                      TCGALabel = 'COREAD',
                      hmsId_decode = CELLector.CFEs.HMSid_decode)
}
\author{
Francesco Iorio (fi9232@gmail.com)
}
\seealso{
\code{CELLector.CellLine.BEMs_v2}, \code{CELLector.PrimTum.BEMs_v2}, \code{CELLector.hms_look_up}
}
\keyword{datasets}
