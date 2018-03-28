\name{CELLector.CFEs.CNAid_decode}
\alias{CELLector.CFEs.CNAid_decode}
\docType{data}
\title{
Decoding table for copy number alteration cancer functional events
}
\description{
A table with identifiers of cancer functional events (CFEs) involving chromosomal regions of recurrent copy number alterations (RACSs, as defined by [1], i.e. identified throgh ADMIRE [2]) and their annotation.
}
\usage{data("CELLector.CFEs.CNAid_decode")}
\format{
  A data frame with 731 observations (one for each CNA CFE) on the following 15 variables.
  \describe{
    \item{\code{Identifier}}{The RACS identifer, as defined in [1]}
    \item{\code{CancerType}}{A TCGA label indicating the cancer type where the RACS has been identified (via ADMIRE [2])}
    \item{\code{Recurrent}}{A string specifying whether the RACS under consideration is frequently amplified (value = \code{Amplification}) or deleted ) (value = \code{deleted})}
    \item{\code{chr}}{Chromosome number of the RACS}
    \item{\code{start}}{Starting position of the RACS}
    \item{\code{stop}}{Ending position of the RACS}
    \item{\code{nGenes}}{Number of protein coding genes included in the RACS}
    \item{\code{locus}}{Genomic locus of the RACS}
    \item{\code{ContainedGenes}}{A string with comma separated symbols of the genes included in the RACS}
    \item{\code{CNA_Identifier}}{A string containing the identifer of the RACS as it appears in the Binary Event Matrix (BEM) of the cancer type specified in the \code{CancerType} field included in the

    \code{CELLector.CellLine.BEMs} and the \code{CELLector.PrimTum.BEMs} data objects}
  }
}
\details{
This data frame contains a comprehensive annotation of the CFEs involving RACSs appearing in the BEMs of cell lines and primary tumours, contained in the \code{CELLector.CellLine.BEMs} and the \code{CELLector.PrimTum.BEMs} data objects.
Please note that the same RACS identifier across multiple cancer types might indicate different chromosomal regions, therefore in order to be decode it should be considered jointly with the TCGA label of the data frame it has been extracted from.

This table is used by the \code{\link{CELLector.cna_look_up}} function to decode the identifier of CFE involving a RACS.
}
\references{
[1] Iorio, F. et al. A Landscape of Pharmacogenomic Interactions in Cancer. Cell 166, 740–754 (2016).

[2] van Dyk, E., Reinders, M. J. T. & Wessels, L. F. A. A scale-space method for detecting recurrent DNA copy number changes with analytical false discovery rate control. Nucleic Acids Res. 41, e100 (2013).

}
\examples{
data(CELLector.CFEs.CNAid_decode)
haed(CELLector.CFEs.CNAid_decode)

data(CELLector.CellLine.BEMs)
colnames(CELLector.CellLine.BEMs$COREAD)[8]

CELLector.cna_look_up(cna_ID = colnames(CELLector.CellLine.BEMs$COREAD)[8],
                      TCGALabel = 'COREAD',
                      cnaId_decode = CELLector.CFEs.CNAid_decode)
}
\seealso{
\code{CELLector.CellLine.BEMs}, \code{CELLector.PrimTum.BEMs}, \code{CELLector.cna_look_up}

}
\keyword{datasets}
