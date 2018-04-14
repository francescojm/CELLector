\name{CELLector.makeSelection}
\alias{CELLector.makeSelection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Genomics guided selection of Cancer Cell lines
}
\description{
Given a CELLector searching space

(outputted by the \code{CELLector.Build_Search_Space} function) with tumour genomic subtypes and matched underlying signatures of cancer functional events (CFEs, as defined in [1]), and a map of human cancer cell lines on it (outputted by the \code{CELLector.buildModelMatrix} function), this function selects \emph{n} most representative cell lines by applying a greedy strategy described in [2] in order to maximise the covered genomic heterogeneity of primary tumours.
}
\usage{
CELLector.makeSelection(modelMat, n, searchSpace)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{modelMat}{
  A named binary matrix with tumour suptypes numerical identifiers on the rows, cell line names on the column and entries specifiyng whether the cell line in the column is representative of the subtype on the row (based on the collective presence/absence of the corresponding signature of CFEs). This is outputted by the \code{CELLector.buildModelMatrix} function starting from a CELLector search space (outputted by the \code{CELLector.Build_Search_Space}) and a cell line binary event matrix (BEM): a data frame in which the first two columns contain the COSMIC [3] identiefiers and names of cell lines (one per row), respectively, and then binary entries indicating the status of each CFEs (one per column) across cell lines. The format is the same of the entries of the list in the built-in \code{CELLector.CellLine.BEMs} object
}
  \item{n}{
An integer specifying the number of cell lines to select
}
  \item{searchSpace}{
  A CELLector searching space, outputted by the \code{CELLector.Build_Search_Space} from a BEM modeling a cohort of cancer patients. With cancer functional events (CFEs) on the columns and sample identifers on the rows. See CELLector.PrimTum.BEMs for further details
}
}

\value{
A data frame with one row per selected cell line and the following columns:
  \item{Tumour.SubType.Index}{The numerical index of the represented tumour subtype (this is the same index that the subtype has in the inputted CELLector searching space)}
  \item{Representative.Cell.Line}{The name of the selected cell line}
  \item{Signature}{The signature of CFEs underlying the subtype under consideration and collectively present in the selected cell line}
  \item{percentage.patients}{The size of the considered represented subtype with respect the whole cohort of cancer patients}

}
\references{
[1] Iorio, F. et al. A Landscape of Pharmacogenomic Interactions in Cancer. Cell 166, 740–754 (2016).

[2] Najgebauer, H. et al. Genomics Guided Selection of Cancer in vitro Models.

https://doi.org/10.1101/275032

[3] Forbes, S. A. et al. COSMIC: exploring the world’s knowledge of somatic mutations in human cancer. Nucleic Acids Res. 43, D805–11 (2015).
}
\author{
Hanna Najgebauer and Francesco Iorio
}
\seealso{
\code{\link{CELLector.Build_Search_Space}},

\code{\link{CELLector.buildModelMatrix}},

\code{\link{CELLector.CellLine.BEMs}},

\code{\link{CELLector.PrimTum.BEMs}}

}
\examples{
data(CELLector.PrimTum.BEMs)
data(CELLector.Pathway_CFEs)
data(CELLector.CFEs.CNAid_mapping)
data(CELLector.CFEs.CNAid_decode)
data(CELLector.HCCancerDrivers)
data(CELLector.CellLine.BEMs)


### Change the following two lines to work with a different cancer type
tumours_BEM<-CELLector.PrimTum.BEMs$COREAD
CELLlineData<-CELLector.CellLine.BEMs$COREAD

### unicize the sample identifiers for the tumour data
tumours_BEM<-CELLector.unicizeSamples(tumours_BEM)

### building a CELLector searching space focusing on three pathways
### and TP53 wild-type patients only
CSS<-CELLector.Build_Search_Space(ctumours = t(tumours_BEM),
                                  verbose = FALSE,
                                  minGlobSupp = 0.05,
                                  cancerType = 'COREAD',
                                  pathwayFocused = c("RAS-RAF-MEK-ERK / JNK signaling",
                                                     "PI3K-AKT-MTOR signaling",
                                                     "WNT signaling"),
                                  pathway_CFEs = CELLector.Pathway_CFEs,
                                  cnaIdMap = CELLector.CFEs.CNAid_mapping,
                                  cnaIdDecode = CELLector.CFEs.CNAid_decode,
                                  cdg = CELLector.HCCancerDrivers,
                                  subCohortDefinition='TP53',
                                  NegativeDefinition=TRUE)


### take all the signatures from the searching space
Signatures <- CELLector.createAllSignatures(CSS$navTable)

### mapping the cell lines on the CELLector searching space
ModelMat<-CELLector.buildModelMatrix(Signatures$ES,CELLlineData,CSS$navTable)


### selecting 10 cell lines
selectedCellLines<-CELLector.makeSelection(modelMat = ModelMat,
                                           n=10,
                                           searchSpace = CSS$navTable)

selectedCellLines

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{analysis}
