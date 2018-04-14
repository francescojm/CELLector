\name{CELLector.buildModelMatrix}
\alias{CELLector.buildModelMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Mapping cell lines on the CELLector searching space
}
\description{
This function maps cell line on the subtypes identified and assembled in the CELLector searching space, based on the collective presence/absence of the signatures of cancer functional events underlying these subtypes. The subtypes lacking representative cell lines are not considered and, in the output, the subtypes (indicated by their numerical id, which matches that in the CELLector searching space) are ranked based on the greedy algorithm described in [1] based on their covered genomic heterogeneity.
}
\usage{
CELLector.buildModelMatrix(Sigs, dataset, searchSpace)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Sigs}{
A vector of string, in which each element represents a signature of cancer functional events (CFEs, defined in [2]) corresponding to a node in the CELLector searching space. This is expressed as a logic formula (rule), which a cancer patient's genome must satisfy in order to be included in the sub-population represented by the node under consideration. This vector is outputted by the \code{CELLector.createAllSignatures} function starting from a CELLector searching space (created by the \code{CELLector.Build_Search_Space}) function
}
  \item{dataset}{
A data frame in which the first two columns contain the COSMIC [3] identiefiers and names of cell lines (one per row), respectively, and then binary entries indicating the status of each CFEs (one per column) across cell lines. The format is the same of the entries of the list in the built-in \code{CELLector.CellLine.BEMs} object.
}
  \item{searchSpace}{
A CELLector searching space encoded as binary tree in a navigable table, as returned by the \code{CELLector.Build_Search_Space} function
}
}
\value{
  A named binary matrix with suptypes numerical identifiers on the rows, cell line names on the column and entries specifiyng whether the cell line in the column is representative of the subtype on the row (based on the collective presence/absence of the corresponding signature of CFEs)
}
\references{
[1] Najgebauer, H. et al. Genomics Guided Selection of Cancer in vitro Models.

https://doi.org/10.1101/275032

[2] Iorio, F. et al. A Landscape of Pharmacogenomic Interactions in Cancer. Cell 166, 740–754 (2016).

[3] Forbes, S. A. et al. COSMIC: exploring the world’s knowledge of somatic mutations in human cancer. Nucleic Acids Res. 43, D805–11 (2015).
}
\author{
Hanna Najgebauer and Francesco Iorio
}
\seealso{
\code{\link{CELLector.createAllSignatures}},

\code{\link{CELLector.Build_Search_Space}},

\code{\link{CELLector.CellLine.BEMs}},

\code{\link{CELLector.Build_Search_Space}}
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

### mapping colorectal cancer cell lines onto the CELLector searching space
ModelMat<-CELLector.buildModelMatrix(Signatures$ES,CELLlineData,CSS$navTable)

head(ModelMat)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{analysis}% use one of  RShowDoc("KEYWORDS")
