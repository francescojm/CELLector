\name{CELLector.createAllSignatures}
\alias{CELLector.createAllSignatures}
\title{
Derive signatures underlying the cancer patients' subtypes in the CELLector search space
}
\description{
This function takes in input the CELLector search space encoded as a binary tree in a \emph{navigable table}. Then, for each individual path (from the root to a node) of this tree it derives a rule (signature), represented as a logic AND of multiple terms (or their negation), one per each node in the path. Negations are added when right branches are encountered.
}
\usage{
CELLector.createAllSignatures(NavTab)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{NavTab}{
    A CELLector searching space encoded as binary tree in a navigable table, as returned by the \code{CELLector.Build_Search_Space} function.
}
}
\value{
  A list with two vectors of strings and a numerical vector. Each element of the first two vectors represent a signature of cancer functional events (CFEs, defined in [1]) corresponding to a node in the CELLector searching space. This is expressed as a logic formula (rule), which a cancer patient's genome must satisfy in order to be included in the sub-population represented by the node under consideration. The first vector (\code{S}) contains decoded signatures, i.e. where the CFEs involving copy number alterations are represented by a genomic loci and contained cancer driver genes. The second vector (\code{ES}) contains signatures of CFEs as they are represented in the binary event matrix containing the patients genomic data used to build the CELLector searching space. Further deatils are provided in [2]. The third vector (\code{STS}) contains the percentage of cancer patients belonging to the subtype represented by the signatures.
}
\references{
[1] Iorio, F. et al. A Landscape of Pharmacogenomic Interactions in Cancer. Cell 166, 740–754 (2016).

[2] Najgebauer, H. et al. Genomics Guided Selection of Cancer in vitro Models.

https://doi.org/10.1101/275032
}
\author{
Hanna Najgebauer and Francesco Iorio
}
\seealso{
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

### derive signatures from searching space
Signatures <- CELLector.createAllSignatures(CSS$navTable)

data.frame(Signatures = Signatures$S,'SubType Size'=Signatures$STS)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{annotation/decoding}
