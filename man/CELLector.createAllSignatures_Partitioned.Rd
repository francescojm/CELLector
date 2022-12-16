\name{CELLector.createAllSignatures_Partitioned}
\alias{CELLector.createAllSignatures_Partitioned}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extract signatures underlying the cancer patients' subtypes in the partitioned CELLector search space
}
\description{
This function takes in input the partitioned CELLector search space encoded as a \emph{navigable table} in a data frame. The signature rule is already derived to defined the paritioned space and extracted here in the same format of \code{\link{CELLector.createAllSignatures}}
}
\usage{
CELLector.createAllSignatures_Partitioned(NavTab)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{NavTab}{
    A CELLector searching space partitioned encoded as a navigable table in data frame format, as returned by the \code{CELLector.Build_Search_Space_Partitioned} function.
}
}
\value{
A list with two vectors of strings and two numerical vector. Each element of the first two vectors represent a signature of cancer functional events (CFEs, defined in [1]) corresponding to a group in the CELLector searching space partitioned. This is expressed as a logic formula (rule), which a cancer patient's genome must satisfy in order to be included in the group under consideration. The first vector (\code{S}) contains decoded signatures, i.e. where the CFEs involving copy number alterations are represented by a genomic loci and contained cancer driver genes. The second vector (\code{ES}) contains signatures of CFEs as they are represented in the binary event matrix containing the patients genomic data used to build the CELLector searching space. Further deatils are provided in [2]. The third vector (\code{STS}) contains the percentage of cancer patients belonging to the group represented by the signatures. The fourth vector (\code{SC}) contains the number of cancer patients from the cohort belonging to the group.
}
\references{
[1] Iorio, F. et al. A Landscape of Pharmacogenomic Interactions in Cancer. Cell 166, 740–754 (2016).

[2] Najgebauer, H. et al. Genomics Guided Selection of Cancer in vitro Models.

https://doi.org/10.1016/j.cels.2020.04.007
}
\author{
Lucia Trastulla and Francesco Iorio
}
\seealso{
\code{\link{CELLector.Build_Search_Space_Partitioned}}
\code{\link{CELLector.createAllSignatures}}
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
CSS_p <- CELLector.Build_Search_Space_Partitioned(ctumours = t(tumours_BEM),
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
Signatures <- CELLector.createAllSignatures_Partitioned(CSS_p$partitioned)
data.frame(Signatures = Signatures$S,'SubType Size'=Signatures$STS, 'SubType Count'=Signatures$SC)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{annotation/decoding}
