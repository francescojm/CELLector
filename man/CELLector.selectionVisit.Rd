\name{CELLector.selectionVisit}
\alias{CELLector.selectionVisit}
\title{
Visiting the CELLector searching space
}
\description{
This function visits the CELLector searching space using the greedy algorithm described in [1] and it returns the cancer patients subtypes identifiers sorted for cell line selection, maximising the genomic heterogeneity observed in the analysed cohort of cancer patients when considering the first \emph{n} subtypes.
}
\usage{
CELLector.selectionVisit(TAV)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{TAV}{
  A CELLector searching space stored as a \emph{Navigable Table}, outputted by the \code{CELLector.Build_Search_Space}
  function
}
}
\value{
  A numerical vector containing the subtypes identifers sorted as explained above
}
\references{
[1] Najgebauer, H. et al. Genomics Guided Selection of Cancer in vitro Models.

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


tumours_BEM<-CELLector.PrimTum.BEMs$COREAD
CELLlineData<-CELLector.CellLine.BEMs$COREAD

### unicize the sample identifiers for the tumour data
tumours_BEM<-CELLector.unicizeSamples(tumours_BEM)

### building the CELLector searching space
CSS<-CELLector.Build_Search_Space(ctumours = t(tumours_BEM),
                                  verbose = FALSE,
                                  minGlobSupp = 0.05,
                                  cancerType = 'COREAD',
                                  pathway_CFEs = CELLector.Pathway_CFEs,
                                  cnaIdMap = CELLector.CFEs.CNAid_mapping,
                                  mutOnly = FALSE,
                                  cnaIdDecode = CELLector.CFEs.CNAid_decode,
                                  cdg = CELLector.HCCancerDrivers)

CELLector.selectionVisit(CSS$navTable)
}

\keyword{analysis}
