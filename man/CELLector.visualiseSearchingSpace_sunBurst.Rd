\name{CELLector.visualiseSearchingSpace_sunBurst}
\alias{CELLector.visualiseSearchingSpace_sunBurst}
\title{
Visualising the CELLector searching space as a sunburst
}
\description{
  This function visualise the CELLector searching space as an interactive sunburst.
}
\usage{
CELLector.visualiseSearchingSpace_sunBurst(searchSpace)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{searchSpace}{A CELLector searching space as returned by the
  \code{CELLector.Build_Search_Space} function
}
}
\author{
Hanna Najgebauer and Francesco Iorio
}
\seealso{
\code{\link{CELLector.Build_Search_Space}},
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

### visualising the Searching space as sunburst
CELLector.visualiseSearchingSpace_sunBurst(searchSpace = CSS)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{visualisation}% use one of  RShowDoc("KEYWORDS")

