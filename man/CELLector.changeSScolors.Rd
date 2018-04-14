\name{CELLector.changeSScolors}
\alias{CELLector.changeSScolors}
\title{
Changing CELLector searching space color scheme
}
\description{
  This function reassigns a set of distinc random colors to the cancer patients' subtypes identified and stored in the
  CELLector searching space (compute by the \code{CELLector.Build_Search_Space} function). These colors are used by the visualisation functions \code{CELLector.visualiseSearchingSpace} and \code{CELLector.visualiseSearchingSpace_sunBurst}}
\usage{
CELLector.changeSScolors(searchSpace)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{searchSpace}{
A CELLector searching space as returned by the CELLector.Build_Search_Space function
}
}
\value{
A CELLector searching space as returned by the CELLector.Build_Search_Space function, with new colors
}
\author{
Hanna Najgebauer and Francesco Iorio
}
\seealso{
\code{\link{CELLector.Build_Search_Space}},

\code{\link{CELLector.Build_Search_Space}},

\code{\link{CELLector.visualiseSearchingSpace_sunBurst}}
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

### visualising the Searching space as collapsible tree
CELLector.visualiseSearchingSpace(searchSpace = CSS,CLdata = CELLlineData)
### visualising the Searching space as sunburst
CELLector.visualiseSearchingSpace_sunBurst(searchSpace = CSS)

### changing the searching space colors
CSS<-CELLector.changeSScolors(CSS)

### visualising the Searching space as collapsible tree
CELLector.visualiseSearchingSpace(searchSpace = CSS,CLdata = CELLlineData)
### visualising the Searching space as sunburst
CELLector.visualiseSearchingSpace_sunBurst(searchSpace = CSS)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{visualisation}% use one of  RShowDoc("KEYWORDS")
