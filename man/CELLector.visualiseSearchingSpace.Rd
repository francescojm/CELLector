\name{CELLector.visualiseSearchingSpace}
\alias{CELLector.visualiseSearchingSpace}
\title{
Visualising the CELLector searching space as a collapsibe tree
}
\description{
  This function visualise the CELLector searching space as an interactive collapsible tree.
  Moving the mouse over the nodes of this tree shows patients' subtype details and (optionally) matched cell lines
}
\usage{
CELLector.visualiseSearchingSpace(searchSpace, CLdata = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{searchSpace}{A CELLector searching space as returned by the
  \code{CELLector.Build_Search_Space} function

}
  \item{CLdata}{
  If different from \code{NULL} (default) this argument must contain a data frame in which the first two columns contain the COSMIC [1] identiefiers and names of cell lines (one per row), respectively, and then binary entries indicating the status of cancer functional events (CFEs, as defined in [2], one per column) across cell lines. The format is the same of the entries of the list in the built-in \code{CELLector.CellLine.BEMs}. Based on this data the cell lines are mapped onto the CELLector searching space and the number of cell lines mapped in each subtype showed when moving the mouse on the corresponding node.
}
}
\references{
[1] Forbes, S. A. et al. COSMIC: exploring the world’s knowledge of somatic mutations in human cancer. Nucleic Acids Res. 43, D805–11 (2015).

[2] Iorio, F. et al. A Landscape of Pharmacogenomic Interactions in Cancer. Cell 166, 740–754 (2016).
}
\author{
Hanna Najgebauer and Francesco Iorio
}
\seealso{
\code{\link{CELLector.Build_Search_Space}},

\code{\link{CELLector.CellLine.BEMs}}
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{visualisation}% use one of  RShowDoc("KEYWORDS")
