\name{CELLector.Score}
\alias{CELLector.Score}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cell line quality score
}
\description{
  This function computes a quality score of each cell line in terms of its ability to represent an entire cohort of considered disease matching patients. This is quantified as a trade-off between two factors. The first factor is the length of the CELLector signatures (in terms of number of composing individual alterations) that are present in the cell line under consideration. This is proportional to the granularity of the representative ability of the cell line, i.e. the longest the signature the more precisely defined is the represented sub-cohort of patients. The second factor is the size of the patient subpopulation represented by the signatures that can be observed in the cell line under consideration, thus accounting for the prevalence of the sub-cohort modeled by that cell line. An input parameter allows for these two factor to be weighted equally or differently.
}
\usage{
CELLector.Score(NavTab, CELLlineData, alpha = 0.75)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{NavTab}{
A CELLector searching space encoded as binary tree in a navigable table, as returned by the \code{CELLector.Build_Search_Space} function.
}
  \item{CELLlineData}{
A data frame in which the first two columns contain the COSMIC [1] identiefiers (or any other alphanumeric identifier) and names of cell lines (one per row), respectively, and then binary entries indicating the status of cancer functional events (CFEs, as defined in [2], one per column) across cell lines. The format is the same of the entries of the list in the built-in \code{CELLector.CellLine.BEMs}. Based on this data the cell lines are mapped onto the CELLector searching space and their quality assessed.
}
  \item{alpha}{
A parameter that weights the contribution of two factors accounted to compute the quality score of a cell line: respectively, the length of a signature that can be observed in a given cell line (weight = alpha) and the ratio of the sub-cohort of patients represented by that signature (1 - alpha). A cell line can be positive for multiple signatures, thus resulting in multiple quality scores. The largest is selected for each cell line.
}
}
\value{
A data frame with one row per cell line and the following columns:
  \item{CellLines}{containing an alphanumerical identifier of each cell line. These are the same as in \code{CELLlineData};}
  \item{GlobalSupport}{the ratio of patients represented by the signature yielding the best
    \code{CELLectorScore} for the cell line under consideration;}
  \item{SignatureLength}{the length of the signature (in terms of number of individual CFEs [2]) yielding the best \code{CELLectorScore} for the cell line under consideration;}
  \item{CELLectorScores}{the CELLector quality scores;}
  \item{Signature}{The signature yielding the best CELLector quality score for the cell line under consideration.}
}
\references{
[1] Forbes, S. A. et al. COSMIC: exploring the world’s knowledge of somatic mutations in human cancer. Nucleic Acids Res. 43, D805–11 (2015).

[2] Iorio, F. et al. A Landscape of Pharmacogenomic Interactions in Cancer. Cell 166, 740–754 (2016).
}
\author{
Francesco Iorio (fi9232@gmail.com)
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

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

### Evaluating cell line quality
CScores<-CELLector.Score(NavTab=CSS$navTable,CELLlineData = CELLlineData)

### Visualising best cell lines
head(CScores)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{analysis}% use one of  RShowDoc("KEYWORDS")
