\name{CELLector.Summary_Projection}
\alias{CELLector.Summary_Projection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary of cell lines mapping on the partitioned CELLector search space
}
\description{
This function builds a summary table of the cell lines projection into groups defined from CELLector, both for hierarchical and partitioned versions. It builds from the derived binary matrix of signatures specific of each group obtained from \code{\link{CELLector.buildModelMatrix}} or \code{\link{CELLector.buildModelMatrix_Partitioned}}. The table includes a entry for each detected sub-population and the corresponding mapping cell lines that satisfy the same genomic signature rule, highlithing patient populations that lack in-vitro representations.
}
\usage{
CELLector.Summary_Projection(Signatures, ModelMat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Signatures}{
  Output of \code{\link{CELLector.createAllSignatures_Partitioned}} or \code{\link{CELLector.createAllSignatures}}, is a list including two vector of strings and one/two numeric verctor. Each entry represents a signature of cancer functional events corresponding to a sub-population in the CELLector searching space. Note that only the partitioned version includes the number of patients in each group that is also added to the final summary table.
}
  \item{ModelMat}{
A named binary matrix with sub-population/groups numerical identifiers on the rows, cell line names on the column and entries specifiyng whether the cell line in the column is representative of the subtype on the row (based on the collective presence/absence of the corresponding signature of CFEs). It is obtained as an output of \code{\link{CELLector.buildModelMatrix}} or  \code{\link{CELLector.buildModelMatrix_Partitioned}} functions
}
}

\value{
A table with each row representing a patients sub-population/group defined by a genomic signature, with the columns indicating
\describe{
  \item{\code{Subtype}}{A numerical index for the sub-population/group as in the Idx entry of the navigable table}
  \item{\code{Signature}}{The combination of presence or absence of CFE, identified from the hierarchical strucutre as described in details}
  \item{\code{SignatureDecoded}}{Same as \code{Signature} but with identifiers of RACSs decoded, i.e. with loci and included driver genes (inputted in the \code{cdg} argument), indicated among brackets}
  \item{\code{N_patients}}{Number of patients described by that genomic signature}
  \item{\code{P_patients}}{Percentage of patients described by that genomic signature over the considered cohort}
  \item{\code{N_CL}}{Number of cell lines mapping to that genomic signature}
  \item{\code{repr_CL}}{Cell lines names mapping to that genomic signature. "Lack of in vitro models" indicates that no cell lines satisfy that genomic signature}
   \item{\code{subpop_score}}{Score assigned to each sub-population/group computed as the product of percentage of patients defined by that genomic signature and the number of CFEs that are required to be present to define that signature.}
  }
}
\author{
Lucia Trastulla and Francessco Iorio
}
\note{
The sub-population score reflects the level of granularity of a group and the coverage defined in terms of percentage of patients. Higher scores indicates a group higly occuring and/or defined from the presence of multiple CFEs.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{CELLector.buildModelMatrix}},
\code{\link{CELLector.buildModelMatrix_Partitioned}},
\code{\link{CELLector.createAllSignatures}},
\code{\link{CELLector.createAllSignatures_Partitioned}}
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

## hierarchical ##
### take all the signatures from the searching space
Signatures <- CELLector.createAllSignatures(CSS_p$hierarchical$navTable)
### mapping colorectal cancer cell lines onto the CELLector searching space
ModelMat <- CELLector.buildModelMatrix(Signatures$ES,CELLlineData,CSS_p$hierarchical$navTable)

CELLector.Summary_Projection(Signatures, ModelMat)

## partitioned ##
### take all the signatures from the searching space
Signatures_p <- CELLector.createAllSignatures_Partitioned(CSS_p$partitioned)
### mapping colorectal cancer cell lines onto the CELLector searching space
ModelMat_p <- CELLector.buildModelMatrix_Partitioned(Signatures$ES,CELLlineData,CSS_p$partitioned)

CELLector.Summary_Projection(Signatures_p, ModelMat_p)


}
\keyword{analysis}
