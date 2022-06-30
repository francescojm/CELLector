\name{CELLector.Build_Search_Space_Partitioned}
\alias{CELLector.Build_Search_Space_Partitioned}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
CELLector partitioned search space derivation
}
\description{
This function assembles a user defined CELLector search space analysing genomic data from a larg cohort of cancer patients (specified in input). It identifies recurrent subtypes with matched genomic signatures (as combination of cancer functional events (CFEs), defined in [1]), linking them into a hierarchical structure shaped as a a binary three with a corresponding navigable table, as detailed in [2] and converting it into K non-overlapping groups defined as partitioned version (see details).
}
\usage{
CELLector.Build_Search_Space_Partitioned(ctumours,
                                         cancerType,
                                         minlen=1,
                                         verbose=TRUE,
                                         mutOnly=FALSE,
                                         cnaOnly=FALSE,
                                         includeHMS=FALSE,
                                         minGlobSupp=0.01,
                                         FeatureToExclude=NULL,
                                         pathway_CFEs = NULL,
                                         pathwayFocused=NULL,
                                         subCohortDefinition=NULL,
                                         NegativeDefinition=FALSE,
                                         cnaIdMap=NULL,
                                         cnaIdDecode=NULL,
                                         hmsIdDecode=NULL,
                                         cdg=NULL,
                                         UD_genomics=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ctumours}{
  A binary event matrix (BEM) modeling a cohort of cancer patients. With cancer functional events (CFEs) on the columns and sample identifers on the rows. See \code{\link{CELLector.PrimTum.BEMs}} for further details
}
  \item{cancerType}{
  The cancer type under consideration (specified via a TCGA label): currently available types = \emph{BLCA, BRCA, COREAD, GBM, HNSC, KIRC, LAML, LGG, LUAD, LUSC, OV, PRAD, SKCM, STAD, THCA, UCEC}
}
  \item{minlen}{
    The minimal length of the genomic signatures (how many indivudal CFEs it is made of) in order to be considered in the analysis (1 by default)
}
  \item{verbose}{
A boolean argument specifying whether step-by-step information on the algorithm progression should be displayed run-time
}
  \item{mutOnly}{
  A boolean argument specifying whether only CFEs involving somatic mutations should be considered in the analysis. If the \code{cnaOnly} argument is equal to \code{TRUE} then this must be \code{FALSE} (default value)
}
  \item{cnaOnly}{
  A boolean argument specifying whether only CFEs involving copy number alterations (CNAs) of chromosomal segments that are recurrently CN altered should be considered in the analysis. If the \code{mutOnly} argument is equal to \code{TRUE} then this must be \code{FALSE} (default value)
}
  \item{includeHMS}{
  A boolean argument specitying whether methylation data should be considered while building the searching space (\code{FALSE} by default).
}
  \item{minGlobSupp}{
Minimal size of the outpputted subtypes, as ratio of the number patients included in the whole cohort (1\% by default).
}
  \item{FeatureToExclude}{
A string (or a vector of strings) with identifiers of CFEs that should be ignored
}
  \item{pathway_CFEs}{
A named list of string vectors, whose elements are CFEs involving genes in a biological pathway (specified by the name of the corresponding entry). A list for 14 key cancer pathways is contained in the \code{\link{CELLector.Pathway_CFEs}} data object (see corresponding help page for further deatails)
}
  \item{pathwayFocused}{
  If different from \code{NULL} (default value), it should be a vector of strings. In this case the analysis will consider only CFEs involving genes in a set of pathways, whose names are contained in this argument and must be present as names of the \code{pathway_CFEs} argument
}
  \item{subCohortDefinition}{
  If different from \code{NULL} (default value), it should be a string containing the identifier of a CFE. In this case the analysis will consider only the primary tumour samples harbouring (or not harbouring, depending on the \code{NegativeDefinition} argument) the specified CFE
}
  \item{NegativeDefinition}{
  If the \code{subCohortDefinition} argument is not \code{NULL} then this paramenter determines whether to consider primary tumour samples that harbour (if equal to \code{FALSE}, default value) or not (if equal to \code{TRUE}) the specified CFE
}
  \item{cnaIdMap}{
A data frame mapping chromosomal regions of recurrent copy number amplifications/deletions in cancer (RACSs, as defined in [1]) identified via ADMIRE [3] in the context of specific cancer types to PanCancer RACSs. The built-in object \code{CELLector.CFEs.CNAid_mapping} (or an alternative data frame with the same format) should be used.
}
  \item{cnaIdDecode}{
A table with identifiers of cancer functional events (CFEs) involving chromosomal regions of recurrent copy number alterations (RACSs, as defined by [1], i.e. identified throgh ADMIRE [3]) and their annotation. The built-in object \code{CELLector.CFEs.CNAid_decode} (or an alternative data frame with the same format) should be used.
}
  \item{hmsIdDecode}{
  Data frame containing annotation for the hypermethylated gene promoters CFEs. The format should be the same of the \code{CELLector.CFEs.HMSid_decode} object.
}
  \item{UD_genomics}{
  A boolean argument specifying whether the analysis is performed on user defined genomic data (\code{TRUE}) or CELLector buit-in genomic data (\code{FALSE}, default value).
}
  \item{cdg}{
A list of genes that are used when decoding the identifiers of cancer functional events (CFEs) involving chromosomal regions of recurrent copy number alterations (RACSs, as defined by [1]). These will be visualised in the signatures containing the RACSs including them. A predefined list of high confidence cancer driver genes (from [1]) is provided as built-in data object (\code{CELLector.HCCancerDrivers})
}
}
\details{
The function builds cancer patient partition, creating \eqn{K} non-overlapping groups from the hierarchical division. It first constructs the CELLector space space (CSS) in the form of \emph{binary tree} and \emph{navigable table} for a assigned value of \code{minGlobSupp}. Each node of the binary tree starting from the root can include a \bold{left child}, i.e. the subset of samples that have the feature described in the parent node AND the feature described in the left child node and a \bold{right child}, i.e. the subset of samples that DO NOT have the feature described in the parent node AND have the feature described in the right child node (considered as complementary).

The partitioned structure is obtained from the hierarchical binary tree as follows. For each node express as a row in the \emph{navigable table}, let \eqn{U} be the set of samples in the considered node and \eqn{S} its signature.
If the node has a left child \eqn{U_l} (that includes feature \eqn{F_l}), the function first defines the set of samples \eqn{U_rm} to be removed from \eqn{U} as \eqn{U_rm} <- \eqn{U_l}. If \eqn{U_l} has a right child \eqn{U_r} (with feature \eqn{F_r}), then \eqn{U_rm} is updated as the union of \eqn{U_rm} and \eqn{U_r}. If \eqn{U_r} has in turn another right child, this procedure is repeated and \eqn{U_rm} is updated as descrived above until the considered node does not have a right child. Finally, the new set of samples is defined as \eqn{U_new = U \ U_rm} and the associated feature is \eqn{S_new = S, ~ F_l, ~ F_r, ...}.
Instead, if the node \eqn{U} does not have a left child, the iterative procedure move to the next node and \eqn{U} together with the associated signature \eqn{S} are kept as they are. Finally, a last node is creating as the remaining samples that are described by  any CELLector signature detected in the hierarchical version. The corresponding signature of this group is created from the negation of the root node and all the right children (recursively as before).
Note that the newly created groups can be composed of a fraction of samples lower than the predefined \code{minGlobSupp} value.

In this way, \eqn{K} non-overlapping groups are created with \eqn{K} equals to the number of signatures defined from the hierarchical version + 1.

}
\value{
A named list with the CELLector search space output for both hierarchical and partitioned output. The first is stored in \code{hierarchical} field and is the output of \code{\link{CELLector.Build_Search_Space}} function. The second is stored in \code{partitioned} and is \emph{navigable table} in the form of a data frame, output of \code{\link{CELLector.from_Hierarchical_to_Partition}} function. Each row represent a group of patients, with the columns indicating
\describe{
  \item{\code{Idx}}{A numerical index for the group}
  \item{\code{Signature}}{The combination of presence or absence of CFE, identified from the hierarchical strucutre as described in details}
  \item{\code{SignatureDecoded}}{Same as \code{Signature} but with identifiers of RACSs decoded, i.e. with loci and included driver genes (inputted in the \code{cdg} argument), indicated among brackets}
  \item{\code{Points}}{The identifiers of the patients in the group satisfying the signature rule}
  \item{\code{Total}}{Number of patients satisfying the signature rule}
  \item{\code{Support}}{Fraction of patients satisfying the signature rule compared to the total cohort}
  \item{\code{COLORS}}{A vector of strings containing hexadecimal color identifiers: one for each node. These are used by the visualisation functions (\code{CELLector.visualiseSearchingSpace}, and \code{CELLector.visualiseSearchingSpace_sunBurst}, and can be changed using the \code{CELLector.changeSScolors} function.
}
}}
\references{
[1] Iorio, F. et al. A Landscape of Pharmacogenomic Interactions in Cancer. Cell 166, 740–754 (2016).

[2] Najgebauer, H. et al. Genomics Guided Selection of Cancer in vitro Models.

https://doi.org/10.1016/j.cels.2020.04.007
}

\author{
Lucia Trastulla and Francesco Iorio
}

\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{CELLector.Build_Search_Space}}
\code{\link{CELLector.from_Hierarchical_to_Partition}}
}
\examples{

data(CELLector.PrimTum.BEMs_v2)

### Change the following two lines to work with a different cancer type
tumours_BEM<-CELLector.PrimTum.BEMs_v2$COREAD

### unicize the sample identifiers for the tumour data
tumours_BEM<-CELLector.unicizeSamples(tumours_BEM)

### loading decoding table for hypermethylation CFE identifiers
data(CELLector.CFEs.HMSid_decode)

### building a CELLector searching space
CSS_p <- CELLector.Build_Search_Space_Partitioned(ctumours = t(tumours_BEM),
                                  verbose = FALSE,
                                  minGlobSupp = 0.05,
                                  cancerType = 'COREAD',
                                  pathway_CFEs = CELLector.Pathway_CFEs,
                                  cnaIdMap = CELLector.CFEs.CNAid_mapping,
                                  cnaIdDecode = CELLector.CFEs.CNAid_decode,
                                  hmsIdDecode = CELLector.CFEs.HMSid_decode,
                                  cdg = CELLector.HCCancerDrivers,
                                  subCohortDefinition='TP53',
                                  NegativeDefinition=TRUE,
                                  includeHMS = TRUE)

### visualising partitioned patients and group-specific feature
CSS_p$partitioned

}
\keyword{analysis}
