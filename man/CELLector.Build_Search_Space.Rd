\name{CELLector.Build_Search_Space}
\alias{CELLector.Build_Search_Space}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
CELLector search space construction
}
\description{
This function assembles a user defined CELLector search space analysing genomic data from a larg cohort of cancer patients (specified in input). It identifies recurrent subtypes with matched genomic signatures (as combination of cancer functional events (CFEs), defined in [1]), linking them into a hierarchical structure shaped as a a binary three with a corresponding navigable table, as detailed in [2].
}
\usage{
CELLector.Build_Search_Space<-function(ctumours,
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
      !@#!@#!@#!@#!@#!@#!@#!@#!#@!
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
!@#!@#!@#!@#!@#!@#!@#!@#!#@!
}
  \item{UD_genomics}{
!@#!@#!@#!@#!@#!@#!@#!@#!#@!
}
  \item{cdg}{
A list of genes that are used when decoding the identifiers of cancer functional events (CFEs) involving chromosomal regions of recurrent copy number alterations (RACSs, as defined by [1]). These will be visualised in the signatures containing the RACSs including them. A predefined list of high confidence cancer driver genes (from [1]) is provided as built-in data object (\code{CELLector.HCCancerDrivers})
}
}
\details{
Starting from an initial cohort of patients affected by a given cancer type and modeled by the inputted binary event matrix (BEM), the most frequent alteration or set of molecular alterations (depending on the \code{minlen} argument) with the largest support (the subpopulation of patients in which these alterations occur simultaneously) is identified using the \code{eclat} function of the \code{arules} R package.

Based on this, the cohort of patients is split into two subpopulations depending on the collective presence or absence of the identified alterations. This process is then executed recursively on the two resulting subpopulations and it continues until all the alteration sets (with a support of minimal size, as specified in the \code{minGlobSupp} argument) are identified.

Each of the alterations sets identified through this recursive process is stored in a tree node. Linking nodes identified in adjacent recursions yields a binary tree: the CELLector search space. Each individual path (from the root to a node) of this tree defines a rule (signature), represented as a logic AND of multiple terms (or their negation), one per each node in the path. If the genome of a given patient in the analysed cohort satisfies the rule then it is contained in the subpopulation represented by the terminal node of that path. Collectively, all the paths in the search space provide a representation of the spectrum of combinations of molecular alterations observed in a given cancer type, and their clinical prevalence in the analysed patient population.
}
\value{
A named list with the CELLector search space stored as a \code{data.tree} object in the \code{TreeRoot} field
and as a \emph{navigable table}: a data frame with a row for each node of the tree and the following columns
\describe{
  \item{\code{Idx}}{A numerical index for the node}
  \item{\code{Item}}{The most supported CFE (or a combination of CFE), identified at the iteration in which the node has been added to the three, (i) in the whole cohort of patients (for the \code{Root}), (ii) in the sub population that satisfies the parent node rule (for \code{Left.Child} nodes) or (iii) its complement (for \code{Right.Child} nodes)}
  \item{\code{ItemsDecoded}}{Same as \code{Item} but with identifiers of RACSs decoded, i.e. with loci and included driver genes (inputted in the \code{cdg} argument), indicated among brackets}
  \item{\code{Type}}{The node type: Root (first node added), Right.Child (a node resulting from the analyses of the complementar population of patients with respect to that satisfifying the Parent node rule), Left.Child  (a node resulting from refining the population of patients satisfifying the Parent node rule)}
  \item{\code{Parent.Idx}}{The numerical index of the parent node (0 for the \code{Root})}
  \item{\code{AbsSupport}}{The number of patients satisfying the node rule}
  \item{\code{CurrentTotal}}{The number of patients included in the population under consideration at the iteration time of the node inclusion in the tree, this is the same of the parent's \code{AbsSupport} for \code{Left.Child} nodes}
  \item{\code{PercSupport}}{The ratio of patients collectively harbouring the combination of CFEs specified in \code{Items} within the subpopulation under consideration at the iteration time of the node inclusion in the tree (whose size is specified in \code{CurrentTotal})}
  \item{\code{GlobalSupport}}{The ratio of patients satisfying the node rule with respect to the total number of patients in the whole cohort}
  \item{\code{Left.Child.Index}}{Numerical index of the left child node (0 indicates absence of a left child node)}
  \item{\code{Right.Child.Index}}{Numerical index of the right child node (0 indicates absence of a right child node)}
  \item{\code{currentPoints}}{The identifiers of the patients in the sub-population under consideration at the iteration time of the node inclusion in the tree}
  \item{\code{currentFeatures}}{The CFEs considered at the at the iteration time of the node inclusion in the tree}
  \item{\code{positivePoints}}{The identifiers of the patients satisfying the node rule}
  \item{\code{COLORS}}{A vector of strings containing hexadecimal color identifiers: one for each node. These are used by the visualisation functions (\code{CELLector.visualiseSearchingSpace}, and \code{CELLector.visualiseSearchingSpace_sunBurst}}, and can be changed using the \code{CELLector.changeSScolors} function.
}

}
\references{
[1] Iorio, F. et al. A Landscape of Pharmacogenomic Interactions in Cancer. Cell 166, 740–754 (2016).

[2] Najgebauer, H. et al. Genomics Guided Selection of Cancer in vitro Models.

https://doi.org/10.1101/275032

[3] van Dyk, E., Reinders, M. J. T. & Wessels, L. F. A. A scale-space method for detecting recurrent DNA copy number changes with analytical false discovery rate control. Nucleic Acids Res. 41, e100 (2013).
}
\author{
Hanna Najgebauer and Francesco Iorio
}
\seealso{
\code{\link{CELLector.PrimTum.BEMs}},

\code{\link{CELLector.Pathway_CFEs}},

\code{\link{CELLector.CFEs.CNAid_mapping}},

\code{\link{CELLector.CFEs.CNAid_decode}},

\code{\link{CELLector.HCCancerDrivers}},

\code{\link{CELLector.visualiseSearchingSpace}},

\code{\link{CELLector.visualiseSearchingSpace_sunBurst}},

\code{\link{CELLector.changeSScolors}}
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

### visualising the CELLector searching space as a binary tree
CSS$TreeRoot

### visualising the first attributes of the tree nodes
CSS$navTable[,1:11]

### visualising the sub-cohort of patients whose genome satisfies the rule of the 4th node
str_split(CSS$navTable$positivePoints[4],',')

}
\keyword{analysis}% use one of  RShowDoc("KEYWORDS")
