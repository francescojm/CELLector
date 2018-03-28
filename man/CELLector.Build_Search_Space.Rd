\name{CELLector.Build_Search_Space}
\alias{CELLector.Build_Search_Space}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
CELLector search space construction
}
\description{
This function assembles a user defined CELLector search space analysing genomic data from a larg cohort of cancer patients (spcified in input). It identifies recurrent subtypes with matched genomic signatures (as combination of cancer functional events (CFEs), defined in [1]), linking them into a hierarchical structure shaped as a a binary three with a corresponding navigable table, as detailed in [2].
}
\usage{
CELLector.Build_Search_Space(ctumours,
                             cancerType,
                             minlen = 1,
                             verbose = TRUE,
                             mutOnly = FALSE,
                             cnaOnly = FALSE,
                             minGlobSupp = 0.01,
                             FeatureToExclude = NULL,
                             pathway_CFEs = NULL,
                             pathwayFocused = NULL,
                             subCohortDefinition = NULL,
                             NegativeDefinition = FALSE,
                             cnaIdMap,
                             cnaIdDecode,
                             cdg)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ctumours}{
  A binary event matrix (BEM) modeling a cohort of cancer patients. With cancer functional events (CFEs) on the columns and sample identifers on the rows. See \code{\link{CELLector.PrimTum.BEMs}} for further details
}
  \item{cancerType}{
  The cancer type under consideration (spcified via a TCGA label): currently available types = \emph{BLCA, BRCA, COREAD, GBM, HNSC, KIRC, LAML, LGG, LUAD, LUSC, OV, PRAD, SKCM, STAD, THCA, UCEC}
}
  \item{minlen}{
    The minimal length of the genomic signatures (of how many indivudal CFE it needs to be composed) in order to be considered in the analysis (1 by default)
}
  \item{verbose}{
A boolean argument specifying whether step-by-step information on the algorithm progression should be displayed run-time
}
  \item{mutOnly}{
  A boolean argument specifying wether only CFE involving somatic mutations should be considered in the analysis. If the \code{cnaOnly} argument is equal to \code{TRUE} then this must be \code{FALSE} (default value)
}
  \item{cnaOnly}{
  A boolean argument specifying wether only CFE involving copy number alterations (CNAs) of chromosomal segments that are recurrently CN altered should be considered in the analysis. If the \code{mutOnly} argument is equal to \code{TRUE} then this must be \code{FALSE} (default value)
}
  \item{minGlobSupp}{
Minimal size of the outpputted subtypes, as ratio with respect to the whole cohort of patients
}
  \item{FeatureToExclude}{
A string (or a vector of strings) with identifiers of CFEs that should be ignored
}
  \item{pathway_CFEs}{
TO BE CONTINUED
}
  \item{pathwayFocused}{
%%     ~~Describe \code{pathwayFocused} here~~
}
  \item{subCohortDefinition}{
%%     ~~Describe \code{subCohortDefinition} here~~
}
  \item{NegativeDefinition}{
%%     ~~Describe \code{NegativeDefinition} here~~
}
  \item{cnaIdMap}{
%%     ~~Describe \code{cnaIdMap} here~~
}
  \item{cnaIdDecode}{
%%     ~~Describe \code{cnaIdDecode} here~~
}
  \item{cdg}{
%%     ~~Describe \code{cdg} here~~
}
}
\details{
Starting from an initial cohort of patients affected by a given cancer type and modeled by the inputted binary event matrix (BEM), the most frequent alteration or set of molecular alterations (depending on the \code{minlen} argument) with the largest support (the subpopulation of patients in which these alterations occur simultaneously) is identified using the \code{eclat} function of the \code{arules} R package.

Based on this, the cohort of patients is split into two subpopulations depending on the collective presence or absence of the identified alterations. This process is then executed recursively on the two resulting subpopulations and it continues until all the alteration sets (with a support of minimal size, as specified in the \code{minGlobSupp} argument) are identified.

Each of the alterations sets identified through this recursive process is stored in a tree node. Linking nodes identified in adjacent recursions yields a binary tree: the CELLector search space. Each individual path (from the root to a node) of this tree defines a rule (signature), represented as a logic AND of multiple terms (which can be also negated), one per each node in the path. If the genome of a given patient in the analysed cohort satisfies the rule then it is contained in the subpopulation represented by the terminal node of that path. Collectively, all the paths in the search space provide a representation of the spectrum of combinations of molecular alterations observed in a given cancer type, and their clinical prevalence in the analysed patient population.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{CELLector.PrimTum.BEMs}},

}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
