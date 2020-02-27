\name{CELLector.Tumours_buildBEM}
\alias{CELLector.Tumours_buildBEM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Building a Genomic Binary Event Matrix (BEM) for primary tumours
}
\description{
This function takes in input a catalogue of somatic genomic variants observed in primary tumours (or it uses a built in catalogue from TCGA, presented in [1]) and it converts it into a presence/absence (binary) matrix, which can be processed by the CELLector package and CELLector shiny app for identifying patient subtypes, and map in vitro models onto these.
}
\usage{
CELLector.Tumours_buildBEM(varCat = NULL,
                           Cancer_Type,
                           GenesToConsider = NULL,
                           VariantsToConsider = NULL)
}
\arguments{
  \item{varCat}{
    A data frame containing a catalogue of somatic genomic variants observed in primary tumours, with one row per variant. The format should be the same of the \code{CELLector.PrimTumVarCatalog} data object (which is used when this parameter is set to its default \code{NULL} value) or at least contain the following column headers \code{SAMPLE}, \code{Cancer.Type}, \code{Gene}, \code{cDNA}, \code{AA}, \code{Recurrence.Filter}

}
  \item{Cancer_Type}{
  A string specifying the cancer type for which individual variants should be extracted from the catalogue and assembled into the final matrix. It must be a value included in the \code{Cancer.Type} column of the \item{varCat} data object
}
  \item{GenesToConsider}{
  A list of strings with HGNC symbols [2] for genes hosting the variants to be extracted from the catalogue and assembled into the final matrix. When set to its default \code{NULL} value, all genes hosting at least one variants are considered.
}
  \item{VariantsToConsider}{
  A list of individual somatic variants to be extracted from the catalogue and assembled into the final matrix. The format should be the same of the \code{CELLector.RecfiltVariants} (which is used when this parameter is set to its default \code{NULL} value)
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
}
\references{
[1] Iorio, F. et al. A Landscape of Pharmacogenomic Interactions in Cancer. Cell 166, 740–754 (2016).

[2] Braschi, B. et al. Genenames.org: the HGNC and VGNC resources in 2019. Nucleic Acids Res. Epub 2018 Oct 10. PMID: 30304474 DOI: 10.1093/nar/gky930

[3] Tate JG, Bamford S, Jubb HC, et al. COSMIC: the Catalogue Of Somatic Mutations In Cancer. Nucleic Acids Res. 2019;47(D1):D941–D947. doi:10.1093/nar/gky1015
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{CELLector.PrimTumVarCatalog}}
}
\examples{
## loading high-confidence cancer driver genes from [1]
data(CELLector.HCCancerDrivers)

## loading COSMIC [3] variants observed it at least two patients from [1]
data(CELLector.RecfiltVariants)

## Assembling a BRCA primary tumour binary event matrix (BEM)
BRCA_tum_BEM<-
  CELLector.Tumours_buildBEM(Cancer_Type = 'BRCA',
                             VariantsToConsider =
                             CELLector.RecfiltVariants)

## showing first 100 entries of the BEM
BRCA_tum_BEM[1:10,1:10]

## showing a bar diagram with mutation frequency of 30 top frequently altered genes
barplot(100*sort(rowSums(BRCA_tum_BEM),
                 decreasing=TRUE)[1:30]/ncol(BRCA_tum_BEM),
                 las=2,ylab='\% patients')

}
\keyword{analysis}
\keyword{annotation/decoding}
