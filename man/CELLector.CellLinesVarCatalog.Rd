\name{CELLector.CellLinesVarCatalog}
\alias{CELLector.CellLinesVarCatalog}
\docType{data}
\title{
Built-in Cell Line variants' catalogue
}
\description{
List of somatic variants found in cell-lines, from Iorio et al, Cell 2016 [1]
}
\usage{data("CELLector.CellLinesVarCatalog")}
\format{
  A data frame with 486243 observations (one for each somatic variants) on the following 13 variables.
  \describe{
    \item{\code{SAMPLE}}{a character vector specifying the name of the cell line in which the variant has been found}
    \item{\code{COSMIC_ID}}{a numeric vector specifying COSMIC [2] identifier of the cell line in which the variant has been found}
    \item{\code{Cancer.Type}}{a character vector specifying the TCGA [3] classification of the tissue of origin of the cell line in which the variant has been
    found}
    \item{\code{Gene}}{a character vector specifying the HUGO [4] symbol of the gene hosting the variant under consideration (from Ensemble v56)}
    \item{\code{Transcript}}{a character vector specifying the identifier of the transcript affected by the variant under consideration (from Ensemble v56)}
    \item{\code{cDNA}}{a character vector specifying variant position and nucleotide change relating to the cDNA}
    \item{\code{AA}}{a character vector specifying aminoacid positon and alteration related to the variant under consideration}
    \item{\code{Classification}}{a character vector specifying a summary of the variant type	(e.g. missense, frameshift, nonsense, etc)}
    \item{\code{Gene_list}}{a character vector: when non empty then the variant is within a high-confidence cancer driver gene (according to [1]), thus resulting in a cancer functional event (CFE) summarised in the Cell line Binary Event Matrices (BEMs)}
    \item{\code{Recurrence.Filter}}{a character vector: if 'Yes' then the variant under consideration is observed in COSMIC [2] (v68) at a minimal threshold frequency (specified in [1])}
    \item{\code{Subs}}{a numeric vector: Missense/substitution variants occurring in codons mutated  in the systematic screen data in COSMIC [2] (v68) (select >=3}
    \item{\code{Truncating}}{a numeric vector: Truncating variant count from the systematic screen data in COSMIC [2] (v68) (select >10)}
    \item{\code{inframe}}{a numeric vector: Inframe indel alterations occurring in codons mutated in the systematic screen data in COSMIC [2] (v68) (select >=3)}
  }
}
\references{
[1] Iorio, F. et al. A Landscape of Pharmacogenomic Interactions in Cancer. Cell 166, 740–754 (2016).

[2] Forbes, S. A. et al. COSMIC: exploring the world’s knowledge of somatic mutations in human cancer. Nucleic Acids Res. 43, D805–11 (2015).

[3] Hutter C and Zenklusen JC. The Cancer Genome Atlas: Creating Lasting Value beyond Its Data. Cell. 2018;173(2):283–285. doi:10.1016/j.cell.2018.03.042

[4] Braschi, B. et al. Genenames.org: the HGNC and VGNC resources in 2019. Nucleic Acids Res. Epub 2018 Oct 10. PMID: 30304474 DOI: 10.1093/nar/gky930
}
\examples{
data(CELLector.CellLinesVarCatalog)
head(CELLector.CellLinesVarCatalog)
}
\author{
Francesco Iorio (fi9323@gmail.com)
}
\seealso{
\code{\link{CELLector.CellLine.BEMs}}, \code{\link{CELLector.CellLine.BEMs_v2}}, \code{\link{CELLector.CFEs}}
}
\keyword{datasets}
