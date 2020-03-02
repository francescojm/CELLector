\name{CELLector.RecfiltVariants}
\alias{CELLector.RecfiltVariants}
\docType{data}
\title{
Somatic mutations recurrently observed in cancer
}
\description{
Somatic variant passing a recurrence based threshold (defined in [1]) in COSMIC [2]
}
\usage{data("CELLector.RecfiltVariants")}
\format{
A vector of string with one entry per variant. Each string is a space delimited list of three elements:
The HUGO symbol of the gene hosting the variant, a string codifying the effect of the variant on cDNA, a string codifying the effect of the variant in terms of
aminoacid substitutions.
}
\references{
[1] Iorio, F. et al. A Landscape of Pharmacogenomic Interactions in Cancer. Cell 166, 740–754 (2016).

[2] Forbes, S. A. et al. COSMIC: exploring the world’s knowledge of somatic mutations in human cancer. Nucleic Acids Res. 43, D805–11 (2015).
}
\examples{
data(CELLector.RecfiltVariants)
head(CELLector.RecfiltVariants)
}
\keyword{datasets}
