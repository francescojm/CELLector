\name{CELLector.MSIstatus}
\alias{CELLector.MSIstatus}
\docType{data}
\title{
Cell lines' Microsatellite status
}
\description{
The microsatellite status of the cell lines in the CELLector collection, which can be stable (MSI-S), lowly instable (MSI-L), or highly instable (MSI-H) from [1]
}
\usage{data("CELLector.MSIstatus")}
\format{
A named vector of string with one entry per cell lines (with COSMIC [2] identifiers as names) specigying the MSI status of each cell line as detailed in the description above.
}
\references{
[1] Iorio, F. et al. A Landscape of Pharmacogenomic Interactions in Cancer. Cell 166, 740–754 (2016).

[2] Forbes, S. A. et al. COSMIC: exploring the world’s knowledge of somatic mutations in human cancer. Nucleic Acids Res. 43, D805–11 (2015
}
\author{
Francesco Iorio (fi9232@gmail.com)
}
\examples{
data(CELLector.MSIstatus)
head(CELLector.MSIstatus)
}
\keyword{datasets}
