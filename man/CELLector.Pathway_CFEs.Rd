\name{CELLector.Pathway_CFEs}
\alias{CELLector.Pathway_CFEs}
\docType{data}
\title{
Cancer functional events in biological pathways
}
\description{
Lists of cancer functional events (CFEs) from [1] involving genes in 14 key cancer biological pathways
}
\usage{data("CELLector.Pathway_CFEs")}
\format{

Named list of string vectors, whose elements are CFEs involving genes in a fixed biological pathway.
  }
\references{
[1] Iorio, F. et al. A Landscape of Pharmacogenomic Interactions in Cancer. Cell 166, 740–754 (2016).
}
\examples{
data(CELLector.Pathway_CFEs)
CELLector.Pathway_CFEs$`RAS-RAF-MEK-ERK / JNK signaling`
}
\keyword{datasets}
