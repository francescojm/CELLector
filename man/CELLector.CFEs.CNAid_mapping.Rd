\name{CELLector.CFEs.CNAid_mapping}
\alias{CELLector.CFEs.CNAid_mapping}
\docType{data}
\title{
Pan-Cancer/Cancer-Specific RACSs map.
}
\description{
A data frame mapping chromosomal regions of recurrent copy number amplifications/deletions in cancer (RACSs, as defined in [1]) identified via ADMIRE [2] in the context of specific cancer types to PanCancer RACSs.
}
\usage{data("CELLector.CFEs.CNAid_mapping")}
\format{
  A data frame with 425 observations (one for each PanCancer RACS) and a column for each of 27 different cancer types (specified by
  TCGA labels). The entry in position \emph{i,j} contains the identifier of the \emph{i}th PanCancer RACS in the context of the \emph{j}th
  cancer type (where available).
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
[1] Iorio, F. et al. A Landscape of Pharmacogenomic Interactions in Cancer. Cell 166, 740–754 (2016).

[2] van Dyk, E., Reinders, M. J. T. & Wessels, L. F. A. A scale-space method for detecting recurrent DNA copy number changes with analytical false discovery rate control. Nucleic Acids Res. 41, e100 (2013).

}
\examples{
data(CELLector.CFEs.CNAid_mapping)
head(CELLector.CFEs.CNAid_mapping)
}
\keyword{datasets}
