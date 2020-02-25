\name{CELLector.hms_look_up}
\alias{CELLector.hms_look_up}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Decoding identifers of hypermethylated gene promoters
}
\description{
This functions shows the annotation for coded hypermethylated gene promoters as defined in [1].
}
\usage{
CELLector.hms_look_up(hms_ID, hmsId_decode, TCGALabel)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hms_ID}{
A string containin the identifier of the hypermethylated gene promoters. Full list available in the \code{\link{CELLector.CFEsV2}} object.
}
  \item{hmsId_decode}{
A data frame containing the hypermethylated gene promoter annotations,
available in the \code{\link{CELLector.CFEs.HMSid_decode}} object
}
  \item{TCGALabel}{
A TCGA label indicating the cancer type under consideration:
\emph{BLCA, BRCA, COREAD, GBM , HNSC, KIRC, LAML, LGG, LUAD, LUSC, OV, PRAD, SKCM, STAD, THCA, UCEC}
available in this version.
}
}
\value{
A data frame with a single line containing the hypermethylated gene promoter indicated in input.
}
\references{
[1] Iorio, F. et al. A Landscape of Pharmacogenomic Interactions in Cancer. Cell 166, 740–754 (2016).
}
\author{
Francesco Iorio (fi9323@gmail.com)
}
\note{
%%  ~~further notes~~
}
\seealso{
\code{\link{CELLector.CFEsV2}}
}
\examples{
data(CELLector.CFEs.HMSid_decode)
CELLector.hms_look_up('hms26',
                      CELLector.CFEs.HMSid_decode,
                      TCGALabel = 'COREAD')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{annotation/decoding}% use one of  RShowDoc("KEYWORDS")
