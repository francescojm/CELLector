\name{CELLector.cna_look_up}
\alias{CELLector.cna_look_up}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Decoding identifers of chromosomal regions of recurrent Copy Number Alterations
}
\description{
This functions shows the annotation for a chromosomal region of recurrent copy number alterations (RACS) as defined in [1].
}
\usage{
CELLector.cna_look_up(cna_ID, cnaId_decode, TCGALabel)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cna_ID}{
A string containin the RACS identifier. Full list available in the \code{\link{CELLector.CFEs}} object.
}
  \item{cnaId_decode}{
A data frame containing the RACSs' annotation, available in the

\code{\link{CELLector.CFEs.CNAid_decode}} object
}
  \item{TCGALabel}{
A TCGA label indicating the cancer type under consideration: \emph{BLCA, BRCA, COREAD, GBM, HNSC, KIRC, LAML, LGG, LUAD, LUSC, OV, PRAD, SKCM, STAD, THCA, UCEC} available in this version.
}
}
\value{
A data frame with a single line containing the annotation of the RACS indicated in input.
}
\references{
[1] Iorio, F. et al. A Landscape of Pharmacogenomic Interactions in Cancer. Cell 166, 740–754 (2016).
}
\author{
Hanna Najgebauer and Francesco Iorio
}
\seealso{
\code{\link{CELLector.CFEs}},

\code{\link{CELLector.CFEs.CNAid_decode}}
}
\examples{
data(CELLector.CFEs.CNAid_decode)
CELLector.cna_look_up(cna_ID='cna26',
                      cnaId_decode = CELLector.CFEs.CNAid_decode,
                      TCGALabel = 'BRCA')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{annotation/decoding}% use one of  RShowDoc("KEYWORDS")
