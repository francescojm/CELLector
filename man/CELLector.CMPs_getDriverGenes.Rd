\name{CELLector.CMPs_getDriverGenes}
\alias{CELLector.CMPs_getDriverGenes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  List of high-confidence cancer driver genes from the Cell Model Passports
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
CELLector.CMPs_getDriverGenes
(URL = "https://cog.sanger.ac.uk/cmp/download/cancer_genes_latest.csv.gz")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{URL}{
The URL specifying the online location of the file containing the high-confidence cancer driver genes. Default value points to the most up-to-date version of this file (which is keptupdated).
}
}
\value{
A vector of strings (containinng HGNC symbols) with one entry per high-confidence cancer driver genes.
}
\references{
van der Meer D, Barthorpe S, Yang W, et al. Cell Model Passports-a hub for clinical, genetic and functional datasets of preclinical cancer models. Nucleic Acids Res. 2019;47(D1):D923–D929. doi:10.1093/nar/gky872
}
\author{
Francesco Iorio (fi9323@gmail.com)
}
\examples{
CMP_cdgenes<-CELLector.CMPs_getDriverGenes()
  head(CMP_cdgenes)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Cell Model Passports Interface}% use one of  RShowDoc("KEYWORDS")

