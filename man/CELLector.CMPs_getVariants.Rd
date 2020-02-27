\name{CELLector.CMPs_getVariants}
\alias{CELLector.CMPs_getVariants}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Somatic variants of in-visto models from the Cell Model Passports
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
CELLector.CMPs_getVariants
(URL = "https://cog.sanger.ac.uk/cmp/download/mutations_2018-08-01_1640.csv.gz")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{URL}{
The URL specifying the online location of the file containing the somatic variants. Default value points to the most up-to-date version of this file (which is keptupdated).

}
}
\value{
A data frame with one entry per variant and columns specifying different annotation entries.
}
\references{
van der Meer D, Barthorpe S, Yang W, et al. Cell Model Passports-a hub for clinical, genetic and functional datasets of preclinical cancer models. Nucleic Acids Res. 2019;47(D1):D923–D929. doi:10.1093/nar/gky872
}
\author{
Francesco Iorio (fi9323@gmail.com)
}
\examples{
CMP_variants<-CELLector.CMPs_getVariants()
  head(CMP_variants)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Cell Model Passports Interface}% use one of  RShowDoc("KEYWORDS")
