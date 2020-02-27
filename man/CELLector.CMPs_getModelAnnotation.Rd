\name{CELLector.CMPs_getModelAnnotation}
\alias{CELLector.CMPs_getModelAnnotation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  In-vitro models' annotation from the Cell Model Passports
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
CELLector.CMPs_getModelAnnotation
(URL = "https://cog.sanger.ac.uk/cmp/download/model_list_latest.csv.gz")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{URL}{
The URL specifying the online location of the annotation file. Default value points to the most up-to-date version of this file (which is keptupdated).
}
}
\details{
}
\value{
  A data frame with one row per model and one colum per annotation entry.
}
\references{
van der Meer D, Barthorpe S, Yang W, et al. Cell Model Passports-a hub for clinical, genetic and functional datasets of preclinical cancer models. Nucleic Acids Res. 2019;47(D1):D923–D929. doi:10.1093/nar/gky872
}
\author{
Francesco Iorio (fi9323@gmail.com)
}
\examples{
  CMP_annotation<-CELLector.CMPs_getModelAnnotation()
  head(CMP_annotation)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Cell Model Passports Interface}% use one of  RShowDoc("KEYWORDS")
