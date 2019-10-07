\name{CELLector.mostSupported_CFEs}
\alias{CELLector.mostSupported_CFEs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Most recurrent combinations of Cancer Functional Events
}
\description{
This function identifies the most frequent combination of cancer functional events (CFEs) in a large
cohort of cancer patients.
}
\usage{
CELLector.mostSupported_CFEs(transactions,
                             minSupport = 0.05,
                             minlen = 1,
                             maxLen = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{transactions}{
  A named binary matrix with CFEs on the rows, samples on the columns and entries specifying the
  presence/absence of a given CFE in a given sample: the \emph{transactions} object.
}
  \item{minSupport}{
  The minimal support that a combination of CFEs must have, i.e. the minimal ratio of samples in which the CFEs must be observed simoultanously, in order to be considered in the analysis.
}
  \item{minlen}{
  The minimal length of a combination of CFEs (of how many indivudal CFE it needs to be composed) in order to be considered in the analysis (1 by default).
}
  \item{maxLen}{
The maximal length of a combination of CFEs (the maximal number of indivudal CFEs) in order to be considered in the analysis (10 by default).
}
}
\details{
  This function uses the \emph{eclat} function from the R package \emph{arules}.
}
\value{
A list with the following fields:
  \item{MSIS}{A string or a vector of strings (depending on the argument \code{minlen}) specifying the CFE (or the combination of individual CFEs) that is the most frequently observed (simultaneously across the samples in input)}
   \item{SUPPORT}{The ratio of samples where the combination of CFEs in \code{MSIS} is obaserved on the total number of samples, i.e. number of columns in the \code{transactions} argument}
   \item{absSUPPORT}{The number of samples where the combination of CFEs in \code{MSIS} is obaserved}
   \item{supportingSamples}{The identifiers of the samples supporting \code{MSIS}, i.e. the names of the columns of \emph{transactions}, in which the entries corresponding to \code{MSIS} rows are equal to 1.}
}
\references{
Najgebauer et al., CELLector: Genomics Guided Selection of Cancer in vitro Models.

doi:10.1101/275032
}
\author{
Hanna Najgebauer and Francesco Iorio
}
\examples{
data(CELLector.PrimTum.BEMs)
RES<-CELLector.mostSupported_CFEs(transactions = t(CELLector.PrimTum.BEMs$COREAD),
                             minlen = 2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{analysis}
