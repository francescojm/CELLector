\name{CELLector.from_Hierarchical_to_Partition}
\alias{CELLector.from_Hierarchical_to_Partition}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert CELLector output to partioned format
}
\description{
This function goes through binary tree output of CELLector search space in the form of navigable table [1] and convert hierachical configuration to partitioned.}
\usage{
CELLector.from_Hierarchical_to_Partition(navTable,
                                         samples_id,
                                         verbose=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{navTable}{
    \code{navTable} entry of output from \code{\link{CELLector.Build_Search_Space}}
}
  \item{samples_id}{
   Sample names of the cancer patient cohort as input of \code{\link{CELLector.Build_Search_Space}}, i.e. \code{rownames(ctumours)}
}
  \item{verbose}{
    A boolean argument specifying whether step-by-step information on the algorithm progression should be displayed run-time
}
}
\details{
This is an auxiliary function used in \code{\link{CELLector.Build_Search_Space_Partitioned}} (for details check documentation).
}
\value{
\emph{Navigable table}: a data frame with a row for each patient group and columns as described in \code{\link{CELLector.Build_Search_Space_Partitioned}}
}
\references{
[1] Najgebauer, H. et al. Genomics Guided Selection of Cancer in vitro Models.

https://doi.org/10.1016/j.cels.2020.04.007
}
\author{
Lucia Trastulla and Francesco Iorio
}

\seealso{
\code{\link{CELLector.Build_Search_Space}}
\code{\link{CELLector.Build_Search_Space_Partitioned}}
}
\examples{
# see example for CELLector.Build_Search_Space_Partitioned
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{analysis}

