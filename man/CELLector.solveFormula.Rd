\name{CELLector.solveFormula}
\alias{CELLector.solveFormula}
\title{
Identify cell lines harbouring a signature of Cancer Functional Events
}
\description{
This function takes in input a signature of Cancer Functional Events (CFEs, defined in [1]) as outputted by the \code{CELLector.createAllSignatures} function, and a binary event matrix (BEM) modeling the presence/absence of all the CFEs across a set of immortalised human cancer cell lines. It returns the set of cell line collectively harbouring the inputted signature.
}
\usage{
CELLector.solveFormula(RULE, dataset, To_beExcluded = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{RULE}{
  A string representing a signature of cancer functional events (CFEs, defined in [1]), i.e. names of CFEs space separated (and possibly negated ~)
}
  \item{dataset}{
  A data frame in which the first two columns contain the COSMIC [2] identiefiers and names of cell lines (one per row), respectively, and then binary entries indicating the status of each CFEs (one per column) across cell lines. The format is the same of the entries of the list in the built-in \code{CELLector.CellLine.BEMs} object.
}
  \item{To_beExcluded}{
  If different from NULL (default value), then this must be a list of strings with cell line names that should be excluded a priori from the output.
}
}
\value{
  A list with the following entries:
    \item{PS}{Positive samples: names of the cell lines collectively harbouring the signature of CFEs provided in input}
    \item{N}{The number of cell lines in \code{PS}}
    \item{PERC}{The number of cell lines collectively harbouring the signatures of CFEs provided in input as ratio of the total number of cell lines in the inputted BEM}
}
\references{
[1] Iorio, F. et al. A Landscape of Pharmacogenomic Interactions in Cancer. Cell 166, 740–754 (2016).

[2] Forbes, S. A. et al. COSMIC: exploring the world’s knowledge of somatic mutations in human cancer.
Nucleic Acids Res. 43, D805–11 (2015).
}
\author{
Hanna Najgebauer and Francesco Iorio
}
\seealso{
\code{\link{CELLector.createAllSignatures}},
\code{\link{CELLector.CellLine.BEMs}}
}
\examples{
data(CELLector.CellLine.BEMs)

### Selecting colorectal cancer cell lines that are
### APC wild type and BRAF mutant
CellLines<-CELLector.solveFormula('~APC BRAF',
                                  CELLector.CellLine.BEMs$COREAD)

CellLines
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{analysis}% use one of  RShowDoc("KEYWORDS")
