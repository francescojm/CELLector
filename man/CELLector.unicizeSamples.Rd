\name{CELLector.unicizeSamples}
\alias{CELLector.unicizeSamples}
\title{
   Unicize patient samples' identifiers
}
\description{
  This function checks if there are multiple samples derived from the same patients in the
  binary event matrix (BEM) modeling the presence/absence of the cancer functional events (CFEs, defined in [1]), in the cancer patients. These can be mantained (and in this case their identifier will be made unique) or
  discarded
}
\usage{
CELLector.unicizeSamples(ctumours,
                         keepReplicates = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ctumours}{
  A binary matrix with entries indicating the status (Present/Absent) of each CFE (one per row) across primary tumors samples (one per column).
}
  \item{keepReplicates}{
  A boolean value indicating whether the duplicated samples should be kept (and their identifier made unique, by adding a progressive numerical suffix) or discarded (in this case only one sample per patient will be kept and identifiers unchanged).
}
}
\value{
   A binary matrix with entries indicating the status (Present/Absent) of each CFE (one per row) across primary tumors samples (one per column), and with unique patients' (column) identifiers
}
\references{
[1] Iorio, F. et al. A Landscape of Pharmacogenomic Interactions in Cancer. Cell 166, 740–754 (2016).
}
\author{
Hanna Najgebauer and Francesco Iorio
}

\seealso{
\code{\link{CELLector.PrimTum.BEMs}}
}
\examples{
data(CELLector.PrimTum.BEMs)

tumours_data<-CELLector.PrimTum.BEMs$COREAD

dim(tumours_data)
length(unique(colnames(tumours_data)))

tumours_data<-CELLector.unicizeSamples(tumours_data)

dim(tumours_data)
length(unique(colnames(tumours_data)))
}
\keyword{annotation/decoding}
