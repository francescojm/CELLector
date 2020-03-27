\name{CELLector.CELLline_buildBEM}
\alias{CELLector.CELLline_buildBEM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Building a Genomic Binary Event Matrix (BEM) for in-vitro models
}
\description{
This function takes in input a catalogue of somatic genomic variants observed in cancer cell lines (or it uses the most recent catalogue of somatic genomic variants from the Cell Model Passports [1]) and it converts it into a presence/absence (binary) matrix, which can be processed by the CELLector package and CELLector shiny app for mapping cell lines onto cancer-patient derived genomic subtypes.
}
\usage{
CELLector.CELLline_buildBEM(varCat = NULL,
                            Tissue,
                            Cancer_Type,
                            Cancer_Type_details = NULL,
                            sample_site = NULL,
                            excludeOrganoids = FALSE,
                            humanonly = TRUE,
                            msi_status_select = NULL,
                            gender_select = NULL,
                            mutational_burden_th = NULL,
                            age_at_sampling = NULL,
                            ploidy_th = NULL,
                            ethnicity_to_exclude = NULL,
                            GenesToConsider = NULL,
                            VariantsToConsider = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{varCat}{
  A data frame containing a catalogue of somatic genomic variants observed in cancer cell lines with one row per variant. At least the following column/headers should be included: \code{model_id} (the Cell Model Passport identifier [1] of the cell line in which the variant is observed), \code{gene_symbol} (HUGO [2] symbol of the gene hosting the variant), \code{cdna_mutation} (variant specification based on the longest transcript of the gene), \code{aa_mutation} (aminoacid substitution). See format of the \code{CELLector.RecfiltVariants} object for details on the syntax of the variant and aminoacid substitution specifications.
When this parameter is different from \code{NULL} all the others parameters are ignored by this function, with the exception of \code{GenesToConsider} and \code{VariantsToConsider}. By default, this parameter is set to \code{NULL} and the most up-to-date variant catalogue from the Cell Model Passports [1] is used insted. In this case the other parameters define which samples to include in the final genomic binary event matrix.
}
  \item{Tissue}{
  String specifying the tissue of origin of the cell lines that should be considering while assembling the final genomic binary event matrix. Possible values are those specified in the \code{tissue} column of the Cell Model Passport [1] annotation object, downloadable through the \code{CELLector.CMPs_getModelAnnotation} function.
      This parameter is used only when \code{varCat} is not set to its default \code{NULL} value, therefore the most-up-to-date variant catalogue from the
  Cell Model Passports [1] is used instead of a variant catalogue provided in input.
}
  \item{Cancer_Type}{
  String (or vector of strings) specifying the cancer type(s) of origin of the cell lines that should be considering while assembling the final genomic binary event matrix. Possible values are those specified in the \code{cancer_type} column of the Cell Model Passport [1] annotation object, downloadable through the \code{CELLector.CMPs_getModelAnnotation} function.
    This parameter is used only when \code{varCat} is not set to its default \code{NULL} value, therefore the most-up-to-date variant catalogue from the
  Cell Model Passports [1] is used instead of a variant catalogue provided in input.
}
  \item{Cancer_Type_details}{
  String (or vector of strings) specifying the cancer type detail(s) of the cell lines that should be considering while assembling the final genomic binary event matrix. Possible values are those specified in the \code{cancer_type_detail} column of the Cell Model Passport [1] annotation object, downloadable through the \code{CELLector.CMPs_getModelAnnotation} function.
    This parameter is used only when \code{varCat} is not set to its default \code{NULL} value, therefore the most-up-to-date variant catalogue from the
  Cell Model Passports [1] is used instead of a variant catalogue provided in input.
}
  \item{sample_site}{
  String (or vector of strings) specifying the sample site(s) of the cell lines that should be considering while assembling the final genomic binary event matrix. Possible values are those specified in the \code{sample_site} column of the Cell Model Passport [1] annotation object, downloadable through the \code{CELLector.CMPs_getModelAnnotation} function.
    This parameter is used only when \code{varCat} is not set to its default \code{NULL} value, therefore the most-up-to-date variant catalogue from the
  Cell Model Passports [1] is used instead of a variant catalogue provided in input.
}
  \item{excludeOrganoids}{
  Boolean parameter specifying whether organoid models should not be considered while assembling the final genomic binary event matrix.
    This parameter is used only when \code{varCat} is not set to its default \code{NULL} value, therefore the most-up-to-date variant catalogue from the
  Cell Model Passports [1] is used instead of a variant catalogue provided in input.
}
  \item{humanonly}{
  Boolean parameter specifying whether human cell lines only should not be considered while assembling the final genomic binary event matrix.
    This parameter is used only when \code{varCat} is not set to its default \code{NULL} value, therefore the most-up-to-date variant catalogue from the
  Cell Model Passports [1] is used instead of a variant catalogue provided in input.
}
  \item{msi_status_select}{
 String parameter specifying whether microsatellite stable or instable cell lines only should not be considered while assembling the final genomic binary event matrix.
   Possible values for this parameters are "MSI" (for instable cell lines) and "MSS" (for stable cell lines).
    This parameter is used only when \code{varCat} is not set its default \code{NULL} value, therefore the most-up-to-date variant catalogue from the
  Cell Model Passports [1] is used instead of a variant catalogue provided in input.
}


\item{gender_select}{
String (or vector of strings) specifying the gendre of the patients the cell lines that should be considering while assembling the final genomic binary event matrix were resected from. Possible values are those specified in the \code{gender} column of the Cell Model Passport [1] annotation object, downloadable through the \code{CELLector.CMPs_getModelAnnotation} function.
    This parameter is used only when \code{varCat} is not set to its default \code{NULL} value, therefore the most-up-to-date variant catalogue from the
  Cell Model Passports [1] is used instead of a variant catalogue provided in input.
}
  \item{mutational_burden_th}{
  A vector of two real numbers specifying an interval of number of mutations per million base pairs the mutational burden of the cell lines that should be considering while assembling the final genomic binary event matrix should fall in.
    This parameter is used only when \code{varCat} is not set to its default \code{NULL} value, therefore the most-up-to-date variant catalogue from the
  Cell Model Passports [1] is used instead of a variant catalogue provided in input.
}
  \item{age_at_sampling}{
A vector of two real numbers specifying an age at sampling interval the cell lines that should be considering while assembling the final genomic binary event matrix should fall in.
    This parameter is used only when \code{varCat} is not set to its default \code{NULL} value, therefore the most-up-to-date variant catalogue from the
  Cell Model Passports [1] is used instead of a variant catalogue provided in input.
}
  \item{ploidy_th}{
  A vector of two real numbers specifying an interval the ploidy of the cell lines that should be considering while assembling the final genomic binary event matrix should fall in.
    This parameter is used only when \code{varCat} is not set to its default \code{NULL} value, therefore the most-up-to-date variant catalogue from the
  Cell Model Passports [1] is used instead of a variant catalogue provided in input.
}
  \item{ethnicity_to_exclude}{
  String (or vector of strings) specifying the ethnic groups of the patients the cell lines that should not be considering while assembling the final genomic binary event matrix were resected from. Possible values are those specified in the \code{ethnicity} column of the Cell Model Passport [1] annotation object, downloadable through the \code{CELLector.CMPs_getModelAnnotation} function.
    This parameter is used only when \code{varCat} is not set to its default \code{NULL} value, therefore the most-up-to-date variant catalogue from the
  Cell Model Passports [1] is used instead of a variant catalogue provided in input.
}
  \item{GenesToConsider}{
 A list of strings with HGNC symbols [2] for genes hosting the variants to be extracted from the catalogue and assembled into the final matrix. When set to its default \code{NULL} value, all genes hosting at least one variants are considered.
}
  \item{VariantsToConsider}{
  A list of individual somatic variants to be extracted from the catalogue and assembled into the final matrix. The format should be the same of the \code{CELLector.RecfiltVariants}.
   }
}
\value{
A data frame with one row per cell lines, Cell model passport identifiers and cell line names in the first two columns followed by a presence/absence (binary) matrix with gene symbols on the columns, specifying in the \eqn{i,j-}entry the status of the \eqn{j}th gene in the \eqn{i}th cell line, i.e. 0 = wild-type, 1 = mutated.
}
\references{
[1] van der Meer D, Barthorpe S, Yang W, et al. Cell Model Passports-a hub for clinical, genetic and functional datasets of preclinical cancer models. Nucleic Acids Res. 2019;47(D1):D923–D929. doi:10.1093/nar/gky872

[2] Braschi, B. et al. Genenames.org: the HGNC and VGNC resources in 2019. Nucleic Acids Res. Epub 2018 Oct 10. PMID: 30304474 DOI: 10.1093/nar/gky930

[3] Iorio, F. et al. A Landscape of Pharmacogenomic Interactions in Cancer. Cell 166, 740–754 (2016).
}
\author{
Francesco Iorio (fi9323@gmail.com)
}
\seealso{
\code{\link{CELLector.CMPs_getModelAnnotation}}, \code{\link{CELLector.Tumours_buildBEM}}
}
\examples{
## loading high-confidence cancer driver genes from [1]
data(CELLector.HCCancerDrivers)

## loading COSMIC [3] variants observed it at least two patients from [3]
data(CELLector.RecfiltVariants)

## Assembling a genomic binary event matrix (BEM) for human derived microsatellite stable diploid colorectal carcinoma cell lines,
resected from male patients using genomic data from the Cell Model Passports [1]

COREAD_cl_BEM <- CELLector.CELLline_buildBEM(
                  Tissue='Large Intestine',
                  Cancer_Type = 'Colorectal Carcinoma',
                  msi_status_select = 'MSI',
                  ploidy_th = c(2,2),
                  GenesToConsider = CELLector.HCCancerDrivers,
                  VariantsToConsider = CELLector.RecfiltVariants)

## showing first entries of the BEM
head(COREAD_cl_BEM)

## showing a bar diagram with mutation frequencies in the BEM
barplot(sort(colSums(COREAD_cl_BEM[,3:ncol(COREAD_cl_BEM)]),decreasing=TRUE)[1:20],
        las=2,ylab='n. mutated cell lines')

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{analysis}% use one of  RShowDoc("KEYWORDS")
\keyword{annotation/decoding}% __ONLY ONE__ keyword per line
