\name{CELLector.CELLline_buildBEM}
\alias{CELLector.CELLline_buildBEM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Building a Genomic Binary Event Matrix (BEM) for in-vitro models
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
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
When this parameter is different from \code{NULL} all the following parameters are ignored by this function, with the exception of \code{GenesToConsider} and \code{VariantsToConsider}.

}
  \item{Tissue}{

}
  \item{Cancer_Type}{
%%     ~~Describe \code{Cancer_Type} here~~
}
  \item{Cancer_Type_details}{
%%     ~~Describe \code{Cancer_Type_details} here~~
}
  \item{sample_site}{
%%     ~~Describe \code{sample_site} here~~
}
  \item{excludeOrganoids}{
%%     ~~Describe \code{excludeOrganoids} here~~
}
  \item{humanonly}{
%%     ~~Describe \code{humanonly} here~~
}
  \item{msi_status_select}{
%%     ~~Describe \code{msi_status_select} here~~
}
  \item{gender_select}{
%%     ~~Describe \code{gender_select} here~~
}
  \item{mutational_burden_th}{
%%     ~~Describe \code{mutational_burden_th} here~~
}
  \item{age_at_sampling}{
%%     ~~Describe \code{age_at_sampling} here~~
}
  \item{ploidy_th}{
%%     ~~Describe \code{ploidy_th} here~~
}
  \item{ethnicity_to_exclude}{
%%     ~~Describe \code{ethnicity_to_exclude} here~~
}
  \item{GenesToConsider}{
%%     ~~Describe \code{GenesToConsider} here~~
}
  \item{VariantsToConsider}{
%%     ~~Describe \code{VariantsToConsider} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
[1] van der Meer D, Barthorpe S, Yang W, et al. Cell Model Passports-a hub for clinical, genetic and functional datasets of preclinical cancer models. Nucleic Acids Res. 2019;47(D1):D923–D929. doi:10.1093/nar/gky872

[2] Braschi, B. et al. Genenames.org: the HGNC and VGNC resources in 2019. Nucleic Acids Res. Epub 2018 Oct 10. PMID: 30304474 DOI: 10.1093/nar/gky930

}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (varCat = NULL, Tissue, Cancer_Type, Cancer_Type_details = NULL,
    sample_site = NULL, excludeOrganoids = FALSE, humanonly = TRUE,
    msi_status_select = NULL, gender_select = NULL, mutational_burden_th = NULL,
    age_at_sampling = NULL, ploidy_th = NULL, ethnicity_to_exclude = NULL,
    GenesToConsider = NULL, VariantsToConsider = NULL)
{
    if (length(varCat) == 0) {
        varCat <- CELLector.CMPs_getVariants()
        clAnnotation <- CELLector.CMPs_getModelAnnotation()
        clAnnotation$cancer_type_detail <- str_sub(clAnnotation$cancer_type_detail,
            3, end = str_length(clAnnotation$cancer_type_detail) -
                3)
        if (!excludeOrganoids) {
            id <- which(clAnnotation$tissue == Tissue & is.element(clAnnotation$cancer_type,
                Cancer_Type))
        }
        else {
            id <- which(clAnnotation$tissue == Tissue & is.element(clAnnotation$cancer_type,
                Cancer_Type) & clAnnotation$model_type != "Organoid")
        }
        cls <- clAnnotation$model_id[id]
        varCat <- varCat[which(is.element(varCat$model_id, cls)),
            ]
        clAnnotation <- clAnnotation[which(is.element(clAnnotation$model_id,
            cls)), ]
        if (length(Cancer_Type_details) > 0) {
            id <- which(is.element(clAnnotation$cancer_type_detail,
                Cancer_Type_details))
            cls <- clAnnotation$model_id[id]
            varCat <- varCat[which(is.element(varCat$model_id,
                cls)), ]
            clAnnotation <- clAnnotation[which(is.element(clAnnotation$model_id,
                cls)), ]
        }
        if (length(sample_site) > 0) {
            id <- which(is.element(clAnnotation$sample_site,
                sample_site))
            cls <- clAnnotation$model_id[id]
            varCat <- varCat[which(is.element(varCat$model_id,
                cls)), ]
            clAnnotation <- clAnnotation[which(is.element(clAnnotation$model_id,
                cls)), ]
        }
        if (length(humanonly) > 0) {
            id <- which(clAnnotation$species == "Homo Sapiens")
            cls <- clAnnotation$model_id[id]
            varCat <- varCat[which(is.element(varCat$model_id,
                cls)), ]
            clAnnotation <- clAnnotation[which(is.element(clAnnotation$model_id,
                cls)), ]
        }
        if (length(msi_status_select) > 0) {
            id <- which(!is.na(clAnnotation$msi_status) & (clAnnotation$msi_status ==
                msi_status_select | (msi_status_select == "MSI-L/H" &
                (clAnnotation$msi_status == "MSI-L" | clAnnotation$msi_status ==
                  "MSI-H"))))
            cls <- clAnnotation$model_id[id]
            varCat <- varCat[which(is.element(varCat$model_id,
                cls)), ]
            clAnnotation <- clAnnotation[which(is.element(clAnnotation$model_id,
                cls)), ]
        }
        if (length(gender_select) > 0) {
            id <- which(is.element(clAnnotation$gender, gender_select))
            cls <- clAnnotation$model_id[id]
            varCat <- varCat[which(is.element(varCat$model_id,
                cls)), ]
            clAnnotation <- clAnnotation[which(is.element(clAnnotation$model_id,
                cls)), ]
        }
        if (length(ethnicity_to_exclude) > 0) {
            id <- which(!is.element(clAnnotation$ethnicity, ethnicity_to_exclude))
            cls <- clAnnotation$model_id[id]
            varCat <- varCat[which(is.element(varCat$model_id,
                cls)), ]
            clAnnotation <- clAnnotation[which(is.element(clAnnotation$model_id,
                cls)), ]
        }
        if (length(mutational_burden_th) > 0) {
            id <- which(round(clAnnotation$mutational_burden) >=
                mutational_burden_th[1] & round(clAnnotation$mutational_burden) <=
                mutational_burden_th[2])
            cls <- clAnnotation$model_id[id]
            varCat <- varCat[which(is.element(varCat$model_id,
                cls)), ]
            clAnnotation <- clAnnotation[which(is.element(clAnnotation$model_id,
                cls)), ]
        }
        if (length(ploidy_th) > 0) {
            id <- which(round(clAnnotation$ploidy) >= ploidy_th[1] &
                round(clAnnotation$ploidy) <= ploidy_th[2])
            cls <- clAnnotation$model_id[id]
            varCat <- varCat[which(is.element(varCat$model_id,
                cls)), ]
            clAnnotation <- clAnnotation[which(is.element(clAnnotation$model_id,
                cls)), ]
        }
        if (length(age_at_sampling) > 0) {
            id <- which(round(clAnnotation$age_at_sampling) >=
                age_at_sampling[1] & round(clAnnotation$age_at_sampling) <=
                age_at_sampling[2])
            cls <- clAnnotation$model_id[id]
            varCat <- varCat[which(is.element(varCat$model_id,
                cls)), ]
            clAnnotation <- clAnnotation[which(is.element(clAnnotation$model_id,
                cls)), ]
        }
        if (length(GenesToConsider) > 0) {
            varCat <- varCat[which(is.element(varCat$gene_symbol,
                GenesToConsider)), ]
        }
        if (length(VariantsToConsider) > 0) {
            sigs <- paste(varCat$gene_symbol, varCat$cdna_mutation,
                paste("p.", varCat$aa_mutation, sep = ""))
            varCat <- varCat[which(is.element(sigs, VariantsToConsider)),
                ]
        }
    }
    allModels <- sort(unique(varCat$model_id))
    allModel_ids <- varCat$model_id[match(allModels, varCat$model_id)]
    allGenes <- sort(unique(varCat$gene_symbol))
    BEM <- do.call(what = cbind, lapply(allModels, function(x) {
        is.element(allGenes, varCat$gene_symbol[varCat$model_id ==
            x]) + 0
    }))
    rownames(BEM) <- allGenes
    cls <- clAnnotation$model_name[match(allModel_ids, clAnnotation$model_id)]
    BEM <- data.frame(CMP_identifier = allModel_ids, CellLine = cls,
        t(BEM))
    return(BEM)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
