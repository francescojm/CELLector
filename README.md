# CELLector

CELLector is a computational tool to assist the selection of the most relevant cancer cell lines to be included in a new in-vitro study (or to be considered in a retrospective study) in a genomic-guided fashion. CELLector combines methods from graph theory and market basket analysis; it leverages tumour genomics data to explore, rank, and select optimal cell line models in a user-friendly way, enabling scientists to make appropriate and informed choices about model inclusion/exclusion in retrospective analyses and future studies. Additionally, it allows the selection of models within user-defined contexts, for example, by focusing on genomic alterations occurring in biological pathways of interest or considering only predetermined sub-cohorts of cancer patients. Finally, CELLector identifies combinations of molecular alterations underlying disease subtypes currently lacking representative cell lines, providing guidance for the future development of new cancer models.\

License: MIT

## Running Modalities

CELLector can be used in three different modalities: (i) as an R package (within R, code available at: https://github.com/francescojm/CELLector), (ii) as an online R shiny App (available at: https://ot-cellector.shinyapps.io/cellector_app/ - temporary deployment), (iii) running the R shiny App locally (within Rstudio, code available at: https://github.com/francescojm/CELLector_App).


## R package: quick start (interactive vignette)

http://rpubs.com/francescojm/CELLector
