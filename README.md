
![alt text](https://github.com/francescojm/CELLector_App/blob/master/www/Cellector-logo-double-size.png)

# CELLector: Genomics Guided Selection of Cancer in vitro models

### News: v2.0.0
From version 2.0.0 CELLector introduces a partitioned search space that is derived from the hierarchical search space implemented in the previous version. In particular, the hierarchical version identifies K recursive co-occurent signature rules. In the partitioned version, these are additionally converted into a K+1 non-overlapping groups, covering the entire sample space.

v1.2.0

Najgebauer, H., Yang, M., Francies, H., Pacini, C., Stronach, E. A., Garnett, M. J., Saez-Rodriguez, J., & Iorio, F. CELLector: Genomics Guided Selection of Cancer in vitro Models. http://doi.org/10.1101/275032

2020-03-03

CELLector is a computational tool for selecting the most clinically relevant cancer cell lines to be included in a new in-vitro study (or to be considered in a retrospective study), in a patient-genomic guided fashion.

CELLector combines methods from graph theory and market basket analysis; it leverages tumour genomics data to explore, rank, and select optimal cell line models in a user-friendly way, through the [CELLector Rshiny App](https://github.com/francescojm/CELLector_app). This enables making appropriate and informed choices about model inclusion/exclusion in retrospective analyses, future studies and it makes possible bridging cancer patient genomics with public available databases froom cell line based functional/pharmacogenomic screens, such as [CRISPR-cas9 dependency datasets](https://score.depmap.sanger.ac.uk/) and [large-scale in-vitro drug screens](https://www.cancerrxgene.org/).

Furthermore, CELLector includes interface functions to synchronise built-in cell line annotations and genomics data to their latest versions from the [Cell Model Passports](https://cellmodelpassports.sanger.ac.uk/) resource. Through this interface, bioinformaticians can quickly generate binary genomic event matrices (BEMs) accounting for hundreds of cancer cell lines, which can be used in systematic statistical inferences, associating patient-defined cell line subgroups with drug-response/gene-essentiality, for example through [GDSC tools](https://gdsctools.readthedocs.io/en/master/).

Additionally, CELLector allows the selection of models within user-defined contexts, for example, by focusing on genomic alterations occurring in biological pathways of interest or considering only predetermined sub-cohorts of cancer patients. 

Finally, CELLector identifies combinations of molecular alterations underlying disease subtypes currently lacking representative cell lines, providing guidance for the future development of new cancer models.

License: GNU GPLv3

![alt text](https://github.com/najha/CELLector_App/blob/master/www/cellcultures.jpg)

## Running Modalities

CELLector can be used in three different modalities:
  - (i) as an R package (within R, code available at: https://github.com/francescojm/CELLector)
  [Package vignette](https://rpubs.com/francescojm/CELLector)
  
  - (ii) as an online R shiny App (available at: https://ot-cellector.shinyapps.io/CELLector_App/)
  
  - (iii) as a local R shiny App (within Rstudio, code available at: https://github.com/francescojm/CELLector_App)

![alt text](www/preview1.png)

## R package: quick start (interactive vignette)

http://rpubs.com/francescojm/CELLector

![alt text](www/preview2.png)

## R package: Reference manual

https://github.com/francescojm/CELLector/blob/master/CELLector.pdf




