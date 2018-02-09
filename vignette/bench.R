#library(CELLector)

library(sunburstR)
library(collapsibleTree)

data(CELLector.PrimTum.BEMs)
data(CELLector.Pathway_CFEs)
data(CELLector.CFEs.CNAid_mapping)
data(CELLector.CFEs.CNAid_decode)
data(CELLector.HCCancerDrivers)


tumours_BEM<-CELLector.PrimTum.BEMs$COREAD

CELLector.mostSupported_CFEs(t(tumours_BEM),minlen = 1)


tmp<-CELLector.Build_Search_Space(ctumours = t(tumours_BEM),
                                  verbose = FALSE,
                             minGlobSupp = 0.01,
                             cancerType = 'COREAD',
                             pathwayFocused = c("RAS-RAF-MEK-ERK / JNK signaling","PI3K-AKT-MTOR signaling","WNT signaling"),
                             mutOnly = FALSE,
                             pathway_CFEs = CELLector.Pathway_CFEs,
                             cnaIdMap = CELLector.CFEs.CNAid_mapping,
                             cnaIdDecode = CELLector.CFEs.CNAid_decode,
                             cdg = CELLector.HCCancerDrivers,
                             subCohortDefinition='TP53',
                             NegativeDefinition=TRUE)

CELLector.solveFormula('cna27',CELLector.CellLine.BEMs[['COREAD']])

seq<-CELLector_App.sunBurstFormat(tmp$navTable)

tmpCol <- Get(Traverse(tmp$TreeRoot,traversal = 'level'),'Colors')

ttmp<-tmpCol

names(ttmp)<-NULL

nvoid<-grep('Others',sequences$V1,value = TRUE)

stpes<-nvoid

colors <- list(
  domain=c('0 TOTAL',names(tmpCol),stpes),
  range=c('black',ttmp,rep('white',length(stpes)))
)



CELLector.createAllSignatures(tmp$navTable)


tmp<-CELLector.Build_Search_Space(ctumours = t(tumours_BEM),
                             minGlobSupp = 0.05,
                             cancerType = 'COREAD',
                             pathwayFocused = c("RAS-RAF-MEK-ERK / JNK signaling","PI3K-AKT-MTOR signaling","WNT signaling"),
                             pathway_CFEs = CELLector.Pathway_CFEs,
                             cnaIdMap = CELLector.CFEs.CNAid_mapping,
                             cnaIdDecode = CELLector.CFEs.CNAid_decode,
                             cdg = CELLector.HCCancerDrivers,
                             FeatureToExclude = 'TP53')

CC <- colors(distinct = TRUE)
CC <- CC[setdiff(1:length(CC),grep('gray',CC))]
CC <- rgb(t(col2rgb(CC)),maxColorValue = 255)

COLORSbyLev <- CC[sample(length(CC))][1:tmp$TreeRoot$totalCount]

RelatesToFatherAs <- rep('-',tmp$TreeRoot$totalCount)
RelatesToFatherAs[which(Get(Traverse(tmp$TreeRoot,traversal = 'level'),attribute = 'NodeType')=='Right.Child')]<-'Complement'
RelatesToFatherAs[which(Get(Traverse(tmp$TreeRoot,traversal = 'level'),attribute = 'NodeType')=='Left.Child')]<-'Refinement'

tmp$TreeRoot$Set(Colors=COLORSbyLev,traversal = 'level')
tmp$TreeRoot$Set(RelatesToFatherAs=RelatesToFatherAs,traversal = 'level')

collapsibleTree(tmp$TreeRoot,
                fill = 'Colors',
                inputId = 'searchSpace',
                tooltip = TRUE,
                attribute = 'RelatesToFatherAs')

tmpCol <- Get(Traverse(tmp$TreeRoot,traversal = 'level'),'Colors')


colors <- list(
     domain=names(tmpCol),
     range=tmpCol
  )


tmp<-CELLector.sunBurstFormat(tmp$navTable)


sb<-sunburst(tmp$sequences,colors = tmp$colors,breadcrumb = list(w = 100),percent = FALSE,count = FALSE,
             explanation = "function(d) { return d.data.name}")

sb$x$tasks <- htmlwidgets::JS(
  'function(){d3.select(this.el).select("#" + this.el.id + "-trail").style("font-size","60%")}'
)

sb



