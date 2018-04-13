


library(CELLector)
data(CELLector.PrimTum.BEMs)
data(CELLector.Pathway_CFEs)
data(CELLector.CFEs.CNAid_mapping)
data(CELLector.CFEs.CNAid_decode)
data(CELLector.HCCancerDrivers)
data(CELLector.CellLine.BEMs)


### Change the following two lines to work with a different cancer type
tumours_BEM<-CELLector.PrimTum.BEMs$COREAD
CELLlineData<-CELLector.CellLine.BEMs$COREAD

### unicize the sample identifiers for the tumour data
tumours_BEM<-CELLector.unicizeSamples(tumours_BEM)

### building a CELLector searching space focusing on three pathways and TP53 wild-type patients only
CSS<-CELLector.Build_Search_Space(ctumours = t(tumours_BEM),
                                  verbose = FALSE,
                                  minGlobSupp = 0.05,
                                  cancerType = 'COREAD',
                                  pathwayFocused = c("RAS-RAF-MEK-ERK / JNK signaling",
                                                     "PI3K-AKT-MTOR signaling",
                                                     "WNT signaling"),
                                  pathway_CFEs = CELLector.Pathway_CFEs,
                                  cnaIdMap = CELLector.CFEs.CNAid_mapping,
                                  cnaIdDecode = CELLector.CFEs.CNAid_decode,
                                  cdg = CELLector.HCCancerDrivers,
                                  subCohortDefinition='TP53',
                                  NegativeDefinition=TRUE)

tumours_BEM<-CELLector.PrimTum.BEMs$SKCM
CELLlineData<-CELLector.CellLine.BEMs$SKCM

### unicize the sample identifiers for the tumour data
tumours_BEM<-CELLector.unicizeSamples(tumours_BEM)


CSS<-CELLector.Build_Search_Space(ctumours = t(tumours_BEM),
                                  verbose = FALSE,
                                  minGlobSupp = 0.05,
                                  cancerType = 'SKCM',
                                  pathwayFocused = 'RAS-RAF-MEK-ERK / JNK signaling',
                                  pathway_CFEs = CELLector.Pathway_CFEs,
                                  cnaIdMap = CELLector.CFEs.CNAid_mapping,
                                  cnaIdDecode = CELLector.CFEs.CNAid_decode,
                                  cdg = CELLector.HCCancerDrivers,
                                  subCohortDefinition='BRAF',
                                  NegativeDefinition=FALSE)

CELLector.solveFormula('~APC BRAF',dataset = CELLlineData)

### take all the signatures from the searching space
Signatures <- CELLector.createAllSignatures(CSS$navTable)

ModelMat<-CELLector.buildModelMatrix(Signatures$ES,CELLlineData,CSS$navTable)

selectedCellLines<-CELLector.makeSelection(modelMat = ModelMat,
                                           n=10,
                                           searchSpace = CSS$navTable)



### map cell lines onto subtypes, based on the signatures
MODELS<-vector()
for (cc in 1:length(Signatures$ES)){
  solved<-CELLector.solveFormula(Signatures$ES[[cc]],dataset = CELLlineData)
  MODELS[cc]<-paste(sort(solved$PS),collapse=', ')
}

### visit the searching space for the selection
visit<-CELLector.selectionVisit(CSS$navTable)

### put the cell lines in the same order in which the corresponding subtypes are
### encountered in the visit of the searching space
sortedModels<-MODELS[visit]

### ignore subtypes with not mached models (but actually this is what is most interesting for you, right?)
NODEidx<-visit[which(sortedModels!='')]
sortedModels<-sortedModels[which(sortedModels!='')]
modelMat<-CELLector.buildModelMatrix(sortedModels)

### select n = 10 cell lines, and put results in res
res<-CELLector.makeSelection(modelMat,n=10)
res$modelAccounted<-NODEidx[res$modelAccounted]
colnames(res)<-c('Tumour SubType Index','Representative Cell Line')



