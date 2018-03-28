


library(sunburstR)
library(igraph)
library(collapsibleTree)

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
colnames(tumours_BEM)<-paste(colnames(tumours_BEM),'_',1:ncol(tumours_BEM),sep='')

### building a CELLector searching space
tmp<-CELLector.Build_Search_Space(ctumours = t(tumours_BEM),
                                  verbose = FALSE,
                                  minGlobSupp = 0.05,
                                  cancerType = 'COREAD',
                                  pathwayFocused = c("RAS-RAF-MEK-ERK / JNK signaling","PI3K-AKT-MTOR signaling","WNT signaling"),
                                  mutOnly = FALSE,
                                  pathway_CFEs = CELLector.Pathway_CFEs,
                                  cnaIdMap = CELLector.CFEs.CNAid_mapping,
                                  cnaIdDecode = CELLector.CFEs.CNAid_decode,
                                  cdg = CELLector.HCCancerDrivers,
                                  subCohortDefinition='TP53',
                                  NegativeDefinition=TRUE)


### take all the signatures from the searching space
S <- CELLector.createAllSignatures(tmp$navTable)
encodedSignatures<-S$ES

### takes just the numerical part of the CELLlineData, ignore it, this will be fixed in the package
r<-CELLlineData[,2]
COSMICids<-CELLlineData[,1]
CELLlineData<-CELLlineData[,3:ncol(CELLlineData)]
rownames(CELLlineData)<-r

### map cell lines onto subtypes, based on the signatures
MODELS<-vector()
for (cc in 1:length(encodedSignatures)){
  solved<-CELLector.solveFormula(encodedSignatures[[cc]],dataset = CELLlineData)
  MODELS[cc]<-paste(sort(solved$PS),collapse=', ')
}

### visit the searching space for the selection
visit<-CELLector.selectionVisit(tmp$navTable)

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



