#library(CELLector)

data(CELLector.PrimTum.BEMs)
data(CELLector.Pathway_CFEs)
data(CELLector.CFEs.CNAid_mapping)
data(CELLector.CFEs.CNAid_decode)
data(CELLector.HCCancerDrivers)

data(CELLector.PathwayList)

tumours_BEM<-CELLector.PrimTum.BEMs$COREAD

CELLector.mostSupported_CFEs(t(tumours_BEM),minlen = 1)


tmp<-CELLector.Build_Search_Space(ctumours = t(tumours_BEM),verbose = FALSE,
                             minGlobSupp = 0.01,
                             cancerType = 'COREAD',
                             pathwayFocused = c("RAS-RAF-MEK-ERK / JNK signaling","PI3K-AKT-MTOR signaling","WNT signaling"),
                             pathway_CFEs = CELLector.Pathway_CFEs,
                             cnaIdMap = CELLector.CFEs.CNAid_mapping,
                             cnaIdDecode = CELLector.CFEs.CNAid_decode,
                             cdg = CELLector.HCCancerDrivers,
                             subCohortDefinition='TP53',
                             NegativeDefinition=TRUE)


tmp<-CELLector.Build_Search_Space(ctumours = t(tumours_BEM),
                             minGlobSupp = 0.01,
                             cancerType = 'COREAD',
                             pathwayFocused = c("RAS-RAF-MEK-ERK / JNK signaling","PI3K-AKT-MTOR signaling","WNT signaling"),
                             pathway_CFEs = CELLector.Pathway_CFEs,
                             cnaIdMap = CELLector.CFEs.CNAid_mapping,
                             cnaIdDecode = CELLector.CFEs.CNAid_decode,
                             cdg = CELLector.HCCancerDrivers,
                             FeatureToExclude = 'TP53')


CELLector.cna_look_up(cna_ID = "cna52",cnaId_decode = CELLector.CFEs.CNAid_decode,TCGALabel = 'COREAD')
