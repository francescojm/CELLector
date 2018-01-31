library(CELLector)

data(CELLector.PrimTum.BEMs)
data(CELLector.PathwayList)

tumours_BEM<-CELLector.PrimTum.BEMs$COREAD

CELLetor.mostSupported_CFEs(tumours_BEM,minlen = 2)

CELLector.Build_Search_Space(tumours = tumours_BEM,
                             pathwayFocused = "RAS-RAF-MEK-ERK / JNK signaling")


