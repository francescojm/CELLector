

## Exported Documented functions
CELLector.mostSupported_CFEs<-function(transactions,minSupport=0.05,minlen=1,maxLen=10){
  res<-
    eclat(transactions,
          parameter=list(supp=minSupport,minlen=minlen,maxlen=maxLen),control=list(verbose=F))

  if(length(res)>0){
    mostSI<-as(items(sort(res)[1]), "list")[[1]]
    if(length(mostSI)>1){
      absSupport<-sum(rowSums(transactions[,mostSI])==length(mostSI))
      SP<-names(which(rowSums(transactions[,mostSI])==length(mostSI)))
      support<-absSupport/nrow(transactions)
    }else{
      absSupport<-sum(transactions[,mostSI])
      SP<-rownames(transactions)[which(transactions[,mostSI]>0)]
      support<-absSupport/nrow(transactions)
    }
  }else{
    mostSI<-NULL
    SP<-NULL
    absSupport<-0
    support<-0
  }

  RES<-list(MSIS=mostSI,SUPPORT=support,absSUPPORT=absSupport,supportingSamples=SP)
  return(RES)
}
CELLector.cna_look_up <- function(cna_ID, cnaId_decode, TCGALabel) {

  cnaKEY16<-cnaId_decode
  CancerSpecificData <- cnaKEY16 %>% filter(CancerType == paste(TCGALabel))

  if (length(cna_ID) == 1) {

    info <- CancerSpecificData %>% filter(CNA_Identifier == paste0(cna_ID)) %>% select(Identifier, Recurrent, chr, start, stop,
                                                                                       locus, nGenes, ContainedGenes)


  } else if (length(cna_ID) > 1) {

    info <- CancerSpecificData %>% filter(CNA_Identifier %in% cna_ID) %>% select(Identifier, Recurrent, chr, start, stop, locus,
                                                                                 nGenes, ContainedGenes)

  }

  return(info)

}

CELLector.Build_Search_Space<-function(ctumours,
                                       cancerType,
                                       minlen=1,
                                       verbose=TRUE,
                                       mutOnly=FALSE,
                                       cnaOnly=FALSE,
                                       minGlobSupp=0.01,
                                       FeatureToExclude=NULL,
                                       pathway_CFEs = NULL,
                                       pathwayFocused=NULL,
                                       subCohortDefinition=NULL,
                                       NegativeDefinition=FALSE,
                                       cnaIdMap,
                                       cnaIdDecode,
                                       cdg){


  # rownames(ctumours)<-paste(rownames(ctumours),'_',1:nrow(ctumours),sep='')

  PANcna_KEY<-cnaIdMap
  cnaKEY16<-cnaIdDecode


  if(length(FeatureToExclude)>0){
    ctumours<-ctumours[,setdiff(colnames(ctumours),FeatureToExclude)]
  }

  if(length(subCohortDefinition)>0){
    if(is.element(subCohortDefinition,colnames(ctumours))){
      if(NegativeDefinition){
        ctumours<-ctumours[which(ctumours[,subCohortDefinition]==0),]
      }else{
        ctumours<-ctumours[which(ctumours[,subCohortDefinition]==1),]
      }
    }
  }


  if (length(pathwayFocused)>0){
    miniPathways<-pathway_CFEs
    events<-unique(unlist(miniPathways[pathwayFocused]))

    ii<-grep('cnaPANCA',events)

    cnaevents<-events[ii]

    cnaCS<-unique(unlist(c(PANcna_KEY[match(cnaevents,PANcna_KEY$Identifier),2:ncol(PANcna_KEY)])))
    cnaCS<-cnaCS[cnaCS!='']
    cnaCS<-grep(cancerType,cnaCS,value = TRUE)
    cnaCS<-unique(cnaKEY16$CNA_Identifier[match(cnaCS,cnaKEY16$Identifier)])
    cnaCS<-cnaCS[!is.na(cnaCS)]

    events<-c(events[-ii],cnaCS)

    ctumours<-ctumours[,intersect(colnames(ctumours),events)]
  }


  if(mutOnly){
    if(cnaOnly){
      stop("only one between mutOnly and cnaOnly can be TRUE", call. = FALSE)
    }
    ctumours<-ctumours[,setdiff(colnames(ctumours),grep('cna',colnames(ctumours),value=TRUE))]
  }
  if(cnaOnly){
    if(mutOnly){
      stop("only one between mutOnly and cnaOnly can be TRUE", call. = FALSE)
    }
    ctumours<-ctumours[,grep('cna',colnames(ctumours),value=TRUE)]
  }


  SysST<-NULL
  SysTREE<-NULL

  ROOT<-createNode(SystemStack = SysST,
                   transactions = ctumours,
                   currentPoints = rownames(ctumours),
                   currentFeatures = colnames(ctumours),
                   Type='root',
                   Parent.Idx=0,
                   maxId = 0,
                   ctype = cancerType,
                   minlen = minlen,
                   globalSupport = minGlobSupp,
                   cnaId_decode = cnaIdDecode)

  if(length(ROOT)>0){

    nROOT <- Node$new(paste(ROOT$Idx,paste(ROOT$decodedIS,collapse=', ')))

    SysST<-stackPush(SysST,ROOT)
    SysTREE<-stackPush(SysTREE,nROOT)

    NT<-addNodeToNavTable(NavTable = NULL,node = ROOT)

    if(verbose){
      print(paste('adding root node:',paste(ROOT$decodedIS,collapse=', ')))
    }

    MD<- -Inf

    while(length(SysST)>0){

      RES<-stackPop(SystemStack = SysST)
      nRES<-stackPop(SystemStack = SysTREE)

      currentNode<-RES$nNode
      currentNnode<-nRES$nNode

      SysST<-RES$SYST
      SysTREE<-nRES$SYST

      if(currentNode$Idx>MD){
        MD<-currentNode$Idx
      }


      RIGHTCHILD<-createNode(SystemStack = SysST,
                             transactions = ctumours,
                             currentPoints = setdiff(currentNode$currentPoints,currentNode$positivePoints),
                             currentFeatures = setdiff(currentNode$currentFeatures,currentNode$ItemSet),
                             Type='Right.Child',
                             Parent.Idx=currentNode$Idx,
                             maxId=MD,
                             ctype = cancerType,
                             minlen=minlen,
                             globalSupport = minGlobSupp,
                             cnaId_decode = cnaIdDecode)

      if(length(RIGHTCHILD)>0){

        if(verbose){
          print(paste('adding right child: ',paste(RIGHTCHILD$decodedIS,collapse=', '),
                      ' to node ',paste(currentNode$decodedIS,collapse=', ')))
        }

        SysST<-stackPush(SysST,RIGHTCHILD)

        RCnode<-currentNnode$AddChild(paste(RIGHTCHILD$Idx,
                                            paste(RIGHTCHILD$decodedIS,collapse=', ')))

        SysTREE<-stackPush(SysTREE,RCnode)

        NT<-addNodeToNavTable(NavTable = NT,node = RIGHTCHILD)
        MD<-MD+1
      }

      LEFTCHILD<-createNode(SystemStack = SysST,
                            transactions = ctumours,
                            currentPoints = currentNode$positivePoints,
                            currentFeatures = setdiff(currentNode$currentFeatures,currentNode$ItemSet),
                            Type='Left.Child',
                            Parent.Idx=currentNode$Idx,
                            maxId=MD,
                            ctype = cancerType,
                            minlen=minlen,
                            globalSupport = minGlobSupp,
                            cnaId_decode = cnaIdDecode)

      if(length(LEFTCHILD)>0){
        if(verbose){
          print(paste('adding left child: ',paste(LEFTCHILD$decodedIS,collapse=', '),
                      ' to node ',paste(currentNode$decodedIS,collapse=', ')))
        }
        SysST<-stackPush(SystemStack = SysST,node = LEFTCHILD)
        LFTnode<-currentNnode$AddChild(paste(LEFTCHILD$Idx,paste(LEFTCHILD$decodedIS,collapse=', ')))

        SysTREE<-stackPush(SysTREE,LFTnode)

        NT<-addNodeToNavTable(NavTable = NT,node = LEFTCHILD)
      }
    }

    tmp<-unlist(str_split(Get(Traverse(nROOT,'pre-order'),'name'),' '))

    suppressWarnings(pre_orderVisit<-as.numeric(tmp)[!is.na(as.numeric(tmp))])

    globalSuppAttrb<-round(100*NT$GlobalSupport[match(pre_orderVisit,NT$Idx)],digits = 2)
    nodeTypeAttrb<-as.character(NT$Type[match(pre_orderVisit,NT$Idx)])

    Set(Traverse(nROOT,traversal = 'pre-order'),NodeType=nodeTypeAttrb)
    Set(Traverse(nROOT,traversal = 'pre-order'),GlobalSupp=globalSuppAttrb)
  }else{
    NT<-matrix(1)
    nROOT<-NULL
  }

  return(list(navTable=NT,TreeRoot=nROOT))
}

CELLector.unicizeSamples<-function(ctumours,keepReplicates=TRUE){
  if (keepReplicates){
    colnames(ctumours)<-paste(colnames(ctumours),'_',1:ncol(ctumours),sep='')
  }else{
    ctumours<-ctumours[,unique(colnames(ctumours))]
  }
  return(ctumours)
}
CELLector.createAllSignatures<-function(NavTab){
  NN<-NavTab$Idx

  signatures<-vector()
  encodedsignatures<-vector()
  subTypeSizes<-vector()
  for (i in 1:length(NN)){
    S<-createRuleFromNode(NavTab,NN[i])
    signatures[i]<-S$S
    encodedsignatures[i]<-S$ES
    subTypeSizes[i]<-NavTab$GlobalSupport[i]
  }

  names(signatures)<-NN
  names(encodedsignatures)<-NN
  return(list(S=signatures,ES=encodedsignatures,STS=100*subTypeSizes))
}
CELLector.solveFormula<-function(RULE,dataset,To_beExcluded=NULL){

  r<-dataset[,2]
  COSMICids<-dataset[,1]
  CELLlineData<-as.matrix(dataset[,3:ncol(dataset)])
  rownames(CELLlineData)<-r

  tdataset<-CELLlineData

  tdataset<-tdataset[setdiff(rownames(tdataset),To_beExcluded),]

  tokenize<-unlist(str_split(RULE,' '))
  tokenize<-tokenize[tokenize!='']

  NegVar<-grep('~',tokenize)
  PosVar<-setdiff(1:length(tokenize),NegVar)
  ortok<-tokenize
  tokenize<-str_replace(tokenize,'~','')

  tdataset<-t(tdataset)

  notPresentPosVar<-setdiff(tokenize[PosVar],rownames(tdataset))
  notPresentNegVar<-setdiff(tokenize[NegVar],rownames(tdataset))

  if(length(notPresentNegVar)){

    toAdd<-matrix(0,length(notPresentNegVar),ncol(tdataset),dimnames = list(notPresentNegVar,colnames(tdataset)))
    tdataset<-rbind(tdataset,toAdd)

  }

  if(length(notPresentPosVar)==0){
    tdataset<-rbind(tdataset[tokenize[PosVar],],1-tdataset[tokenize[NegVar],])
    rownames(tdataset)<-c(ortok[PosVar],ortok[NegVar])

    positiveSamples<-names(which(colSums(tdataset)==length(ortok)))
    nsamples<-length(positiveSamples)
    frac<-nsamples/nrow(dataset)

    return(list(PS=positiveSamples,N=nsamples,PERC=frac))
  }else{
    return(NULL)
  }

}
CELLector.buildModelMatrix<-function(Sigs,dataset,searchSpace){

  encodedSignatures<-Sigs$ES

  ordataset<-dataset
  ### takes just the numerical part of the CELLlineData
  r<-dataset[,2]
  COSMICids<-dataset[,1]
  dataset<-dataset[,3:ncol(dataset)]
  rownames(dataset)<-r

  ### map cell lines onto subtypes, based on the signatures
  MODELS<-vector()
  for (cc in 1:length(encodedSignatures)){
    solved<-CELLector.solveFormula(encodedSignatures[[cc]],dataset = ordataset)
    MODELS[cc]<-paste(sort(solved$PS),collapse=', ')
  }

  ### visit the searching space for the selection
  visit<-CELLector.selectionVisit(searchSpace)

  ### put the cell lines in the same order in which the corresponding subtypes are
  ### encountered in the visit of the searching space
  sortedModels<-MODELS[visit]

  ### ignore subtypes with not mached models (but actually this is what is most interesting for you, right?)
  NODEidx<-visit[which(sortedModels!='')]
  sortedModels<-sortedModels[which(sortedModels!='')]

  modellist<-sortedModels

  cls_<-list()
  for (i in 1:length(modellist)){

    cls_[[i]]<-unlist(str_split(modellist[i],', '))
  }

  mappedCLS<-sort(unique(unlist(cls_)))

  modelMatrix<-matrix(0,length(modellist),length(mappedCLS),dimnames = list(1:length(modellist),mappedCLS))

  for (i in 1:length(modellist)){
    modelMatrix[i,cls_[[i]]]<-1
  }

  return(modelMatrix)
}
CELLector.makeSelection<-function(modelMat,n,searchSpace){


  selectedCLS<-vector()
  modelAccounted<-vector()

  modelMatIdx<-0

  flag<-1

  TOinclude<-colnames(modelMat)

  while(length(selectedCLS)<n & length(TOinclude)>0){

    modelMatIdx<-modelMatIdx+1
    if(modelMatIdx>nrow(modelMat)){
      modelMatIdx<-1
    }

    if(length(TOinclude)>1){
      possibleSelections<-names(which(modelMat[modelMatIdx,setdiff(colnames(modelMat),selectedCLS)]>0))
    }else{
      possibleSelections<-TOinclude
    }

    if (length(possibleSelections)>0){
      if (modelMatIdx<nrow(modelMat) & length(possibleSelections)>1){
        remaining<-modelMat[(modelMatIdx+1):nrow(modelMat),possibleSelections]


        if (is.matrix(remaining)){
          remaining<-colSums(remaining)
        }

        MINC<-NULL
        if(length(remaining)){
          MINC<-names(which(remaining==min(remaining)))
        }

        if(length(MINC)>0){
          selection<-sample(MINC,1)
        }else{
          selection<-sample(possibleSelections,1)
        }
      }else{
        selection<-sample(possibleSelections,1)
      }

      selectedCLS[flag]<-selection
      modelAccounted[flag]<-modelMatIdx
      flag<-flag+1

      TOinclude<-setdiff(TOinclude,selectedCLS)


    }


  }


  signatures<-CELLector.createAllSignatures(searchSpace)

  RES<-data.frame('Tumour SubType Index'=modelAccounted,
                  'Representative Cell Line'=selectedCLS,
                  'Signature'=signatures$S[modelAccounted],
                  'percentage patients'=signatures$STS[modelAccounted],
                  stringsAsFactors = FALSE)
  return(RES)
}

## Other Exported functions
CELLector.selectionVisit<-function(TAV){
  reducedTab<-TAV[,c(1,4,5,10,11)]
  currentNode<-1
  pileIdx<-1
  pile<-currentNode
  nodeType<-reducedTab[currentNode,2]

  while(pileIdx<=length(pile)){

    #print(pile)

    pile<-c(pile,rightMostPath(reducedTab,pile[pileIdx]))
    nodeType<-reducedTab[pile,2]
    #print(pile[pileIdx:length(pile)])

    pile<-c(pile,setdiff(leftChildPattern(reducedTab,pile[(pileIdx):length(pile)]),pile))
    nodeType<-reducedTab[pile,2]

    # print(pile[pileIdx:length(pile)])

    dd <- which(nodeType=='Left.Child')
    pileIdx<-dd[dd>pileIdx][1]

    if(is.na(pileIdx)){
      break
    }
  }
  return(pile)
}


## not Exported functions

createRuleFromNode<-function(NavTab,nodeIdx){

  RULES<-list()

  pos<-match(nodeIdx,NavTab$Idx)
  orpos<-pos

  SIGNATURE<-''
  encodedSIGNATURE<-''

  while(pos>0){

    currentType<-NavTab$Type[pos]

    pos<-NavTab$Parent.Idx[pos]

    if(pos>0){
      if (currentType=='Right.Child'){
        prefix<-'~'
      }else{
        prefix<-''
      }
      currentTerm<-paste(prefix,NavTab$ItemsDecoded[pos],sep='')
      SIGNATURE<-paste(currentTerm,SIGNATURE)
      EcurrentTerm<-paste(prefix,NavTab$Items[pos],sep='')
      encodedSIGNATURE<-paste(EcurrentTerm,encodedSIGNATURE)
    }

  }

  SIGNATURE<-paste(SIGNATURE,NavTab$ItemsDecoded[orpos],sep='')
  encodedSIGNATURE<-paste(encodedSIGNATURE,NavTab$Items[orpos],sep='')

  SIGNATURE<-str_trim(SIGNATURE)
  encodedSIGNATURE<-str_trim(encodedSIGNATURE)
  return(list(S=SIGNATURE,ES=encodedSIGNATURE))
}
createNode<-function(SystemStack,
                     transactions,
                     currentPoints,
                     currentFeatures,
                     Index,
                     Type,
                     Parent.Idx,maxId,
                     globalSupport=0.02,
                     minlen=1,ctype,
                     cnaId_decode,
                     cdg=NULL){

  nodeIdx<-maxId+1

  if(length(currentFeatures)==1 | length(currentPoints)==1){
    currentDataset<-matrix(transactions[currentPoints,currentFeatures],
                           length(currentPoints),length(currentFeatures),
                           dimnames = list(currentPoints,currentFeatures))
  }else{
    currentDataset<-transactions[currentPoints,currentFeatures]
  }


  nsamples<-ceiling(globalSupport*nrow(transactions))

  minSupport<-nsamples/nrow(currentDataset)

  if(minSupport<1 & sum(unlist(c(currentDataset)))>0){

    RES<-CELLector.mostSupported_CFEs(transactions = currentDataset,
                              minSupport = minSupport,
                              minlen = minlen)

    if (length(RES$MSIS)==0){
      return(NULL)
    }else{


      IS<-RES$MSIS

      dIS<-decodeCNAs(ctype = ctype,codedCNA = IS, cnaId_decode = cnaId_decode, cdg = cdg)

      if(length(IS)==1){
        gs<-sum(transactions[,IS])/nrow(transactions)
      }else{
        gs<-sum(rowSums(transactions[,IS])==length(IS))/nrow(transactions)
      }

      nNODE<-
        list(Idx=nodeIdx,
             currentPoints=currentPoints,
             positivePoints=RES$supportingSamples,
             currentFeatures=currentFeatures,
             ItemSet=RES$MSIS,
             decodedIS=dIS,
             Type=Type,
             Parent.Idx=Parent.Idx,
             AbsSupport=RES$absSUPPORT,
             CurrentTotal=length(currentPoints),
             PercSupport=RES$SUPPORT,
             GlobalSupport=RES$absSUPPORT/nrow(transactions))

      return(nNODE)
    }
  }else{
    return(NULL)
  }
}

decodeCNAs<-function(ctype,codedCNA,cnaId_decode,cdg){
  IS<-codedCNA
  icna<-grep('cna',IS)

  #if(length(icna)>0){

  noncna<-setdiff(1:length(IS),icna)


  if (length(icna)>0){


    cnalu<-CELLector.cna_look_up(IS[icna],TCGALabel = ctype,cnaId_decode = cnaId_decode)

    altType<-cnalu$Recurrent
    AT<-rep(NA,length(altType))
    AT[which(as.character(altType)=='Amplification')]<-'G'
    AT[which(as.character(altType)=='Deletion')]<-'L'

    loci<-as.character(cnalu$locus)
    genes<-as.character(cnalu$ContainedGenes)

    for (i in 1:length(genes)){
      cgenes<-unlist(str_split(genes[i],','))
      if (length(cgenes)>5){
        cgenes<-intersect(cgenes,cdg)
        cgenes<-paste(paste(cgenes,collapse=','),'...',sep='')
      }else{
        cgenes<-genes[i]
      }
      genes[i]<-cgenes
    }


    IS[icna]<-paste(AT,loci,'(',genes,')',sep='')}

  IS[noncna]<-paste(IS[noncna],'mut',sep='')
  #}
  return(IS)
}

stackPush<-function(SystemStack,node){
  currentLength<-length(SystemStack)
  if(currentLength>0){
    SystemStack[currentLength+1]<-list(node)
  }else{
    SystemStack<-list(node)
  }
  return(SystemStack)
}
stackTop<-function(SystemStack){
  return(SystemStack[[length(SystemStack)]])
}
stackPop<-function(SystemStack){

  if(length(SystemStack)>0){
    nNode<-SystemStack[[length(SystemStack)]]

    if(length(SystemStack)>1){
      SystemStack<-SystemStack[1:(length(SystemStack)-1)]
    }else{
      SystemStack<-NULL
    }
  }

  return(list(nNode=nNode,SYST=SystemStack))
}

addNodeToNavTable<-function(NavTable,node){
    if(length(NavTable)==0){
      NavTable<-data.frame(c(node[1],paste(node[[5]],collapse=', '),paste(node[[6]],collapse=', '),node[7:12],
                             0,0,paste(node[[2]],collapse=','),paste(node[[4]],collapse=','),paste(node[[3]],collapse=',')))
      colnames(NavTable)[2]<-'Items'
      colnames(NavTable)[3]<-'ItemsDecoded'
      colnames(NavTable)[10]<-'Left.Child.Index'
      colnames(NavTable)[11]<-'Right.Child.Index'
      colnames(NavTable)[12]<-'currentPoints'
      colnames(NavTable)[13]<-'currentFeatures'
      colnames(NavTable)[14]<-'positivePoints'

    }else{
      newChunk<-data.frame(c(node[1],paste(node[[5]],collapse=', '),paste(node[[6]],collapse=', '),node[7:12],
                             0,0,paste(node[[2]],collapse=','),paste(node[[4]],collapse=','),paste(node[[3]],collapse=',')))

      colnames(newChunk)[2]<-'Items'
      colnames(newChunk)[3]<-'ItemsDecoded'
      colnames(newChunk)[10]<-'Left.Child.Index'
      colnames(newChunk)[11]<-'Right.Child.Index'
      colnames(newChunk)[12]<-'currentPoints'
      colnames(newChunk)[13]<-'currentFeatures'
      colnames(newChunk)[14]<-'positivePoints'

      NavTable<-rbind(NavTable,newChunk)

      if(node$Type=='Left.Child'){
        NavTable$Left.Child.Index[which(NavTable$Idx==node$Parent.Idx)]<-node$Idx
      }

      if(node$Type=='Right.Child'){
        NavTable$Right.Child.Index[which(NavTable$Idx==node$Parent.Idx)]<-node$Idx
      }
    }

    return(NavTable)
  }

rightMostPath<-function(Tab,node){

  currentNode<-node

  rPath<-vector()

  flag<-1
  while(Tab[currentNode,'Right.Child.Index']>0){
    rPath[flag]<-Tab[currentNode,'Right.Child.Index']
    currentNode<-rPath[flag]
    flag<-flag+1
  }

  return(rPath)
}

leftChildPattern<-function(Tab,nodePattern){
  lc<-Tab[nodePattern,'Left.Child.Index']
  lc<-lc[lc>0]
  return(lc)
}


# library(igraph)
# library(BBmisc)
# library(dplyr)
# library(stringr)
# library(collapsibleTree)
# library(data.tree)
# library(plotly)
#
# source('R/CELLector.cna_look_up.R')
#
# cdg<-
#   union(union(cancerDrivers$ActingDriver_Symbol,InigoList),'MYC')
#
#
# ## internal functions
#
# buildModelMatrix<-function(modellist){
#
#
#   cls_<-list()
#   for (i in 1:length(modellist)){
#
#     cls_[[i]]<-unlist(str_split(modellist[i],', '))
#   }
#
#   mappedCLS<-sort(unique(unlist(cls_)))
#
#   modelMatrix<-matrix(0,length(modellist),length(mappedCLS),dimnames = list(1:length(modellist),mappedCLS))
#
#   for (i in 1:length(modellist)){
#     modelMatrix[i,cls_[[i]]]<-1
#   }
#
#   return(modelMatrix)
# }
#
# makeSelection<-function(modelMat,n){
#
#
#   selectedCLS<-vector()
#   modelAccounted<-vector()
#
#   modelMatIdx<-0
#
#   flag<-1
#
#   TOinclude<-colnames(modelMat)
#
#   while(length(selectedCLS)<n & length(TOinclude)>0){
#
#     modelMatIdx<-modelMatIdx+1
#     if(modelMatIdx>nrow(modelMat)){
#       modelMatIdx<-1
#     }
#
#     if(length(TOinclude)>1){
#       possibleSelections<-names(which(modelMat[modelMatIdx,setdiff(colnames(modelMat),selectedCLS)]>0))
#     }else{
#       possibleSelections<-TOinclude
#     }
#
#     if (length(possibleSelections)>0){
#       if (modelMatIdx<nrow(modelMat) & length(possibleSelections)>1){
#         remaining<-modelMat[(modelMatIdx+1):nrow(modelMat),possibleSelections]
#
#
#         if (is.matrix(remaining)){
#           remaining<-colSums(remaining)
#         }
#
#         MINC<-NULL
#         if(length(remaining)){
#           MINC<-names(which(remaining==min(remaining)))
#         }
#
#         if(length(MINC)>0){
#           selection<-sample(MINC,1)
#         }else{
#           selection<-sample(possibleSelections,1)
#         }
#       }else{
#         selection<-sample(possibleSelections,1)
#       }
#
#       selectedCLS[flag]<-selection
#       modelAccounted[flag]<-modelMatIdx
#       flag<-flag+1
#
#       TOinclude<-setdiff(TOinclude,selectedCLS)
#
#       print(selectedCLS)
#       print(modelAccounted)
#       print(length(selectedCLS))
#
#     }
#
#
#   }
#
#
#   RES<-data.frame(modelAccounted,selectedCLS,stringsAsFactors = FALSE)
#
# }
#
# getPositiveSamples<-function(transactions,itemset){
#   if(length(itemset)==1){
#     return(which(transactions[,itemset]>0))
#   }else{
#     return(which(rowSums(transactions[,itemset])==length(itemset)))
#   }
# }

#

#

#
# selectionVisit<-function(TAV){
#   reducedTab<-TAV[,c(1,4,5,10,11)]
#   currentNode<-1
#   pileIdx<-1
#   pile<-currentNode
#   nodeType<-reducedTab[currentNode,2]
#
#   while(pileIdx<=length(pile)){
#
#     #print(pile)
#
#     pile<-c(pile,rightMostPath(reducedTab,pile[pileIdx]))
#     nodeType<-reducedTab[pile,2]
#     #print(pile[pileIdx:length(pile)])
#
#     pile<-c(pile,setdiff(leftChildPattern(reducedTab,pile[(pileIdx):length(pile)]),pile))
#     nodeType<-reducedTab[pile,2]
#
#     # print(pile[pileIdx:length(pile)])
#
#     dd <- which(nodeType=='Left.Child')
#     pileIdx<-dd[dd>pileIdx][1]
#
#     if(is.na(pileIdx)){
#       break
#     }
#   }
#   return(pile)
# }
# rightMostPath<-function(Tab,node){
#
#   currentNode<-node
#
#   rPath<-vector()
#
#   flag<-1
#   while(Tab[currentNode,'Right.Child.Index']>0){
#     rPath[flag]<-Tab[currentNode,'Right.Child.Index']
#     currentNode<-rPath[flag]
#     flag<-flag+1
#   }
#
#   return(rPath)
# }
# leftChildPattern<-function(Tab,nodePattern){
#   lc<-Tab[nodePattern,'Left.Child.Index']
#   lc<-lc[lc>0]
#   return(lc)
# }
#
# ## exported functions
#
#
#
# createRuleFromNode<-function(NavTab,nodeIdx){
#
#   RULES<-list()
#
#   pos<-match(nodeIdx,NavTab$Idx)
#   orpos<-pos
#
#   SIGNATURE<-''
#   encodedSIGNATURE<-''
#
#   while(pos>0){
#
#     currentType<-NavTab$Type[pos]
#
#     pos<-NavTab$Parent.Idx[pos]
#
#     if(pos>0){
#       if (currentType=='Right.Child'){
#         prefix<-'~'
#       }else{
#         prefix<-''
#       }
#       currentTerm<-paste(prefix,NavTab$ItemsDecoded[pos],sep='')
#       SIGNATURE<-paste(currentTerm,SIGNATURE)
#       EcurrentTerm<-paste(prefix,NavTab$Items[pos],sep='')
#       encodedSIGNATURE<-paste(EcurrentTerm,encodedSIGNATURE)
#     }
#
#   }
#
#   SIGNATURE<-paste(SIGNATURE,NavTab$ItemsDecoded[orpos],sep='')
#   encodedSIGNATURE<-paste(encodedSIGNATURE,NavTab$Items[orpos],sep='')
#
#   SIGNATURE<-str_trim(SIGNATURE)
#   encodedSIGNATURE<-str_trim(encodedSIGNATURE)
#   return(list(S=SIGNATURE,ES=encodedSIGNATURE))
# }
# createAllSignatures<-function(NavTab){
#   NN<-NavTab$Idx
#
#   signatures<-vector()
#   encodedsignatures<-vector()
#
#   for (i in 1:length(NN)){
#     S<-createRuleFromNode(NavTab,NN[i])
#     signatures[i]<-S$S
#     encodedsignatures[i]<-S$ES
#
#   }
#   names(signatures)<-NN
#   names(encodedsignatures)<-NN
#   return(list(S=signatures,ES=encodedsignatures))
#
# }
#
# solveFormula<-function(RULE,dataset){
#
#   tokenize<-unlist(str_split(RULE,' '))
#   tokenize<-tokenize[tokenize!='']
#
#   NegVar<-grep('~',tokenize)
#   PosVar<-setdiff(1:length(tokenize),NegVar)
#   ortok<-tokenize
#   tokenize<-str_replace(tokenize,'~','')
#
#
#
#   tdataset<-t(dataset)
#
#   notPresentPosVar<-setdiff(tokenize[PosVar],rownames(tdataset))
#   notPresentNegVar<-setdiff(tokenize[NegVar],rownames(tdataset))
#
#   if(length(notPresentNegVar)){
#
#     toAdd<-matrix(0,length(notPresentNegVar),ncol(tdataset),dimnames = list(notPresentNegVar,colnames(tdataset)))
#     tdataset<-rbind(tdataset,toAdd)
#
#   }
#
#   if(length(notPresentPosVar)==0){
#     tdataset<-rbind(tdataset[tokenize[PosVar],],1-tdataset[tokenize[NegVar],])
#     rownames(tdataset)<-c(ortok[PosVar],ortok[NegVar])
#
#     positiveSamples<-names(which(colSums(tdataset)==length(ortok)))
#     nsamples<-length(positiveSamples)
#     frac<-nsamples/nrow(dataset)
#
#     return(list(PS=positiveSamples,N=nsamples,PERC=frac))
#   }else{
#     return(NULL)
#   }
#
# }
#
# complementarPieChart<-function(Tree,NavTab,nodeIdx){
#
#   supports<-vector()
#   supports[1]<-NavTab$GlobalSupport[[nodeIdx]]
#
#   flag<-2
#   names(supports)<-paste('SubT.',nodeIdx,sep='')
#
#   NIDX<-nodeIdx
#
#   while(NavTab$Right.Child.Index[nodeIdx]>0){
#
#     nodeIdx<-NavTab$Right.Child.Index[nodeIdx]
#     supports[flag]<-NavTab$GlobalSupport[[nodeIdx]]
#     names(supports)[flag]<-paste('SubT.',nodeIdx,sep='')
#     flag<-flag+1
#     NIDX<-c(NIDX,nodeIdx)
#
#   }
#
#   tmpCol<-Get(Traverse(Tree,traversal = 'level'),'Colors')
#   nn<-names(tmpCol)
#   nn<-str_split(nn,' ')
#   nn<-as.numeric(unlist(lapply(nn,function(x){x[[1]][1]})))
#
#   id<-match(NIDX,nn)
#   COLORS<-tmpCol[id]
#
#   supports<-c(100*supports,100-100*sum(supports))
#   names(supports)[flag]<-'Others'
#
#   return(list(supports=supports,COLORS=COLORS))
# }
#
#
# # NT<-CELLector_search_space(cancerType = 'COREAD',
# #                             verbose = TRUE,
# #                             mutOnly = FALSE,
# #                             cnaOnly=FALSE,
# #                             FeatureToExclude = NULL,
# #                             minlen = 1,
# #                             minGlobSupp = 0.01,
# #                             pathwayFocused = c("RAS-RAF-MEK-ERK / JNK signaling"))
#
# # COLORSbyLev<-rep('gray',NT$TreeRoot$totalCount)
# # COLORSbyLev[which(Get(Traverse(NT$TreeRoot,traversal = 'level'),attribute = 'NodeType')=='Right.Child')]<-'orange'
# # COLORSbyLev[which(Get(Traverse(NT$TreeRoot,traversal = 'level'),attribute = 'NodeType')=='Left.Child')]<-'purple'
# #
# # collapsibleTree(NT$TreeRoot,fill = COLORSbyLev)
# #
# # save(NT,file='appData/ct_example.RData')
