



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
CELLector.hms_look_up <- function(hms_ID, hmsId_decode, TCGALabel) {

  cnaKEY16<-hmsId_decode
  CancerSpecificData <- cnaKEY16 %>% filter(Cancer.Types == paste(TCGALabel))

  if (length(hms_ID) == 1) {

    info <- CancerSpecificData %>% filter(hms_id == paste0(hms_ID))

  } else if (length(hms_ID) > 1) {

    info <- CancerSpecificData %>% filter(hms_id %in% hms_ID)
  }

  return(info)

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
    S<-createRuleFromNode(NavTab = NavTab,nodeIdx = NN[i])
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


  #dRULE<-str_replace_all(RULE,', ','-X-X-X-')

  dRULE<-RULE
  tokenize<-unlist(str_split(dRULE,', '))
  #tokenize<-tokenize[tokenize!='']
  ortok<-tokenize
  #ortok<-str_replace_all(ortok,'-X-X-X-',', ')

  #Id_of_multipleVar<-grep('-X-X-X-',tokenize)

  NegVar<-grep('~',tokenize)
  PosVar<-setdiff(1:length(tokenize),NegVar)

  #NegVarMultiple<-intersect(NegVar,Id_of_multipleVar)
  #PosVarMultiple<-intersect(PosVar,Id_of_multipleVar)

  #NegVarMultiple<-NULLintersect(NegVar,Id_of_multipleVar)
  #PosVarMultiple<-NULLintersect(PosVar,Id_of_multipleVar)

  tokenize<-str_replace_all(tokenize,'~','')

  #NegVarIndividual<-tokenize[setdiff(NegVar,NegVarMultiple)]
  #PosVarIndividual<-tokenize[setdiff(PosVar,PosVarMultiple)]
  #
  # if(length(PosVarMultiple)>0){
  #   currentMultiple<-tokenize[PosVarMultiple]
  #   individualMultiple<-unlist(str_split(tokenize[PosVarMultiple],'-X-X-X-'))
  #   NewTokenizePos<-setdiff(c(setdiff(tokenize,tokenize[PosVarMultiple]),individualMultiple),NegVarIndividual)
  # }else{
  #   NewTokenizePos<-NULL
  # }
  #
  # if(length(NegVarMultiple)>0){
  #   currentMultiple<-tokenize[NegVarMultiple]
  #   individualMultiple<-unlist(str_split(tokenize[NegVarMultiple],'-X-X-X-'))
  #   NewTokenizeNeg<-setdiff(c(setdiff(tokenize,tokenize[NegVarMultiple]),individualMultiple),PosVarIndividual)
  # }else{
  #   NewTokenizeNeg<-NULL
  # }

  # if(length(grep('-X-X-X-',tokenize))>0){
  #   tokenize<-tokenize[-grep('-X-X-X-',tokenize)]
  # }
  #
  # if(length(grep('-X-X-X-',NewTokenizePos))>0){
  #   NewTokenizePos<-NewTokenizePos[-grep('-X-X-X-',NewTokenizePos)]
  # }
  #
  # if(length(grep('-X-X-X-',NewTokenizeNeg))>0){
  #   NewTokenizeNeg<-NewTokenizeNeg[-grep('-X-X-X-',NewTokenizeNeg)]
  # }
  #
  # newNegVar<-sort(union(setdiff(NewTokenizeNeg,NewTokenizePos),NegVarIndividual))
  # newPosVar<-sort(union(setdiff(NewTokenizePos,NewTokenizeNeg),PosVarIndividual))
  #
  # tokenize<-union(newNegVar,newPosVar)
  #
  # tokenize<-sort(tokenize)
  #
  #
  # PosVar<-match(newPosVar,tokenize)
  # NegVar<-match(newNegVar,tokenize)

  tdataset<-t(tdataset)

  notPresentPosVar<-setdiff(tokenize[PosVar],rownames(tdataset))
  notPresentNegVar<-setdiff(tokenize[NegVar],rownames(tdataset))

  if(length(notPresentNegVar)>0){
    toAdd<-matrix(0,length(notPresentNegVar),ncol(tdataset),dimnames = list(notPresentNegVar,colnames(tdataset)))
    tdataset<-rbind(tdataset,toAdd)
  }

  if(length(notPresentPosVar)==0){
    tdataset<-rbind(tdataset[tokenize[PosVar],],1-tdataset[tokenize[NegVar],])

    positiveSamples<-names(which(colSums(tdataset)==length(tokenize)))
    nsamples<-length(positiveSamples)
    frac<-nsamples/nrow(dataset)

    return(list(PS=positiveSamples,N=nsamples,PERC=frac))
  }else{
    return(list(PS=NULL,N=0,PERC=0))
  }
}
CELLector.buildModelMatrix<-function(Sigs,dataset,searchSpace){

  encodedSignatures<-Sigs

  ordataset<-dataset
  ### takes just the numerical part of the CELLlineData
  r<-dataset[,2]
  COSMICids<-dataset[,1]
  dataset<-dataset[,3:ncol(dataset)]
  rownames(dataset)<-r

  ### map cell lines onto subtypes, based on the signatures
  MODELS<-vector()
  for (cc in 1:length(encodedSignatures)){
    solved<-CELLector.solveFormula(RULE = encodedSignatures[[cc]],dataset = ordataset)
    suppressWarnings(MODELS[cc]<-paste(sort(solved$PS),collapse=', '))
  }

  ### visit the searching space for the selection
  suppressWarnings(visit<-CELLector.selectionVisit(searchSpace))

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

  rownames(modelMatrix)<-as.character(NODEidx)
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
CELLector.visualiseSearchingSpace<-function(searchSpace,CLdata=NULL){

  RelatesToFatherAs <- rep('-',searchSpace$TreeRoot$totalCount)
  RelatesToFatherAs[which(Get(Traverse(searchSpace$TreeRoot,traversal = 'level'),
                              attribute = 'NodeType')=='Right.Child')]<-'Complement'
  RelatesToFatherAs[which(Get(Traverse(searchSpace$TreeRoot,traversal = 'level'),
                              attribute = 'NodeType')=='Left.Child')]<-'Refinement'

  searchSpace$TreeRoot$Set(RelatesToFatherAs=RelatesToFatherAs,traversal = 'level')


  levelVisitOrder<-as.numeric(unlist(lapply(str_split(Get(Traverse(searchSpace$TreeRoot,
                                                                   traversal = 'level'),'name'),' '),function(x){x[1]})))


  NPs<-createHtmlNodeProperties(LocalSearchSpace = searchSpace,
                                CLdataset = CLdata)

  searchSpace$TreeRoot$Set(size=searchSpace$navTable$GlobalSupport[levelVisitOrder],traversal = 'level')
  #searchSpace$TreeRoot$Set(tthm=NPs[levelVisitOrder],traversal='level')

  searchSpace$TreeRoot$Set(tthm=NPs,traversal='level')


  collapsibleTree(searchSpace$TreeRoot,
                  nodeSize = 'size',
                  fill = 'Colors',
                  inputId = 'searchSpace',
                  tooltip = TRUE,
                  tooltipHtml = 'tthm',
                  attribute = 'RelatesToFatherAs')

}
CELLector.visualiseSearchingSpace_sunBurst<-function(searchSpace){

  SBF<-sunBurstFormat(searchSpace = searchSpace)

  sequences<-SBF
  tmpCol <- Get(Traverse(searchSpace$TreeRoot,traversal = 'level'),'Colors')
  ttmp<-tmpCol

  names(ttmp)<-NULL

  nvoid<-grep('Others',unique(unlist(strsplit(sequences$V1,'-'))),value = TRUE)

  stpes<-nvoid

  colors <- list(
    domain=c('0 TOTAL',names(tmpCol),stpes),
    range=c('black',ttmp,rep('white',length(stpes)))
  )

  names(ttmp)<-NULL

  nvoid<-grep('Others',unique(unlist(strsplit(SBF$V1,'-'))),value = TRUE)

  stpes<-nvoid

  colors <- list(
    domain=c('0 TOTAL',names(tmpCol),stpes),
    range=c('white',ttmp,rep('white',length(stpes)))
  )



  #sunburst(SBF,percent = TRUE,count = FALSE,colors=colors)

  htmlwidgets::onRender(
    sunburst(SBF,percent = FALSE,count = FALSE,colors=colors,

             explanation = "function(d) {     var ssr = d.data.name
             if (!ssr.match(/Others/gi)){
             return ssr
             }
}"),
        "
    function(el,x){
    d3.select(el).select('.sunburst-sidebar').remove()
    }
    "
    )
  }
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
CELLector.changeSScolors<-function(searchSpace){
  CC <- colors(distinct = TRUE)
  CC <- CC[setdiff(1:length(CC),c(grep('gray',CC),'black'))]
  CC <- rgb(t(col2rgb(CC)),maxColorValue = 255)

  COLORSbyLev <- CC[sample(length(CC))][1:searchSpace$TreeRoot$totalCount]

  names(COLORSbyLev)<-names(Get(Traverse(searchSpace$TreeRoot),'Names'))
  searchSpace$TreeRoot$Set(Colors=COLORSbyLev,traversal = 'level')

  treeLabels<-unlist(lapply(str_split(Get(Traverse(searchSpace$TreeRoot,'level'),'name'),'[(]'),
                            function(x){x[1]}))

  nodeIdx<-as.numeric(unlist(lapply(str_split(treeLabels,' '),function(x){x[1]})))

  COLORS<-rep(NA,length(nodeIdx))
  COLORS[nodeIdx]<-COLORSbyLev

  searchSpace$COLORS<-COLORS

  return(searchSpace)
}
CELLector.Score <- function(NavTab, CELLlineData,alfa=0.75){

  if(alfa>=0 & alfa<=1){
    beta<-1-alfa

    Signatures <- CELLector.createAllSignatures(NavTab)

    SignaturesES<-Signatures$S

    MODELS<-vector()
    for(cc in 1:length(Signatures$ES)){
      solved<-CELLector.solveFormula(Signatures$ES[[cc]],dataset = CELLlineData)
      MODELS[cc]<-paste(sort(solved$PS),collapse=', ')
    }

    tabCSS <- cbind(NavTab[ ,c(1,9)], SignaturesES, MODELS)

    SortedSubpop <- CELLector.selectionVisit(TAV = NavTab)
    tabCSS <- tabCSS[tabCSS$Idx[SortedSubpop],]

    vecSigLength <- c()
    for(i in 1:length(SignaturesES)){
      x <- length(unlist(str_split(SignaturesES[i], " ")))
      vecSigLength[i] <- x
    }
    longestSigniture <- max(vecSigLength)

    vecScores <- c()
    sigLength <- c()


    for(g in 1:nrow(tabCSS)){

      signiture_length <- length(unlist(str_split(tabCSS[g,3], " ")))
      sigLength[g] <- signiture_length

      normalizedSigLenght <- signiture_length/longestSigniture

      CELLectorScore <- alfa*normalizedSigLenght+beta*tabCSS[g,2]

      CELLectorScore <- alfa*normalizedSigLenght+beta*tabCSS[g,2]

      vecScores[g] <- CELLectorScore

    }

    tabCSS$CELLectorScores <- vecScores
    tabCSS$SignitureLength <- sigLength

    # Remove Subtypes with no representative cell lines
    sub_tabCSS <- tabCSS %>% filter(MODELS!="")

    # Deconvolute sub_tabCSS - create table with a cell line per row
    Score_perCL <- data.frame(matrix(nrow=0, ncol=5))
    colnames(Score_perCL) <- c("RepCellLines", "Idx GlobalSupport", "SignaturesES", "CELLectorScores", "SignitureLength")

    for(i in 1:nrow(sub_tabCSS)){

      x <- unlist(str_split(sub_tabCSS[i,4], " "))
      y <- sub_tabCSS[i,c(1,2,3,5,6)]
      RepCellLines <- gsub("\\,*", "", x)
      df <- cbind(data.frame(RepCellLines),y, row.names = NULL)
      Score_perCL <- rbind(Score_perCL, df)
    }

    ucL<-unique(Score_perCL$RepCellLines)

    Scores<-do.call(rbind,
                    lapply(ucL,function(x){
                      id<-which(Score_perCL$RepCellLines==x)
                      Score_perCL[id[order(Score_perCL$CELLectorScores[id],decreasing = TRUE)[1]],]
                    }))

    Scores<-Scores[order(Scores$CELLectorScores,decreasing = TRUE),]
    rownames(Scores)<-NULL

    Scores<-data.frame(CellLines=as.character(Scores$RepCellLines),
                       GlobalSupport=Scores$GlobalSupport,
                       SignatureLength=as.character(Scores$SignitureLength),
                       CELLectorScores=Scores$CELLectorScores,
                       Signature=as.character(Scores$SignaturesES),
                       stringsAsFactors = FALSE)
    nonRepCelLines<-as.character(setdiff(as.character(CELLlineData$CellLine),Scores$RepCellLines))

    remainin<-data.frame(CellLines=nonRepCelLines,
                         GlobalSupport=rep(0,length(nonRepCelLines)),
                         SignatureLength=rep(0,length(nonRepCelLines)),
                         CELLectorScores=rep(0,length(nonRepCelLines)),
                         Signature=rep('-',length(nonRepCelLines)))

    Scores<-rbind(Scores,remainin)
    return(Scores)
  }else{
    print('Error: alfa needs to be >=0 and <=1')
  }
}

CELLector.CMPs_getModelAnnotation <- function(URL='https://cog.sanger.ac.uk/cmp/download/model_list_latest.csv.gz'){
  if(url.exists(URL)){
    X <- read_csv(URL)
  }else{
    X <- NULL
  }
  return(X)
}

## documentation to be updated

CELLector.Build_Search_Space<-function(ctumours,
                                       cancerType,
                                       minlen=1,
                                       verbose=TRUE,
                                       mutOnly=FALSE,
                                       cnaOnly=FALSE,
                                       includeHMS=FALSE,
                                       minGlobSupp=0.01,
                                       FeatureToExclude=NULL,
                                       pathway_CFEs = NULL,
                                       pathwayFocused=NULL,
                                       subCohortDefinition=NULL,
                                       NegativeDefinition=FALSE,
                                       cnaIdMap=NULL,
                                       cnaIdDecode=NULL,
                                       hmsIdDecode=NULL,
                                       cdg=NULL,
                                       UD_genomics=FALSE){


  if(!UD_genomics){
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
                   cnaId_decode = cnaIdDecode,
                   hmsId_decode = hmsIdDecode)

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
                             cnaId_decode = cnaIdDecode,hmsId_decode = hmsIdDecode)

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
                            cnaId_decode = cnaIdDecode,hmsId_decode = hmsIdDecode)

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

  CC <- colors(distinct = TRUE)
  CC <- CC[setdiff(1:length(CC),c(grep('gray',CC),'black'))]
  CC <- rgb(t(col2rgb(CC)),maxColorValue = 255)

  COLORSbyLev <- CC[sample(length(CC))][1:nROOT$totalCount]

  names(COLORSbyLev)<-names(Get(Traverse(nROOT),'Names'))
  nROOT$Set(Colors=COLORSbyLev,traversal = 'level')


  tmpLabels<-Get(Traverse(nROOT,'level'),'name')


  treeLabels<-lapply(tmpLabels,function(x){
    y<-x
    start_<-str_locate_all(y,'[(]')[[1]][,1]
    end_<-str_locate_all(y,'[)]')[[1]][,1]

    for (mm in 1:length(start_)){
      n_start_<-str_locate_all(y,'[(]')[[1]][,1]
      n_end_<-str_locate_all(y,'[)]')[[1]][,1]
      y<-str_remove(y,paste('[(]',str_sub(y,n_start_[1]+1,n_end_[1]-1),'[)]',sep=''))
    }
    return(y)
  }
  )

  nROOT$Set(name=treeLabels,traversal='level')

  nodeIdx<-as.numeric(unlist(lapply(str_split(treeLabels,' '),function(x){x[1]})))

  COLORS<-rep(NA,length(nodeIdx))
  COLORS[nodeIdx]<-COLORSbyLev

  NT<-cbind(NT,COLORS)
  return(list(navTable=NT,TreeRoot=nROOT))
}

## Exported non Documented functions

CELLector.CELLline_buildBEM <- function(varCat=NULL,
                                        Tissue,
                                        Cancer_Type,
                                        Cancer_Type_details=NULL,
                                        sample_site=NULL,
                                        excludeOrganoids=FALSE,
                                        humanonly=TRUE,
                                        msi_status_select=NULL,
                                        gender_select=NULL,
                                        mutational_burden_th=NULL,
                                        age_at_sampling=NULL,
                                        ploidy_th=NULL,
                                        ethnicity_to_exclude=NULL,
                                        GenesToConsider=NULL,
                                        VariantsToConsider=NULL){

  if(length(varCat)==0){
    varCat<-CELLector.CMPs_getVariants()
    clAnnotation<-CELLector.CMPs_getModelAnnotation()
    clAnnotation$cancer_type_detail<-
      str_sub(clAnnotation$cancer_type_detail,3,end = str_length(clAnnotation$cancer_type_detail)-3)

    if(!excludeOrganoids){
      id<-which(clAnnotation$tissue==Tissue & is.element(clAnnotation$cancer_type,Cancer_Type))
    }else{
      id<-which(clAnnotation$tissue==Tissue & is.element(clAnnotation$cancer_type,Cancer_Type) & clAnnotation$model_type!='Organoid')
    }

    cls<-clAnnotation$model_id[id]
    varCat<-varCat[which(is.element(varCat$model_id,cls)),]
    clAnnotation<-clAnnotation[which(is.element(clAnnotation$model_id,cls)),]

    if(length(Cancer_Type_details)>0){
      id<-which(is.element(clAnnotation$cancer_type_detail,Cancer_Type_details))

      cls<-clAnnotation$model_id[id]
      varCat<-varCat[which(is.element(varCat$model_id,cls)),]
      clAnnotation<-clAnnotation[which(is.element(clAnnotation$model_id,cls)),]
    }




    if(length(sample_site)>0){
      id<-which(is.element(clAnnotation$sample_site,sample_site))

      cls<-clAnnotation$model_id[id]
      varCat<-varCat[which(is.element(varCat$model_id,cls)),]
      clAnnotation<-clAnnotation[which(is.element(clAnnotation$model_id,cls)),]

    }

    if(length(humanonly)>0){
      id<-which(clAnnotation$species=="Homo Sapiens")

      cls<-clAnnotation$model_id[id]
      varCat<-varCat[which(is.element(varCat$model_id,cls)),]
      clAnnotation<-clAnnotation[which(is.element(clAnnotation$model_id,cls)),]
    }

    if(length(msi_status_select)>0){
      id<-which(!is.na(clAnnotation$msi_status) & (clAnnotation$msi_status==msi_status_select |
                                                     (msi_status_select=='MSI-L/H' & (clAnnotation$msi_status=='MSI-L' | clAnnotation$msi_status=='MSI-H'))))
      cls<-clAnnotation$model_id[id]
      varCat<-varCat[which(is.element(varCat$model_id,cls)),]
      clAnnotation<-clAnnotation[which(is.element(clAnnotation$model_id,cls)),]
    }

    if(length(gender_select)>0){
      id<-which(is.element(clAnnotation$gender,gender_select))
      cls<-clAnnotation$model_id[id]
      varCat<-varCat[which(is.element(varCat$model_id,cls)),]
      clAnnotation<-clAnnotation[which(is.element(clAnnotation$model_id,cls)),]
    }

    if(length(ethnicity_to_exclude)>0){
      id<-which(!is.element(clAnnotation$ethnicity,ethnicity_to_exclude))
      cls<-clAnnotation$model_id[id]
      varCat<-varCat[which(is.element(varCat$model_id,cls)),]
      clAnnotation<-clAnnotation[which(is.element(clAnnotation$model_id,cls)),]
    }

    if(length(mutational_burden_th)>0){
      id<-which(round(clAnnotation$mutational_burden)>=mutational_burden_th[1] & round(clAnnotation$mutational_burden)<=mutational_burden_th[2])
      cls<-clAnnotation$model_id[id]
      varCat<-varCat[which(is.element(varCat$model_id,cls)),]
      clAnnotation<-clAnnotation[which(is.element(clAnnotation$model_id,cls)),]
    }



    if(length(ploidy_th)>0){
      id<-which(round(clAnnotation$ploidy)>=ploidy_th[1] & round(clAnnotation$ploidy)<=ploidy_th[2])
      cls<-clAnnotation$model_id[id]
      varCat<-varCat[which(is.element(varCat$model_id,cls)),]
      clAnnotation<-clAnnotation[which(is.element(clAnnotation$model_id,cls)),]
    }

    if(length(age_at_sampling)>0){
      id<-which(round(clAnnotation$age_at_sampling)>=age_at_sampling[1] & round(clAnnotation$age_at_sampling)<=age_at_sampling[2])
      cls<-clAnnotation$model_id[id]
      varCat<-varCat[which(is.element(varCat$model_id,cls)),]
      clAnnotation<-clAnnotation[which(is.element(clAnnotation$model_id,cls)),]
    }

    if(length(GenesToConsider)>0){
      varCat<-varCat[which(is.element(varCat$gene_symbol,GenesToConsider)),]
    }

    if(length(VariantsToConsider)>0){
      sigs<-paste(varCat$gene_symbol,varCat$cdna_mutation,paste('p.',varCat$aa_mutation,sep=''))
      varCat<-varCat[which(is.element(sigs,VariantsToConsider)),]
    }
  }

 # print(sort(cls))

  allModels<-sort(unique(varCat$model_id))
  allModel_ids<-varCat$model_id[match(allModels,varCat$model_id)]

  allGenes<-sort(unique(varCat$gene_symbol))

  BEM<-do.call(what = cbind,lapply(allModels,function(x){
    is.element(allGenes,varCat$gene_symbol[varCat$model_id==x])+0
  }))
  rownames(BEM)<-allGenes

  cls<-clAnnotation$model_name[match(allModel_ids,clAnnotation$model_id)]

  BEM<-data.frame(CMP_identifier=allModel_ids,
                  CellLine=cls,
                  t(BEM))



  return(BEM)
}

CELLector.Tumours_buildBEM <- function(varCat=NULL,
                                       Cancer_Type,
                                       GenesToConsider=NULL,
                                       VariantsToConsider=NULL){


  if(length(varCat)==0){
    data(CELLector.PrimTumVarCatalog)
    varCat<-CELLector.PrimTumVarCatalog

    sampleN<-varCat$SAMPLE[varCat$Cancer.Type==Cancer_Type]

    varCat<-varCat[which(is.element(varCat$SAMPLE,sampleN)),]
  }

  if(length(GenesToConsider)>0){
    varCat<-varCat[which(is.element(varCat$Gene,GenesToConsider)),]
  }

  if(length(VariantsToConsider)>0){
    sigs<-paste(varCat$Gene,varCat$cDNA,varCat$AA)
    varCat<-varCat[which(is.element(sigs,VariantsToConsider)),]
  }

  allSamples<-sort(unique(varCat$SAMPLE))

  allGenes<-sort(unique(varCat$Gene))

  BEM<-do.call(what = cbind,lapply(allSamples,function(x){
    is.element(allGenes,varCat$Gene[varCat$SAMPLE==x])+0
  }))

  rownames(BEM)<-allGenes
  colnames(BEM)<-allSamples

  return(BEM)
}

CELLector.CMPs_getDriverGenes <- function(URL='https://cog.sanger.ac.uk/cmp/download/cancer_genes_latest.csv.gz'){
  if(url.exists(URL)){
    X <- read_csv(URL)
    X <- X$gene_symbol
  }else{
    X <- NULL
  }
  return(X)
}
CELLector.CMPs_getVariants <- function(URL='https://cog.sanger.ac.uk/cmp/download/mutations_2018-08-01_1640.csv.gz'){
  if(url.exists(URL)){
    X<-read_csv(URL)

  }else{
    X <- NULL
  }
  return(X)
}



### not documented data objects:

## CELLector.CellLine.BEMs_v2
## CELLector.PrimTum.BEMs_v2
## CELLector.CFEsV2
## CELLector.CFEs.HMSid_decode


## not Exported functions
sunBurstFormat<-function(searchSpace){
  table_tree<-data.frame(lapply(searchSpace$navTable[,1:11], as.character), stringsAsFactors=FALSE)

  table_tree$Left.Child.Index[table_tree$Left.Child.Index==0]<- -1
  table_tree$Right.Child.Index[table_tree$Right.Child.Index==0]<- -1

  table_tree<-rbind(c(0,'TOTAL','TOTAL','root','-1',table_tree$CurrentTotal[1],table_tree$CurrentTotal[1],1,1,0),table_tree)
  table_tree$Type[2]<-'Left.Child'

  table_tree$Left.Child.Index[1]<-1
  table_tree$Right.Child.Index[1]<- -1

  table_tree$Idx<-as.numeric(table_tree$Idx)+1
  table_tree$Parent.Idx<-as.numeric(table_tree$Parent.Idx)+1
  table_tree$Left.Child.Index<-as.numeric(table_tree$Left.Child.Index)+1
  table_tree$Right.Child.Index<-as.numeric(table_tree$Right.Child.Index)+1

  leaves<-which(table_tree$Right.Child.Index==0)
  leaves<-leaves[-1]

  stable_tree<-table_tree

  for (i in 1:length(leaves)){

    CurrentTotal<-as.numeric(table_tree$CurrentTotal[leaves[i]])-
      as.numeric(table_tree$AbsSupport[leaves[i]])



    stable_tree<-rbind(stable_tree,
                       c(nrow(stable_tree)+1,
                         'Others',
                         'Others',
                         'Right.Child',
                         leaves[i],
                         CurrentTotal,
                         CurrentTotal,
                         1,
                         CurrentTotal/as.numeric(stable_tree$CurrentTotal[1]),0,0))
    stable_tree$Right.Child.Index[leaves[i]]<-nrow(stable_tree)
  }


  nnodes<-nrow(stable_tree)
  edgeList<-NULL

  for (i in 1:nnodes){

    currentNode<-i

    if (stable_tree$Type[currentNode]=='Left.Child'){
      edgeList<-c(edgeList,as.numeric(c(stable_tree$Parent.Idx[currentNode],stable_tree$Idx[currentNode])))
      #  print(c(stable_tree$Parent.Idx[currentNode],stable_tree$Idx[currentNode]))
    }else{
      startingNode<-currentNode
      while (stable_tree$Type[startingNode]=='Right.Child'){
        startingNode<-as.numeric(stable_tree$Parent.Idx[startingNode])

      }
      if(currentNode!=1){
        edgeList<-c(edgeList,as.numeric(c(stable_tree$Parent.Idx[startingNode],currentNode)))
      }
    }
  }

  G<-make_graph(edgeList)
  leaves<-V(G)[which(degree(G,v=V(G),'out')==0)]
  paths<-all_simple_paths(G,1,leaves)
  nleaves<-length(paths)

  chainP<-vector()
  npat<-vector()
  for (i in 1:nleaves){

    currentId<-as.numeric(paths[[i]][length(paths[[i]])])

    currentItemDec<-stable_tree$ItemsDecoded[as.numeric(as.character(paths[[i]]))]

    currentItemDec<-unlist(lapply(str_split(currentItemDec,'[(]'),
                                  function(x){x[1]}))

    chainP[i]<-paste(paste(as.numeric(paths[[i]])-1,currentItemDec),
                     collapse='-')

    npat[i]<-as.numeric(stable_tree$AbsSupport[currentId])
  }

  sequences<-data.frame(V1=chainP,V2=npat,stringsAsFactors = FALSE)



  return(sequences)

}
createHtmlNodeProperties<-function(LocalSearchSpace,CLdataset=NULL){

  tree<-LocalSearchSpace$TreeRoot

  nodeIds<-as.numeric(unlist(lapply(str_split(names(Get(Traverse(LocalSearchSpace$TreeRoot,traversal = 'level'),'NodeType')),' '),
                                    function(x){x[1]})))

  nnodes<-length(nodeIds)

  signatures<-CELLector.createAllSignatures(NavTab = LocalSearchSpace$navTable)

  SS<-signatures

  signatures<-signatures$S[nodeIds]

  if(length(CLdataset)>0){
    modelMat<-CELLector.buildModelMatrix(Sigs = SS$ES,dataset = CLdataset,searchSpace = LocalSearchSpace$navTable)


  }


  nodeTypes<-as.character(LocalSearchSpace$navTable$Type[nodeIds])
  parents<-LocalSearchSpace$navTable$Parent.Idx[nodeIds]

  NT<-nodeTypes
  nodeTypes[nodeTypes=='Right.Child']<-paste('Complement of SubType ',parents[nodeTypes=='Right.Child'])
  nodeTypes[nodeTypes=='Left.Child']<-paste('Refinement of SubType ',parents[nodeTypes=='Left.Child'])

  typeColors<-rep('Gray',tree$GlobalSupp)
  typeColors[which(NT=='Right.Child')]<-'Tomato'
  typeColors[which(NT=='Left.Child')]<-'MediumSeaGreen'

  npatients<-LocalSearchSpace$navTable$AbsSupport[nodeIds]
  percOnTotal<-LocalSearchSpace$navTable$GlobalSupport[nodeIds]
  percOnPop<-LocalSearchSpace$navTable$PercSupport[nodeIds]

  html_node_summaries<-vector()
  for (i in 1:nnodes){
    header<-'<!DOCTYPE html><html><head><title>'
    TITLE<-paste('Patient SubType id:',nodeIds[i])
    postTitle<-'</title></head><body>'
    pageContent<-paste('<p style="font-size:15px;"><b>Patient SubType:',nodeIds[i],'</b></p>')
    pageContent<-paste(pageContent,'<b>Underlying signature:</b><br /><i>',signatures[i],'</i><br /><br />')
    pageContent<-paste(pageContent,'<p style="background-color:',typeColors[i],';">',nodeTypes[i],'</p><br />')
    pageContent<-paste(pageContent,'<b>N of Patients:</b>',npatients[i],'<br />')
    pageContent<-paste(pageContent,format(100*percOnTotal[i],digits=3),'% of total <br />',sep='')
    if (NT[i]!='root'){
      if(NT[i]=='Right.Child'){
        pageContent<-paste(pageContent,format(100*percOnPop[i],digits=3),'% of subType',parents[i],' complement<br />',sep='')
      }else{
        pageContent<-paste(pageContent,format(100*percOnPop[i],digits=3),'% of subType',parents[i],'<br />',sep='')
      }

    }

    if(length(CLdataset)>0){
      if (!is.element(nodeIds[i],rownames(modelMat))){
        pageContent<-paste(pageContent,'<p style="color:Tomato;"><b>No cell lines</b></p><br />')
      }else{
        pageContent<-paste(pageContent,'<p style="color:MediumSeaGreen;"><b>',
                           sum(modelMat[as.character(nodeIds[i]),]),'cell lines</b></p><br />')
      }
    }
    tail<-'</body></html>'

    html_node_summaries[i]<-paste(header,TITLE,postTitle,pageContent,tail,sep='')
  }

  return(html_node_summaries)

}
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

      currentTerm<-str_trim(as.character(NavTab$ItemsDecoded[pos]))
      currentTerm<-paste(paste(prefix,unlist(str_split(currentTerm,', ')),sep=''),collapse=', ')

      #currentTerm<-paste(prefix,str_trim(as.character(NavTab$ItemsDecoded[pos])),sep='')

      SIGNATURE<-paste(currentTerm,SIGNATURE,sep=', ')

      EcurrentTerm<-str_trim(as.character(NavTab$Items[pos]))
      EurrentTerm<-paste(paste(prefix,unlist(str_split(currentTerm,', ')),sep=''),collapse=', ')

      EcurrentTerm<-paste(prefix,NavTab$Items[pos],sep='')
      encodedSIGNATURE<-paste(EcurrentTerm,encodedSIGNATURE,sep=', ')
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
                     hmsId_decode,
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

      dIS<-decodeSIG(ctype = ctype,codedSIG = IS, cnaId_decode = cnaId_decode, hmsId_decode = hmsId_decode,cdg = cdg)


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

decodeSIG<-function(ctype,codedSIG,cnaId_decode,hmsId_decode=NULL,cdg){
  IS<-codedSIG
  icna<-grep('cna',IS)

  #if(length(icna)>0){



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

  ihms<-grep('hms',IS)

  if (length(ihms)>0){

    cnalu<-CELLector.hms_look_up(IS[ihms],TCGALabel = ctype,hmsId_decode = hmsId_decode)

    genes<-cnalu$GN

    IS[ihms]<-paste('HypMet_',genes,sep='')


  }

  noncna<-setdiff(1:length(IS),c(icna,ihms))

  IS[noncna]<-paste(IS[noncna],'mut',sep='')

  return(IS)
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


