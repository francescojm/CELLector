CELLector_App.sunBurstFormat<-function(table_tree){

  table_tree<-data.frame(lapply(table_tree[,1:11], as.character), stringsAsFactors=FALSE)
  table_tree$Left.Child.Index[table_tree$Left.Child.Index==0]<- -1
  table_tree$Right.Child.Index[table_tree$Right.Child.Index==0]<- -1

  table_tree<-rbind(c(0,'TOTAL','TOTAL','root','-1',table_tree$CurrentTotal[1],table_tree$CurrentTotal[1],1,1,0),table_tree)
  table_tree$Type[2]<-'Left.Child'

  table_tree$Idx<-as.numeric(table_tree$Idx)+1
  table_tree$Parent.Idx<-as.numeric(table_tree$Parent.Idx)+1
  table_tree$Left.Child.Index<-as.numeric(table_tree$Left.Child.Index)+1
  table_tree$Right.Child.Index<-as.numeric(table_tree$Right.Child.Index)+1

  leaves<-which(table_tree$Left.Child.Index==0 & table_tree$Right.Child.Index==0)

  stable_tree<-table_tree

  for (i in 1:length(leaves)){

    CurrentTotal<-as.numeric(table_tree$CurrentTotal[leaves[i]])-as.numeric(table_tree$AbsSupport[leaves[i]])



    stable_tree<-rbind(stable_tree,c(nrow(stable_tree)+1,'Others','Others','Right.Child',leaves[i],CurrentTotal,CurrentTotal,
                                     1,CurrentTotal/as.numeric(stable_tree$CurrentTotal[1]),0,0))
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
        #print(c(stable_tree$Parent.Idx[startingNode],currentNode))
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


    chainP[i]<-paste(paste(as.numeric(paths[[i]])-1,stable_tree$ItemsDecoded[as.numeric(as.character(paths[[i]]))]),
                     collapse='-')

    npat[i]<-as.numeric(stable_tree$AbsSupport[currentId])
  }

  sequences<-data.frame(V1=chainP,V2=npat,stringsAsFactors = FALSE)


  # labes<-unique(unlist(strsplit(sequences$V1,"-")))
  #
  # n <- length(labes)
  # qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  # col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  #
  # # define specific colors
  # colors <- list(
  #   domain=unique(unlist(strsplit(sequences$V1,"-"))),
  #   range=sample(col_vector, n)
  # )
  #
  #
  # colors$range[colors$domain=='Others']<-'white'
  #

  return(sequences)

}

#library(CELLector)

library(sunburstR)
library(collapsibleTree)

data(CELLector.PrimTum.BEMs)
data(CELLector.Pathway_CFEs)
data(CELLector.CFEs.CNAid_mapping)
data(CELLector.CFEs.CNAid_decode)
data(CELLector.HCCancerDrivers)


tumours_BEM<-CELLector.PrimTum.BEMs$COREAD

#CELLector.mostSupported_CFEs(t(tumours_BEM),minlen = 1)


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

sb<-CELLector_App.sunBurstFormat(tmp$navTable)

domain<-unique(unlist(strsplit(sb$V1,'-')))
tmpCol<-rep('black',length(domain))

colors <- list(
  domain=domain,
  range=tmpCol
)

tsb<-sunburst(sb,breadcrumb = list(w = 100),colors = colors,percent = FALSE,count = FALSE,
             explanation = "function(d) { return d.data.name}")
tsb


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






sb$x$tasks <- htmlwidgets::JS(
  'function(){d3.select(this.el).select("#" + this.el.id + "-trail").style("font-size","60%")}'
)

sb



