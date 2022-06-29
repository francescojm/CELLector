## Exported Documented functions
CELLector.from_Hierarchical_to_Partition <- function(navTable,
                                                     samples_id,
                                                     verbose = TRUE){

  # go through binary tree to convert hierachical configuration to partition (auxiliary)

  N_tot <- length(samples_id)
  navTable_partition <- data.frame(Idx = navTable$Idx, Signature = NA, SignatureDecoded = NA,
                                   Points = NA, Total = NA, Support = NA, COLORS = NA)

  for(id_row in 1:nrow(navTable)){

    if (verbose) {
      print(paste("adapting node index", navTable$Idx[id_row],
                  ": feature", navTable$ItemsDecoded[id_row]))
    }

    positivePoints <- navTable$positivePoints[id_row]
    samples_tmp <- strsplit(positivePoints, split = ',')[[1]]
    negative_feat <- c()
    negative_feat_dec <- c()
    samples_to_remove <- c()

    # check if there is a left child
    if(navTable$Left.Child.Index[id_row]!=0){

      new_id <- navTable$Left.Child.Index[id_row]

      left_child_ids <- navTable$positivePoints[navTable$Idx == new_id]
      samples_to_remove <- strsplit(left_child_ids, split = ',')[[1]]
      negative_feat <- c(negative_feat, navTable$Items[navTable$Idx == new_id])
      negative_feat_dec <- c(negative_feat_dec, navTable$ItemsDecoded[navTable$Idx == new_id])

      # add all the right nodes
      exists_right <- navTable$Right.Child.Index[navTable$Idx == new_id] != 0

      while(exists_right){

        new_id <- navTable$Right.Child.Index[navTable$Idx == new_id]

        right_child_ids <- navTable$positivePoints[navTable$Idx == new_id]
        samples_to_remove <- c(samples_to_remove, strsplit(right_child_ids, split = ',')[[1]])
        exists_right <- navTable$Right.Child.Index[navTable$Idx == new_id] != 0
        negative_feat <- c(negative_feat, navTable$Items[navTable$Idx == new_id])
        negative_feat_dec <- c(negative_feat_dec, navTable$ItemsDecoded[navTable$Idx == new_id])

      }
    }

    samples_partition <- setdiff(samples_tmp, samples_to_remove)
    navTable_partition$Points[id_row] <- paste0(samples_partition, collapse = ',')
    navTable_partition$Total[id_row] <- length(samples_partition)
    navTable_partition$Support[id_row] <- length(samples_partition)/N_tot
    navTable_partition$COLORS[id_row] <- navTable$COLORS[id_row]

    CL_rule <-  createRuleFromNode(navTable, nodeIdx = id_row)
    navTable_partition$Signature[id_row] <- paste(c(CL_rule$ES,
                                                    sprintf(' ~%s', negative_feat)),
                                                  collapse = ',')

    navTable_partition$SignatureDecoded[id_row] <- paste(c(CL_rule$S,
                                                           sprintf(' ~%s', negative_feat_dec)),
                                                         collapse = ',')


  }

  tot_samples_mapped <- unlist(str_split(navTable_partition$Points, ','))

  # add class with no frequent features
  if(any(!samples_id %in% tot_samples_mapped)){

    if (verbose) {
      print(sprintf("adapting node index 0: no frequent feature"))
    }

    samples_partition <- setdiff(samples_id, tot_samples_mapped)
    sign_negative <- sprintf('~%s', navTable$Items[1])
    sign_negative_dec <- sprintf('~%s', navTable$ItemsDecoded[1])
    exists_right <- navTable$Right.Child.Index[1] != 0
    new_id <- 1
    while(exists_right){
      new_id <- navTable$Right.Child.Index[navTable$Idx == new_id]
      sign_negative <- c(sign_negative, sprintf('~%s', navTable$Items[navTable$Idx == new_id]))
      sign_negative_dec <- c(sign_negative_dec, sprintf('~%s', navTable$ItemsDecoded[navTable$Idx == new_id]))
      exists_right <- navTable$Right.Child.Index[navTable$Idx == new_id] != 0
    }

    group_negative <- data.frame(Idx = 0,
                                 Signature = paste(sign_negative, collapse = ', '),
                                 SignatureDecoded = paste(sign_negative_dec, collapse = ', '),
                                 Points = paste0(samples_partition, collapse = ','),
                                 Total = length(samples_partition),
                                 Support = length(samples_partition)/N_tot,
                                 COLORS = "#000000")

    navTable_partition <- rbind(group_negative, navTable_partition)

  }

  return(navTable_partition)
}
