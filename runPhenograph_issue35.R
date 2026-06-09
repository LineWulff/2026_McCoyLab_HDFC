runPhenograph_issue35 <- function (fcd,
                                   input_type,
                                   data_slot,
                                   k,
                                   seed = 91,
                                   prefix = NULL,
                                   nPC = ncol(fcd$pca[[data_slot]]),
                                   markers = colnames(fcd$expr[[data_slot]]),
                                   discard = FALSE){
  set.seed(seed)
  
  # see if selected markers are present in condor_object , if input_type = "expr"
  if (input_type == "expr"){
    for (single in markers){
      if (!single %in% colnames(fcd[["expr"]][["orig"]])){
        stop(paste("ERROR:",single, "not found in expr markers."))
      }
    }
    
    # define markers to use
    if (discard == FALSE){              # (discard == F -> keep specified markers (default = all))
      
      phe_markers <- markers
    }
    
    else if (discard == TRUE) {       # (discard == T -> discard specified markers, error if no markers are specified)
      
      if (length(markers) == length(colnames(fcd$expr[[data_slot]]))){
        
        stop("ERROR: No markers specified. Specify markers to be removed or set 'discard = F'.")
      }
      
      else {
        
        phe_markers <- setdiff(colnames(fcd$expr[[data_slot]]), markers)
        
      }
    }
    
    #define fcd subset for calculation
    data1 <- fcd$expr[[data_slot]][, colnames(fcd$expr[[data_slot]]) %in% phe_markers, drop = F]
    
    
    
  }
  if (input_type == "pca"){
    
    #define fcd subset for calculations and get used markers of PCA analysis
    
    data1 <- fcd$pca[[data_slot]][,1:nPC]
    phe_markers <- used_markers(fcd,  input_type = "pca", data_slot = data_slot, mute = T)
  }
  
  
  
  
  Rphenograph_out <- Rphenograph::Rphenograph(data1, k = k)
  Rphenograph_out <- as.matrix(igraph::membership(Rphenograph_out[[2]]))
  Rphenograph_out <- as.data.frame(matrix(ncol = 1,
                                          data = Rphenograph_out,
                                          dimnames = list(rownames(fcd$expr$orig), "Phenograph")))
  Rphenograph_out$Phenograph <- as.factor(Rphenograph_out$Phenograph)
  Rphenograph_out$Description <- paste(input_type, "_", data_slot, "_k", k, sep = "")
  
  phe_name <- paste("phenograph", sub("^_", "", paste(prefix, input_type, data_slot, "k", k, sep = "_")), sep = "_")
  
  if (input_type == "pca"){
    if (nPC < ncol(fcd[[input_type]][[data_slot]])) {
      
      suffix <- paste0("top", nPC)
      
      phe_name <- paste(phe_name, suffix, sep = "_")
    }
  }
  
  fcd[["clustering"]][[phe_name]] <- Rphenograph_out
  
  #save used markers in "extras"-slot
  fcd[["extras"]][["markers"]][[paste(phe_name,"markers", sep = "_")]] <- phe_markers
  
  
  return(fcd)
}