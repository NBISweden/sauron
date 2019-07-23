
fast_ScaleData <- function(DATA, assay="RNA",vars.to.regress=NULL,do.scale=T,do.center=T,scale.max=NULL){
  
  #Check spelling to see if the variables to regress are in the dataset
  vars.to.regress <- vars.to.regress[vars.to.regress %in% colnames(DATA@meta.data)]
  cat("The following variables were found in your data and will be used to regress out the counts\n")
  print(vars.to.regress)
  
  #plan(strategy = "multicore", workers = nbrOfWorkers())
  
  #Regress factors, if any
  if(!is.null(vars.to.regress)){
    cat("Your data contains ",ncol(eval(parse(text=paste0("DATA@assays$",assay,"@data")))),"samples and",nrow(eval(parse(text=paste0("DATA@assays$",assay,"@data"))))," features\n")
    cat("fast_ScaleData regression started running, please wait ...\n")
    ttt <- Sys.time()
    l1 <- apply(eval(parse(text=paste0("DATA@assays$",assay,"@data"))),1,function (x){
      model <- paste0("x ~ ", paste(paste0("DATA$",vars.to.regress),collapse = "+"))
      m <- glm(eval(parse(text=model)))
      return(scale(m$residuals,T,T))
    })
    cat("fast_ScaleData ran in ",(Sys.time() - ttt)," minutes\n")
  }
  invisible(gc())
  
  #Scale data
  l2 <- apply(l1,ifelse(!is.null(vars.to.regress),2,1),function (x){  return(scale(x,do.scale,do.center)) })
  rownames(l2) <- colnames(DATA@assays[[assay]]@data)
  if(!is.null(scale.max)){ l2[l2 >= scale.max] <- scale.max }
  rm(l1)
  invisible(gc())
  
  #Assign it back to the DATA object
  DATA@assays[[assay]]@scale.data <-  as.matrix(t(l2))
  return(DATA)
}

