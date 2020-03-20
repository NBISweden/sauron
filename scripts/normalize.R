install.packages("glmer")
library(future)
library(MASS)

gene_chunk <- cut(1:nrow(a),breaks = 50)
gene_chunk <- lapply( levels(gene_chunk) , function(x){ rownames(a)[gene_chunk==x] })

plan(multiprocess, workers=(availableCores()-1) )
f <- list()
for(ind in 1:length(gene_chunk) ){
  f[[ ind ]] <- future({
      normlzdata <- apply(a@assays$RNA@counts[gene_chunk[[ind]],],1,function(x){
      try(mod <- glm(x ~ a$Plate+a$perc_rpl+a$perc_rps+a$perc_mito+a$nFeature_RNA+a$nCount_RNA+a$perc_b2+a$perc_tms,
                 family = negative.binomial(theta=1)))
      return(log2(round((mod$residuals+1)*(1-mod$coefficients[1]),1)+1))
      })
      return( Matrix::t( Matrix::Matrix(normlzdata,sparse = T) ) )
  })
}

norm_data <- lapply(f, FUN = value)
norm_data2 <- do.call(rbind,norm_data)
dim(norm_data2)
#######

