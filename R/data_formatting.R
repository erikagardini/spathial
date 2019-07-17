# myfile<-system.file("extdata", "liver_tcga.csv", package = "spathial")
# data<-read.csv(myfile,as.is=TRUE,header=TRUE)
# X <- data[,2:(dim(data)[2]-1)]
# X_labels <- data[,dim(data)[2]]
# rownames(X)<-data[,1]
# rownames(X)<-gsub("\\.","-",rownames(X))
#
# #Variance filtering
# varFiltering <- TRUE
# nvars <- 100
# if(varFiltering){
#   vars <- apply(X,2,var)
#   topfeatures <- names(sort(vars, dec=TRUE))[1:nvars]
#   X <- X[,topfeatures]
# }
#
# X$Y <- X_labels
# write.csv(X, "liver_tcga_example1.csv")
