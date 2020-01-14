# #Load data
# myfile<-system.file("extdata", "luad_tcga.rda", package = "spathial")
# data<-load(myfile)
# #X <- data[,2:ncol(data)]
# #rownames(X)<-data[,1]
# #X_labels<-ifelse(substr(rownames(X),14,15)%in%c("01"),"tumor","normal")
# #names(X_labels)<-rownames(X)
#
# #Compute centroids
# normal_centroid <- colMeans(X[which(X_labels == "normal"),], na.rm = TRUE)
# tumor_centroid <- colMeans(X[which(X_labels == "tumor"),], na.rm = TRUE)
#
# #Subdivide normal samples from tumor samples
# normal <- X[which(X_labels == "normal"),]
# tumor <- X[which(X_labels == "tumor"),]
#
# startPoint <- getMaxDistancePoint(tumor_centroid, normal)
# endPoint <- getMaxDistancePoint(normal_centroid, tumor)
#
# getMaxDistancePoint <- function(centroid, samples){
#   library(pracma)
#
#   centroid <- t(centroid)
#   dst<-pracma::distmat(
#     as.matrix(centroid),
#     as.matrix(samples)
#   )
#   ord <- order(-dst)
#   return(rownames(samples)[ord[1]])
# }
#
# #Set centroids with the rownames previously computed
# boundaries <- spathial::spathialBoundaryIds(X, X_labels, mode=3, from=startPoint, to=endPoint)
# boundary_ids<-boundaries$boundary_ids
# X <- boundaries$X
# X_labels <- boundaries$X_labels
#
# spathial_res <- spathial::spathialWay(X, boundary_ids, NC = 50)
#
# spathial::spathialPlot(X, X_labels, boundary_ids, spathial_res, perplexity_value = 30, mask = NULL)
#
#
# #COMPARISON
# #Spathial correlation
# res_stat <- spathial::spathialStatistics(spathial_res)
# pathcor_model_dist<-sort(setNames(res_stat$correlations,colnames(spathial_res)),dec=TRUE)
#
# ### Cancer Gene Census
# load(system.file("extdata", "cgc_genelists.rda", package = "spathial"))
# if(TRUE){ # Remove genes that are both TSG and oncogenes
#   tmp<-setdiff(oncogenes,tsg)
#   tsg<-setdiff(tsg,oncogenes)
#   oncogenes<-tmp
#   rm(tmp)
# }
#
# library(edgeR)
#
# #labels<-ifelse(substr(rownames(X),14,15)%in%c("01"),"tumor","normal")
# #labels[substr(rownames(X),14,15)=="06"]<-"metastasis"
# #names(labels)<-rownames(X)
#
# rawcounts <- t(X)
#
# group<-as.factor(X_labels)
# group<-relevel(group,ref="normal")
# edger<-DGEList(rawcounts,group=group,genes=rownames(rawcounts))
# design<-model.matrix(~0+group)
# colnames(design)<-levels(group)
# keep <- filterByExpr(edger, design) # Filtering to remove low counts
# table(keep)
# edger<-edger[keep,,keep.lib.sizes=FALSE]
# edger<-calcNormFactors(edger) # Normalization for composition bias
# edger<-estimateDisp(edger,design,robust=TRUE) # Dispersion estimation
# fit<-glmQLFit(edger,design,robust=TRUE) # Quasi-Likelihood dispersion estimate
# contrast<-makeContrasts(tumor-normal,levels=design)
# res<-glmQLFTest(fit,contrast=contrast)
# res<-res$table
# res$padj<-p.adjust(res$PValue,method="BH")
# colnames(res)<-c("logFC","logCPM","F","pvalue","padj")
# edger<-res
#
# compareTool(oncogenes, tsg, edger, pathcor_model_dist, "edger_vs_spathial_dist_new")
#
# compareTool <- function(oncogenes, tsg, edger, pathcor, filename){
#   ### Make them comparable
#   common<-intersect(names(pathcor),rownames(edger))
#   x<-setNames(-log10(edger$padj)*sign(edger$logFC),rownames(edger))
#   x<-x[common]
#   y<-pathcor[common]
#   xrank<-rank(x,ties.method="random")
#   yrank<-rank(y,ties.method="random")
#   revxrank<-rank(-x,ties.method="random")
#   revyrank<-rank(-y,ties.method="random")
#
#   ### Most different genes Path vs. DE
#   #ponc<-intersect(names(x[x<quantile(y,p=0.999)]),ponc)
#   #ptsg<-intersect(names(x[x>quantile(y,p=0.001)]),ptsg)
#
#   ponc<-intersect(names(y[y>quantile(y,p=0.7)]),oncogenes)
#   ptsg<-intersect(names(y[y<quantile(y,p=0.3)]),tsg)
#   ratio<-xrank[ponc]^2-yrank[ponc]^2
#   ponc<-names(sort(ratio))[1:30]
#   ratio<-revxrank[ptsg]^2-revyrank[ptsg]^2
#   ptsg<-names(sort(ratio))[1:30]
#
#   plotComparison(x, y, ptsg, ponc, xrank, yrank, filename)
#   #plotComparisonBN(x, y, ptsg, ponc, xrank, yrank, paste0(filename, "_BN"))
# }
#
# plotComparison <- function(x, y, ptsg, ponc, xrank, yrank, filename){
#   if(!is.na(ptsg) && !is.na(ponc)){
#     set.seed(1)
#     png(paste0(filename, ".png"),w=3000,h=3000,res=500,pointsize = 12)
#     plot(xrank,yrank,col="lightgrey",
#          xlim=c(1,20000),
#          # ylim=c(-3000,22000),
#          xlab="edgeR Rank",pch=20,
#          ylab="Path Rank"
#     )
#     points(xrank[ponc],yrank[ponc],col="#FF000033",pch=20)
#     points(xrank[ptsg],yrank[ptsg],col="#0000FF33",pch=20)
#     textplot3(xrank[ponc],yrank[ponc],words=ponc,col="red3",font=2,cex=1.1,line.col="#FF000033",line.width=2)
#     textplot3(xrank[ptsg],yrank[ptsg],words=ptsg,col="navy",font=2,cex=1.1,line.col="#0000FF33",line.width=2)
#     abline(h=yrank[which.min(abs(y))])
#     abline(v=xrank[which.min(abs(x))])
#     legend("bottomright",legend=c("Oncogene","TSG"),text.col=c("red3","navy"),text.font=2,cex=1.3)
#     dev.off()
#   }else{
#     message("No common values found")
#   }
# }
#
# plotComparisonBN <- function(x, y, ptsg, ponc, xrank, yrank, filename){
#   if(!is.na(ptsg) && !is.na(ponc)){
#     set.seed(1)
#     png(paste0(filename, ".png"),w=3000,h=3000,res=400)
#     plot(xrank,yrank,col="#D8D8D8",
#          xlim=c(1,20000),
#          # ylim=c(-3000,22000),
#          xlab="edgeR Rank",pch=20,
#          ylab="Principal Path Rank"
#     )
#     points(xrank[ponc],yrank[ponc],col="gray",pch=20)
#     points(xrank[ptsg],yrank[ptsg],col="gray",pch=20)
#     textplot3(xrank[ponc],yrank[ponc],words=ponc,col="black",font=2,cex=1.1,line.col="#55555B",line.width=2)
#     textplot3(xrank[ptsg],yrank[ptsg],words=ptsg,col="#4D4D4F",font=2,cex=1.1,line.col="#55555B",line.width=2)
#     abline(h=yrank[which.min(abs(y))])
#     abline(v=xrank[which.min(abs(x))])
#     legend("bottomright",legend=c("Oncogene","TSG"),text.col=c("black","#4D4D4F"),text.font=2,cex=1.3)
#     dev.off()
#   }else{
#     message("No common values found")
#   }
# }
#
# library(wordcloud)
# textplot3<-function(x, y, words, cex = 1, pch = 16, pointcolor = "#FFFFFF00",line.width=1,
#                     new = FALSE, show.lines = TRUE, line.col="black",rstep=0.5,...)
# {
#   if (new) {
#     plot(x, y, type = "n", ...)
#   }
#   lay <- wordlayout(x, y, words, cex, rstep=rstep, ...)
#   if (show.lines) {
#     for (i in seq_len(length(x))) {
#       xl <- lay[i, 1]
#       yl <- lay[i, 2]
#       w <- lay[i, 3]
#       h <- lay[i, 4]
#       if (x[i] < xl || x[i] > xl + w || y[i] < yl || y[i] >
#           yl + h) {
#         points(x[i], y[i], pch = pch, col = pointcolor,
#                cex = 0.5)
#         nx <- xl + 0.5 * w
#         ny <- yl + 0.5 * h
#         lines(c(x[i], nx), c(y[i], ny), col = line.col,lwd=line.width)
#       }
#     }
#   }
#   text(lay[, 1] + 0.5 * lay[, 3], lay[, 2] + 0.5 * lay[, 4],
#        words, cex = cex, ...)
# }
