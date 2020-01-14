# myfile<-system.file("extdata", "karlsson-fpkms.rda", package = "spathial")
# load(myfile)
#
# genevars<-names(sort(apply(fpkms,1,var),dec=TRUE))
# expmat3<-fpkms[genevars,]
# X <- t(expmat3)
# X_labels <- annotation
#boundaries <- spathial::spathialBoundaryIds(X, X_labels, mode=2, from="G1", to="G2/M")
# boundary_ids<-boundaries$boundary_ids
# X <- boundaries$X
# X_labels <- boundaries$X_labels
#
# spathial_res <- spathial::spathialWay(X, boundary_ids, NC = 50)
# spathial::spathialPlot(X, X_labels, boundary_ids, spathial_res, perplexity_value = 10, mask = NULL)
