spathial_way <-
function(X, boundaries, NC, prefiltering){
  if(prefiltering){
    ### Prefilter the data (function pp.rkm_prefilter)
    prefiltered<-rkm_prefilter(X,boundary_ids,plot_ax=TRUE)
    X<-prefiltered$X_filtered
    boundary_ids<-prefiltered$boundary_ids_filtered
    X_g<-prefiltered$X_garbage
    rm(prefiltered)
  }

  ### Initialize waypoints
  waypoint_ids<-initMedoids(X, NC, 'kpp', boundary_ids)
  waypoint_ids<-c(boundary_ids[1],waypoint_ids,boundary_ids[2])
  init_W<-X[waypoint_ids,]

  ### Annealing with rkm
  s_span<-pracma::logspace(5,-5,n=NC)
  s_span<-c(s_span,0)
  #models<-array(data=NA,dim=c(length(s_span),NC+2,ncol(X)))
  #s<-s_span[1]

  models<-list()
  for(i in 1:length(s_span)){
    s<-s_span[i]
    W<-rkm(X,init_W,s,plot_ax=TRUE)
    init_W<-W
    models[[as.character(s)]]<-W
    #models[i,,]<-W
  }
  W_dst_var <- rkm_MS_pathvar(models, s_span, X)
  s_elb_id <- find_elbow(cbind(s_span, W_dst_var))
  return(models[[s_elb_id]])
}
