#dsmartR

###Based on the C5 tree modelling the function determines:
#1. The pixel based probability to each class
#2. The n most probable classes at each pixel
##All outputs are save to file in raster format (It is preferable is the working directory is set to the )
# Class probability maps
# n most probable maps
#R pbjects are saved:
#1. The probability rasters (rasterStack) is requested
#2. The n most probable classes (raster Stack)
#FUnction Requires:
# rLocs: A rasterStack of the rasters generated from the dsmart algorithm
# nprob: the n most probable class maps to be produced
# sepP: logical of whether class probability maps should be produced
# lookup: lookup table produced from dsmart that numerically links soil class codes to a number

dsmartR<- function(rLocs= NULL, nprob = NULL, sepP=FALSE, lookup= NULL, cpus=1){
  beginCluster(cpus)
  #setwd(rLocs)
  dir.create("counts/",showWarnings = F)
  dir.create("probabilities/",showWarnings = F)
  dir.create("nProbable/",showWarnings = F)
  strc<- paste(getwd(),"/counts/",sep="")
  strp<- paste(getwd(),"/probabilities/",sep="")
  strn<- paste(getwd(),"/nProbable/",sep="")
  
  s1<- rLocs
  param<- nrow(lookup)
  param2<-nlayers(s1)
  
  #counts
  nme1<- paste(strc,"countOuts.tif" ,sep="") 
  f1<- function(x) {
    tabulate(x, nbins=param)}
  assign("param", param, envir=.GlobalEnv)
  counts<-clusterR(s1, calc, args=list(fun= f1), export = "param",filename=nme1,format="GTiff",overwrite=T)
  
  
  #probabilities
  nme2<- paste(strp,"countOutsPropbs.tif" ,sep="") 
  param2 = param2
  f2<- function(x) (x/param2)
  assign("param2", param2, envir=.GlobalEnv)
  probs= clusterR(counts, calc,  args=list(fun=f2), export= "param2",filename=nme2,format="GTiff",overwrite=T )
  
  if (sepP==TRUE) {s3<- stack()
                   for (np in 1:nlayers(probs)){
                     nme5<- paste(paste(strp, as.character(lookup[np,1]), sep=""), "_probs.tif", sep="")
                     names(probs[[np]])<- as.character(lookup[np,1])
                     s3<- stack(s3,probs[[np]])
                     writeRaster(probs[[np]],filename=nme5,format="GTiff",overwrite=T)}}
  
  #Most probable
  nme3<- paste(strn,"nProbable.tif",sep="")
  f3<- function(x) order(x, decreasing=TRUE, na.last=TRUE)
  ordered.indices= clusterR(counts, calc,  args=list(fun=f3), filename=nme3,format="GTiff",overwrite=T )
  s4<- stack()
  for (zz in 1:nprob){
    nme4<- paste(strn,paste(zz,"_probable.tif",sep=""),sep="")
    s4<- stack(s4,ordered.indices[[zz]])
    writeRaster(ordered.indices[[zz]],filename=nme4,format="GTiff",overwrite=T)}
  endCluster()
  if (sepP==TRUE){retval<-list(s3,s4)} else{ retval<-list(s4)}
  message(paste("DSMART outputs can be located at:",getwd(), sep=" "))
  return(retval)}

#END