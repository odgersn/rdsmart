library(rgdal); library(raster); library(sp); library(gtools); library(C50)


#Files location
setwd("/home/brendo/myWork/dsmart/data")
# Identify map unit polygon shapefile (attribute table structured the same way as for Python DSMART)
shapefile = "dlr_polys_alb_7km.shp"
# Map unit composition file
compositionFile = "formatted_attribute_table2.txt"


# Load data
# Load shapefile
shp = readOGR(dsn=paste0("./polygons/", shapefile), layer=substr(shapefile, 1, nchar(shapefile) - 4))
polygons<- shp

# Load covariates as raster stack
covariates=stack()
for(filename in dir(path=paste0(getwd(),"/covariates/"), pattern=".img$"))
{
  r = raster(paste0(getwd(),"/covariates/",filename))
  # Load raster to stack
  covariates=stack(covariates, raster(paste0(getwd(),"/covariates/",filename)))
}

# Load map unit composition text file
composition = read.table(paste0(getwd(), "/polygons/", compositionFile), sep=",", header=TRUE)
colnames(composition) = c("poly","mapunit","soil_class","proportion")


####### DSMART tree creation #######
dsmart(covariates = covariates, polygons = shp, composition = composition, n=15, reals = 20, cpus=8)

####### DSMART probability rasters #######
setwd("/home/brendo/myWork/dsmart/data/dsmartOuts/rasters")
list.files(getwd(),  pattern="tif$", full.names=FALSE)
files<- list.files(getwd(), pattern='tif$',full.names=T)

s8<- stack()
for (i in 1:length(files)){
  r1<- raster(files[i])
  s8<- stack(s8, r1)}
lookup <- read.table("classLookupTable.txt",sep=",", header=TRUE)

test1<- dsmartR(rLocs= s8, nprob = 3, sepP=TRUE, lookup= lookup,cpus=2)



