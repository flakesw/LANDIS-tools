library(rgdal)
library(rgeos)
library(picante)
library(betapart)
library(CommEcol)
library(raster)
library(terra)
library(dplyr)
library(sf)


##Copy and paste beta diversity function. This chunk of code sets up the function to be used later on.

betagrid<-function(gridshp, comp, xfeature, yfeature, radius, phylotree, phylobeta=F, index="sorensen"){
  data<-data.frame(gridshp[xfeature],gridshp[yfeature],comp)
  mean_turnover<-numeric(length(comp[,1]))
  mean_nestedness<-numeric(length(comp[,1]))
  mean_beta<-numeric(length(comp[,1]))
  for(i in 1:length(gridshp[[2]])){
    adj<-select.window(xf=data[i,1], yf=data[i,2], radius, xydata=data)[,-c(1,2)]
    if(phylobeta==F){
      ifelse(sum(nrow(adj))==0 || ncol(adj)==0, res<-0 , res<-beta.pair(adj, index.family=index))
    }else if(phylobeta==T){
      ifelse(sum(nrow(adj))==0 || ncol(adj)==0, res<-0 , res<-phylo.beta.pair(adj, phylotree, index.family=index))
    }
    ifelse(sum(nrow(adj))==0 || ncol(adj)==0, mean_turnover[i]<-0 , mean_turnover[i]<-mean(as.matrix(res[[1]])[2:length(as.matrix(res[[1]])[,1]),1],na.rm=TRUE) )
    ifelse(sum(nrow(adj))==0 || ncol(adj)==0, mean_nestedness[i]<-0 , mean_nestedness[i]<-mean(as.matrix(res[[2]])[2:length(as.matrix(res[[2]])[,1]),1],na.rm=TRUE) )
    ifelse(sum(nrow(adj))==0 || ncol(adj)==0, mean_beta[i]<-0 , mean_beta[i]<-mean(as.matrix(res[[3]])[2:length(as.matrix(res[[3]])[,1]),1],na.rm=TRUE) )
  }
  return(data.frame(cell=row.names(comp), mean_turnover, mean_nestedness, mean_beta))
}


## Assign input directory. This directory leads to all the CWHR suitability rasters for all species for all time steps
in.dir<- "E:/TCSI LANDIS/CWHR_outputs_HS/"


# read in csv with species codes for functional groups
fnlg<-read.csv("./Continuous Beta Code/Functional_Groups_v3.csv")

# For each folder/scenario going into the loop
scenarios<-list.dirs(in.dir, recursive = F)

# Assign output directory
out.dir<-"E:/TCSI LANDIS/CWHR_diversity"



##This loop rolls through each folder and modifies suitability rasters for all species to a format
##which the betagrid function can work with and sets up a polygon to write results into
for (i in 1:length(scenarios)){
  #grab all files in a folder
  grids <- list.files(scenarios[i] , pattern = "*.tif$")
  
  #create a raster stack from the input raster files 
  s <- rast(paste0(scenarios[i],"/", grids))
  my_stack <- classify(s, rbind(c(0,0.2,0),  c(0.3,1,1)))

  # binarize each raster in stack and subset to seed dispersers
  fg <- fnlg$Cavity_excavators 
  fg <- fg[fg!= ""]
  
  cgridstack <- raster::subset(my_stack, grep(paste(fg, collapse="|"), grids))
  
  ##Convert raster stack to dataframe for betagrid function and create polygon to write results to  
  occurrences<-as.data.frame(cgridstack,xy=TRUE)
  occurrences$id<-1:nrow(occurrences)
  nrow(occurrences)
  occurrences$longitude<-occurrences$x
  occurrences$latitude<-occurrences$y
  df <- subset(occurrences, select=c(x,y,id,longitude,latitude))
  b <- raster::rasterFromXYZ(df,res=c(180,180))
  poly<- raster::rasterToPolygons(b,dissolve=FALSE) #should already in the right CRS
  newdata <- occurrences[c(3:14)]
  poly2 <- poly
  
  ##Run betagrid function with 540 m radius moving window
  results <- betagrid(gridshp=poly2, comp=newdata, xfeature=2, yfeature=3, radius=540,index="sorensen")
  
  ##summary(results)
  
  poly2 <- st_as_sf(poly)
  
  ##Extract the turnover result and write/export to a raster
  poly2$betadiv <- results[,2]
  emptyraster <- terra::rast(ext(vect(poly2)), resolution = 180, crs = crs(poly2))
  beta <- terra::rasterize(vect(poly2), field="betadiv", emptyraster) #terra is much faster at rasterizing

  terra::writeRaster(beta, paste0(out.dir, gsub(in.dir,"",scenarios[i]),"_cavityexcavators_540_turnover"), filetype="GTiff")
  
  ##Extract the nestedness result and write/export to a raster
  poly2$betadiv <- results[,3]
  emptyraster <- terra::rast(ext(vect(poly2)), resolution = 180, crs = crs(poly2))
  beta <- terra::rasterize(vect(poly2), field="betadiv", emptyraster) #terra is much faster at rasterizing
  beta <- rasterize(poly2, field="betadiv", emptyraster)
  
  terra::writeRaster(beta, paste0(out.dir,gsub(in.dir,"",scenarios[i]),"_cavityexcavators_540_nestedness"), filetype="GTiff")
  
  ##Extract the mean beta result and write/export to a raster
  poly2$betadiv <- results[,4]
  emptyraster <- terra::rast(ext(vect(poly2)), resolution = 180, crs = crs(poly2))
  beta <- terra::rasterize(vect(poly2), field="betadiv", emptyraster) #terra is much faster at rasterizing
  
  terra::writeRaster(beta, paste0(out.dir,gsub(in.dir,"",scenarios[i]),"_cavityexcavators_540_meanbeta"), filetype="GTiff")

}

