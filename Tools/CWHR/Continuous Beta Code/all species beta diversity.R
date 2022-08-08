library(rgdal)
library(rgeos)
library(picante)
library(betapart)
library(CommEcol)
# library(raster)
library(terra)


##Copy and paste the following chunk of code. This sets up the betagrid function which you will use to calculate
##beta diversity on a continuous surface using a moving window technique.

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


# For each folder/scenario going into the loop
scenarios<-list.dirs(in.dir, recursive = F)

# out dir where you will export files
out.dir<-"C:/Users/User/Documents/Postdoc/TCSI Blueprint Project/out/"


## This loop rolls through each folder and converts rasters to data type which the betagrid function can work with,
##calculates the betagrid function, and exports final surfaces to rasters.

#TODO Make a polygonal grid
#we only need to do this once, outside the loop. Values are calculated inside the loop

# poly <- 

i <- 1
for (i in 1:length(scenarios)){
  #grab all files in a folder
  grids <- list.files(scenarios[i] , pattern = "*.tif$")
  
  #create a raster stack from the input raster files 
  s <- rast(paste0(scenarios[i],"/", grids))
  
  #reclassify suitability <= 0.2 as 0, > 0.2 as 1. SF changed 0.3 to 0.2 to avoid values being ignored
  my_stack <- classify(s, rbind(c(0,0.2,0),  c(0.2,1,1)), right = TRUE)
  
    
  occurrences<-as.data.frame(my_stack,xy=TRUE) 
  occurrences$id<-1:nrow(occurrences)
  # nrow(occurrences)
  # occurrences$longitude<-occurrences$x
  # occurrences$latitude<-occurrences$y
  df <- subset(occurrences, select=c(x,y,id))
  b <- terra::rast(df, crs = crs(my_stack))
  poly<- as.polygons(b,dissolve=FALSE) #SF this is very slow ┤-- can it be avoided?
  newdata <- occurrences[c(3:205)]


  # poly2<-poly
  # crs(poly2)<-""
  # crs(poly2)
  # crs(poly2) <- "+proj=longlat +datum=NAD83 +units=m +no_defs"
  
  ##function to calculate beta diversity using moving window with a radius of 540 m and using sorensen index
  results <- betagrid(gridshp=poly2, comp=newdata, xfeature=2, yfeature=3, radius=540,index="sorensen")
  
  summary(results)
  
  ##The next lines of code extract results from the betagrid function and overlay them on a polygon surface, rasterize
  ##them, and then export the raster surfaces
  poly2$betadiv <- results[,2]
  emptyraster <- raster(extent(poly2))
  res(emptyraster)=180
  beta <- rasterize(poly2, field="betadiv", emptyraster) #SF TODO this is sooooo slow
  
    ##turnover surface
  writeRaster(beta, paste0(out.dir,gsub(in.dir,"",scenarios[i]),"_all species_540_turnover"), format="GTiff", overwrite=TRUE)
  
  poly2$betadiv <- results[,3]
  emptyraster <- raster(extent(poly2))
  res(emptyraster)=180
  beta <- rasterize(poly2, field="betadiv", emptyraster)
  
    
    ##nestedness surface
    writeRaster(beta, paste0(out.dir,gsub(in.dir,"",scenarios[i]),"_all species_540_nestedness"), format="GTiff", overwrite=TRUE)
  
  poly2$betadiv <- results[,4]
  emptyraster <- raster(extent(poly2))
  res(emptyraster)=180
  beta <- rasterize(poly2, field="betadiv", emptyraster)
  
    
    ##mean beta diversity surface
    writeRaster(beta, paste0(out.dir,gsub(in.dir,"",scenarios[i]),"_all species_540_meanbeta"), format="GTiff", overwrite=TRUE)

}


