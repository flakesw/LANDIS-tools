library(picante)
library(betapart)
library(CommEcol)
library(terra)
library(dplyr)
library(doSNOW)

##Copy and paste beta diversity function. This chunk of code sets up the function to be used later on.

betagrid<-function(gridstack, radius, phylotree, phylobeta=F, index="sorensen"){
  
  crds <- terra::crds(gridstack)
  xydata = terra::values(gridstack, dataframe = TRUE) %>%
    filter(complete.cases(.)) %>%
    cbind(crds, .)
  cells <- cellFromXY(gridstack, crds)
  
  mean_turnover<-numeric(length(crds[,1]))
  mean_nestedness<-numeric(length(crds[,1]))
  mean_beta<-numeric(length(crds[,1]))
  
  for(i in 1:nrow(xydata)){
    adj <- select.window(xf= crds[i, 1], yf = crds[i, 2], radius, xydata=xydata)[,-c(1,2)]
    
    if(phylobeta==F){
      ifelse(sum(nrow(adj))==0 || ncol(adj)==0, res<-0 , res<-beta.pair(adj, index.family=index))
      
    }else if(phylobeta==T){
      ifelse(sum(nrow(adj))==0 || ncol(adj)==0, res<-0 , res<-phylo.beta.pair(adj, phylotree, index.family=index))
    }
    
    ifelse(sum(nrow(adj))==0 || ncol(adj)==0, mean_turnover[i]<-0, mean_turnover[i]<-mean(as.matrix(res[[1]])[2:length(as.matrix(res[[1]])[,1]),1],na.rm=TRUE) )
    ifelse(sum(nrow(adj))==0 || ncol(adj)==0, mean_nestedness[i]<-0 , mean_nestedness[i]<-mean(as.matrix(res[[2]])[2:length(as.matrix(res[[2]])[,1]),1],na.rm=TRUE) )
    ifelse(sum(nrow(adj))==0 || ncol(adj)==0, mean_beta[i]<-0 , mean_beta[i]<-mean(as.matrix(res[[3]])[2:length(as.matrix(res[[3]])[,1]),1],na.rm=TRUE) )
  }
  results <- data.frame(cell=cells, mean_turnover, mean_nestedness, mean_beta)
  return(data.frame(cell=cells, mean_turnover, mean_nestedness, mean_beta))
}


## Assign input directory. This directory leads to all the CWHR suitability rasters for all species for all time steps
in.dir<- "E:/TCSI LANDIS/CWHR/CWHR_outputs_species_suitability_new/"

# read in csv with species codes for functional groups
fnlg<-read.csv("./Continuous Beta Code/Functional_Groups_v3.csv")

# For each folder/scenario going into the loop
scenarios<-list.dirs(in.dir, recursive = F)

# Assign output directory
out.dir<-"E:/TCSI LANDIS/CWHR/CWHR_diversity_new"

# scenarios <- scenarios[c(281, 289, 264) ]

cl <- makeCluster(3)
registerDoSNOW(cl)
iterations <- length(scenarios)
pb <- txtProgressBar(max = iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

foreach(i = 1:length(scenarios),
        .combine = "c",
        .verbose=T, 
        .packages = c("terra", "betapart", "CommEcol", "dplyr"),
        .options.snow = opts) %dopar% {
  ##This loop rolls through each folder and modifies suitability rasters for all species to a format
  ##which the betagrid function can work with and sets up a polygon to write results into

  # this takes 40 minutes per iteration
          
  #grab all files in a folder
  grids <- list.files(scenarios[i] , pattern = "*.tif$")
  
  # binarize each raster in stack and subset to predators
  fg <- fnlg$Predators 
  fg <- fg[fg != ""]
  
  #create a raster stack from the input raster files 
  cgridstack <- rast(paste0(scenarios[i],"/", grids)) %>%
    classify(rbind(c(0,0.2,0),  c(0.3,1,1))) %>%
    raster::subset(grep(paste(fg, collapse="|"), grids)) 
  
  ##Run betagrid function with 540 m radius moving window
  results <- betagrid(gridstack = cgridstack, radius=540,index="sorensen")

  
  emptyraster <- terra::rast(ext(cgridstack), resolution = 180, crs = crs(cgridstack))
  
  ##Extract the turnover result and write/export to a raster
  r <- emptyraster
  r[results$cell] <- results$mean_turnover
  terra::writeRaster(r, paste0(out.dir, gsub(in.dir,"",scenarios[i]),"_predators_540_turnover.tif"), filetype="GTiff")
  
  ##Extract the turnover result and write/export to a raster
  r <- emptyraster
  r[results$cell] <- results$mean_nestedness
  terra::writeRaster(r, paste0(out.dir, gsub(in.dir,"",scenarios[i]),"_predators_540_nestedness.tif"), filetype="GTiff")

  
  ##Extract the mean beta result and write/export to a raster
  r <- emptyraster
  r[results$cell] <- results$mean_beta
  terra::writeRaster(r, paste0(out.dir, gsub(in.dir,"",scenarios[i]),"_predators_540_beta.tif"), filetype="GTiff")
  
  return(NULL)
  
}
