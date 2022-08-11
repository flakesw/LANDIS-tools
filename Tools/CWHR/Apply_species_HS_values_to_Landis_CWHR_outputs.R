########code to produce species-specific suitability values based on CHWR maps for TCSI
######this file is specifically to send to Kathy Zeller to manipulate based on Tahoe West code/analysis

# Revised by Sam Flake 8/10/2022 for new TCSI LANDIS model runs

library(raster)
library(tidyr)

## in/out dirs

filepath1 <- "E:/TCSI LANDIS/CWHR_outputs2/" ###file were cwhr tiffs live
out.dir <- "E:/TCSI LANDIS/CWHR_outputs_HS2/"###empty file for output

###TCSI_SPP_suitability... csv with habitat suitability values for all species organized by cover type
coefs <- read.csv("TCSI_Spp_suitability_values_1_28_2021.csv", sep=",", header=TRUE)

###create a list of each scenario time step map
veg.files <- list.files(path=filepath1, pattern= "tif$") %>%
  subset(grepl("CWHR_ID", .))

###gets list of species
species<-names(coefs[-c(1:2)])



library(doSNOW)

cl <- makeCluster(4)
registerDoSNOW(cl)
iterations <- length(species)
pb <- txtProgressBar(max = iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

foreach(i = 1:length(species),
        .combine = "c",
        .verbose=T, 
        .packages = "raster",
        .options.snow = opts) %dopar% {

  
  sp <- species[i]
  
  temp.lib<-as.data.frame(cbind(coefs$CWHR_ID, coefs[,sp]))
  
  for (k in veg.files){
    ###base model
    
    temp.rast <-  raster(paste0(filepath1,k)) 
    name.rast<-gsub("\\..*","",k)
    name.rast <- gsub("CWHR_ID-", "", name.rast)
    name.rast <- gsub(" ", "", name.rast)
    # # replace any raster value <100 with NA
    # values(temp.rast)[values(temp.rast<100)]<-NA
    # # replace any raster value >1400 with 1400
    # values(temp.rast)[values(temp.rast>1400)]<-1400
    
    base<-subs(temp.rast, temp.lib, by=1, which=2)
    if(!dir.exists(paste0(out.dir,name.rast))){
     dir.create(paste0(out.dir,name.rast))
    }
    
    writeRaster(base, paste0(out.dir,name.rast,"/", sp, "_",name.rast,'.tif') ,overwrite=TRUE)
    
  }
  
  return(NULL)
}


