library("raster")
library("tidyverse")

in_folder <- "E:/TCSI LANDIS/CWHR/CWHR_diversity_new/"
out_folder <- "E:/TCSI LANDIS/CWHR/CWHR_diversity_new_reprojected/"

project_to_template <- function(input_raster, template){
  #function to project landis input rasters with no CRS to an input file
  
  #replace values of template with values from input raster
  values(template) <- (values(input_raster))
  
  
  return(template)
}


template <- raster("mask_9311.tif")

raster_list <- list.files(in_folder, pattern = ".tif")
# raster_list <- raster_list[!grepl("CWHR_ID", raster_list)]
raster_list <- raster_list[extension(raster_list) %in% c(".img", ".tif")]

rasters_stripped <- sub('\\..*$', '', basename(raster_list))



library(doSNOW)

cl <- makeCluster(4)
registerDoSNOW(cl)
iterations <- length(raster_list)
pb <- txtProgressBar(max = iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

foreach(i = 1:length(raster_list),
        .combine = "c",
        .verbose=T, 
        .packages = "raster",
        .options.snow = opts) %dopar% {

  oldrast <- raster(paste0(in_folder, raster_list[i]))
  newrast <- project_to_template(oldrast, template)
  writeRaster(newrast, 
              paste0(out_folder,rasters_stripped[i], ".tif"),
              filetype = "GEOTiff",
              datatype = dataType(oldrast),
              overwrite = TRUE)
  
  return(NULL)
}

