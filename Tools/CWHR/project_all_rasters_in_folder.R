library("terra")
library("tidyverse")

project_to_template <- function(input_raster, template){
  #function to project landis input rasters with no CRS to an input file
  
  if(ncol(input_raster) == ncol(template) & nrow(input_raster) == nrow(template)){
    #replace values of template with values from input raster
    out_raster <- template %>%
      `values<-`(values(input_raster))
  } else print("Rasters have different numbers of rows or cols")
  
  # plot(out_raster)
  
  return(out_raster)
}

input_folder <- "E:/TCSI LANDIS/Outputs_to_DHSVM_new/" #folder containing all the rasters to project
output_folder <- input_folder #change this if you want to!

if(!dir.exists(paste0(output_folder, "projected"))){
  dir.create(paste0(output_folder, "projected"))
}

template <- rast("C:/Users/swflake/Documents/TCSI-conservation-finance/Models/Inputs/masks_boundaries/mask_9311.tif") #raster to replace values of

raster_list <- list.files(input_folder)
raster_list <- raster_list[grepl(".img", raster_list) | grepl(".tif", raster_list)]

rasters_stripped <- sub('//..*$', '', basename(raster_list))

for(i in 1:length(raster_list)){
  oldrast <- suppressWarnings(rast(paste0(input_folder, raster_list[i])))
  newrast <- project_to_template(oldrast, template) %>%
    terra::project("epsg:32610", method = "near") #CHANGE THIS for expected output
  terra::writeRaster(newrast, 
                     paste0(output_folder, "projected/",rasters_stripped[i]),
                     filetype = "GTiff",
                     datatype = ifelse(is.int(oldrast), "INT2S", "FLT4S"),
                     overwrite = TRUE)
}






#old version, need to update
# library("raster")
# library("tidyverse")
# 
# in_folder <- "E:/TCSI LANDIS/CWHR/CWHR_diversity_new/"
# out_folder <- "E:/TCSI LANDIS/CWHR/CWHR_diversity_new_reprojected/"
# 
# project_to_template <- function(input_raster, template){
#   #function to project landis input rasters with no CRS to an input file
#   
#   #replace values of template with values from input raster
#   values(template) <- (values(input_raster))
#   
#   
#   return(template)
# }
# 
# 
# template <- raster("mask_9311.tif")
# 
# raster_list <- list.files(in_folder, pattern = ".tif")
# # raster_list <- raster_list[!grepl("CWHR_ID", raster_list)]
# raster_list <- raster_list[extension(raster_list) %in% c(".img", ".tif")]
# 
# rasters_stripped <- sub('\\..*$', '', basename(raster_list))
# 
# 
# 
# library(doSNOW)
# 
# cl <- makeCluster(4)
# registerDoSNOW(cl)
# iterations <- length(raster_list)
# pb <- txtProgressBar(max = iterations, style = 3)
# progress <- function(n) setTxtProgressBar(pb, n)
# opts <- list(progress = progress)
# 
# foreach(i = 1:length(raster_list),
#         .combine = "c",
#         .verbose=T, 
#         .packages = "raster",
#         .options.snow = opts) %dopar% {
# 
#   oldrast <- raster(paste0(in_folder, raster_list[i]))
#   newrast <- project_to_template(oldrast, template)
#   writeRaster(newrast, 
#               paste0(out_folder,rasters_stripped[i], ".tif"),
#               filetype = "GEOTiff",
#               datatype = dataType(oldrast),
#               overwrite = TRUE)
#   
#   return(NULL)
# }

