#Project LANDIS output rasters to original projection of input files

#Sam Flake, 2023-3-10
#swflake@ncsu.edu

#This script takes a folder full of output rasters, e.g. the outputs from NECN,
# which are in the native LANDIS crs (i.e., no CRS), and writes a new set
# of rasters with a projection matching your input files. Optionally, 
# it can also reproject to a new CRS. By default, it will create a new folder
# inside your input folder (i.e., if "./Landis run/NECN/" has your input files,
# it will create "./Landis run/NECN/projected/"). 


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

#where are the rasters stored?
input_folder <- 'C:/Users/Sam/Downloads/initial layers/'
output_folder <- input_folder #change this if you want to

#if you want to reproject after translating from LANDIS -> template, set reproject = TRUE
# and put the proper crs (e.g. "epsg:26910") for new_crs
reproject <- FALSE
new_crs <- ""

#if input folder and output folder are the same, create a subfolder to store the projected outputs
if(!dir.exists(paste0(output_folder, "projected"))){
  dir.create(paste0(output_folder, "projected"))
}

#the template raster can be any of your input rasters. It should have a CRS, and the same number of
#rows and columns as your outputs. Instead of reprojecting our LANDIS outputs,
#we'll just reassign values of the template from the outputs we want to project
# template <- rast("C:/Users/swflake/Documents/TCSI-conservation-finance/Models/Inputs/masks_boundaries/mask_9311.tif") #raster to replace values of
template <- rast("C:/Users/Sam/Documents/Research/TCSI conservation finance/Models/Inputs/masks_boundaries/mask_9311_new.tif")

raster_list <- list.files(input_folder)
raster_list <- raster_list[grepl(".img", raster_list) | grepl(".tif", raster_list)] #select .img or .tif files

#strip the directory and extension
rasters_stripped <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(raster_list))

for(i in 1:length(raster_list)){
  oldrast <- suppressWarnings(rast(paste0(input_folder, raster_list[i])))
  newrast <- project_to_template(oldrast, template) 
  if(reproject){
    newrast <- terra::project(newrast, new_crs, method = "near") #"near" is safer in case of categorical rasters
  }
  terra::writeRaster(newrast, 
                     paste0(output_folder, "projected/",rasters_stripped[i], ".tif"),
                     filetype = "GTiff",
                     datatype = ifelse(is.int(oldrast), "INT4S", "FLT4S"),
                     overwrite = TRUE)
}