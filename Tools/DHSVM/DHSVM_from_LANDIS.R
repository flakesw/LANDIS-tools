# This script takes a LANDIS model run and exports the DHSVM veg type,
# canopy cover, LAI by overstory and understory

# This process requires three steps:
#   1. Classify forest type
#   2. Estimate stand seral stage
#   3. Estimate canopy cover
#   4. Estimate shrub cover to decide where understory exists

# Inputs:
# Folder with landis model run
# timestep
# template (mask) raster
# .RDS file with model to predict CC
# .RDS file with models to convert AGE -> DIA
# .RDS file with models to convert AGE -> HT
# .RDS file to model to predict shrub cover

# Outputs:
# Forest type raster
# Stand seral stage raster
# Canopy cover raster
 

library("tidyverse")
library("terra")
library("earth") #for the predictions
library("betareg") #for estimating shrub cover

source("LANDIS_DHSVM_functions.R")

#where are the regressions to relate age to diameter and biomass to canopy cover?
can_regresison_rds_loc <- "canopy_cover_with_CWR_type_lm.RDS"
can_regression_rds_no_sp_loc <- "canopy_cover_without_CWR_type_lm.RDS"
dia_regression_rds_loc = "linear_models_diam_from_age.RDS"
dia_regression_rds_no_sp_loc = "linear_models_diam_from_age_no_sp.RDS"
shrub_regression_loc = "beta_model_shrub_cover.RDS"

#What models and timesteps to use?
timesteps <- seq(0, 80, by = 10)
landis_folders <- list.dirs("E:/TCSI LANDIS/LANDIS runs", recursive = FALSE)
landis_folders <- paste0(landis_folders[grep("Scenario", landis_folders)], "/")[-c(97:102, 16, 17)]

template <- terra::rast("C:/Users/swflake/Documents/TCSI-conservation-finance/Models/Inputs/masks_boundaries/mask_9311.tif")

#where should the outpuits go?
output_folder <- "E:/TCSI LANDIS/Outputs_to_DHSVM_new/"

#Combinations of timesteps and landis runs
input_mods_times <- expand.grid(landis_folders, timesteps) %>%
  rename(landis_folder = Var1,
         timestep = Var2)

error_list <- data.frame(landis_folder = character(0), timestep = numeric(0), iter = numeric(0))

#process DHSVM layers for all model runs
for(i in 1:nrow(input_mods_times)){
  timestep <- input_mods_times[i, "timestep"]
  landis_folder <- as.character(input_mods_times[i, "landis_folder"])
  
  #how should the outputs be named? Timestep and layer name are automatically added
  scenario_name <- paste0(strsplit(landis_folder, "/")[[1]][[4]], " - ")
  
  
  error_flag <- FALSE
  
  start <- Sys.time()
  tryCatch(
    {
      process_DHSVM_layers_and_write_rasters(landis_folder = landis_folder,
                                     template = template,
                                     timestep = timestep,
                                     output_folder = output_folder,
                                     prefix = scenario_name,
                                     dia_regression_rds_loc = dia_regression_rds_loc,
                                     dia_regression_rds_no_sp_loc = dia_regression_rds_no_sp_loc,
                                     can_regresison_rds_loc = can_regresison_rds_loc,
                                     can_regression_rds_no_sp_loc = can_regression_rds_no_sp_loc,
                                     shrub_regression_loc = shrub_regression_loc,
                                     class = FALSE)
    },
    error = function(cond) {
      
      message(paste("Error processing inputs for ", landis_folder," timestep ", timestep))
      message("Here's the original error message:")
      message(cond)
      error_flag <<- TRUE
    }
  )
  
  if(error_flag){
    error_list <- rbind(error_list, cbind(landis_folder, timestep, i))
    next()
  } 
  
  
  end <- Sys.time()
  end - start
  
}
