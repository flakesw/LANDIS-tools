# This script takes a LANDIS model run and exports the CWHR classification

# This process requires three steps:
#   1. Classify forest type
#   2. Estimate stand seral stage
#   3. Estimate canopy cover

# Inputs:
# Folder with landis model run
# timestep
# template (mask) raster
# .RDS file with model to predict CC

# Outputs:
# Forest type raster
# Stand seral stage raster
# Canopy cover raster
# 

timestep <- ttt

landis_folder <- xxx

biomass_folder <- yyy

community_table <- ccc

community_raster <- rrr

#*******************************************************************************
#Forest type classification
#*******************************************************************************

forest_raster <- classify_forest_type(biomass_folder, timestep)
writeRaster()

#*******************************************************************************
# Stand seral stage
#*******************************************************************************







