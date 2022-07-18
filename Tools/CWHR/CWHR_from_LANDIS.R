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
# .RDS file with models to convert AGE -> DIA

# Outputs:
# Forest type raster
# Stand seral stage raster
# Canopy cover raster
# 

library("tidyverse")
library("terra")

source("LANDIS_CWHR_functions.R")

timestep <- 0

landis_folder <- "C:/Users/Sam/My Documents/GlobusEndpoint/"

# biomass_folder <- yyy

comm_table <- paste0(landis_folder, "community-input-file-",timestep, ".csv")

comm_raster <- terra::rast(paste0(landis_folder, "output-community-",timestep, ".img"))

#*******************************************************************************
#Forest type classification
#*******************************************************************************

# this takes about 5 minutes to run. It might take a lot of RAM if you're doing
# a big landscape. A place to improve might be to break this into chunks and do 
# each chunk individually, then mosaic back together.


#TODO update these
type_to_number <- as.data.frame(
          rbind(c("asp", 1),
          c("mhw", 2),
          c("mri", 3),
          c("wfr", 4),
          c("rfr", 5),
          c("jpn", 6),
          c("ppn", 7),
          c("dfr", 8),
          c("mhc", 9),
          c("lpn", 10),
          c("smc", 11),
          c("scn", 12),
          c("juo", 13),
          c("mch", 14)
          )) %>%
           dplyr::rename(type = V1,
                         num = V2)

# mutate(forname =
#          case_when( forname == 0 ~ "Non-forested",
#                     forname == 1 ~ "Aspen",
#                     forname == 2 ~ "Montane Hardwood",
#                     forname == 3 ~ "Montane Riparian",
#                     forname == 4 ~ "White fir",
#                     forname == 5 ~ "Red fir",
#                     forname == 6 ~ "Jeffrey pine",
#                     forname == 7 ~ "Ponderosa pine",
#                     forname == 8 ~ "Douglas-fir",
#                     forname == 9 ~ "Mixed hardwood conifer",
#                     forname == 10 ~ "Lodgepole pine",
#                     forname == 11 ~ "Sierra mixed conifer",
#                     forname == 12 ~ "Sierra high elevation mixed conifer",
#                     forname == 13 ~ "Juniper",
#                     forname == 14 ~ "Chapparal"
#          ))

#TODO fix the NAs in forest type
forest_types <- classify_forest_type_comm_output(comm_table) 

forest_types2 <- ungroup(forest_types )%>%
  dplyr::mutate(MapCode = as.numeric(as.character(MapCode))) %>%
  dplyr::mutate(CWHR_type = as.numeric(type_to_number[match(CWHR_type, type_to_number$type), "num"]))
forest_types2$MapCode <- as.integer(forest_types2$MapCode)
forest_types2$CWHR_type <- as.integer(forest_types2$CWHR_type)

# test <- values(comm_raster[!(values(comm_raster) %in%  forest_types2$MapCode)])

#TODO this is very slow -- speed up?
forest_raster <- terra::classify(comm_raster, rcl = forest_types2) %>%
  terra::clamp(upper = max(as.integer(type_to_number$num)), values = FALSE) 
  #if there are values greater than 22, then they will be assigned as NA

plot(forest_raster)

writeRaster(forest_raster, "test_output_forest_type.tif", overwrite = TRUE)

#*******************************************************************************
# Stand seral stage
#*******************************************************************************
comm <- read.csv(comm_table) %>%
  rename(TOTAGE = CohortAge)
comm[comm$SpeciesName %in% c("FX_R_SEED", "NOFX_R_SEED", "NOFX_NOR_SEED"), "SpeciesName"] <- "Shrub"
# comm_test <- comm[comm$MapCode %in% c(1:100), ]

# TODO check whether adding 10 to the ages is a good idea; it makes the young
# cohorts very small when making predictions

#import regression equations fit from FIA data, see regression_age_diam.R
age_to_dia_reg <- readRDS("linear_models_diam_from_age.RDS")
age_to_dia_reg_no_sp <- readRDS("linear_models_diam_from_age_no_sp.RDS")

# get predictions for diameter from age
# TODO check whether we can include biomass in regressions

sp_in_comm <- unique(comm$SpeciesName)

#this takes about a minute
for(i in 1:length(sp_in_comm)){
  sp = sp_in_comm[i]
  if(sp %in% age_to_dia_reg$species_code){
    comm[comm$SpeciesName == sp, "dia"] <- exp(predict(age_to_dia_reg$model[age_to_dia_reg$species_code == sp][[1]],
                    newdata = comm[comm$SpeciesName == sp, ]))
  } else{
    comm[comm$SpeciesName == sp, "dia"] <- exp(predict(age_to_dia_reg_no_sp, newdata = comm[comm$SpeciesName == sp, ]))
  }
}

comm$dia <- ifelse(comm$dia < 1, 0.1, comm$dia) #increase size of tiny trees; TODO is this right?

#aggregate to plot
plot_age_dia <- comm %>%
  dplyr::group_by(MapCode) %>%
  dplyr::summarise(mean_age = weighted.mean(TOTAGE, CohortBiomass),
                   weighted_mean_diam = weighted.mean(dia, CohortBiomass),
                   qmd = sqrt(weighted.mean((dia)^2, CohortBiomass)),
                   max_dia = max(dia)) #is this right?

plot(plot_age_dia$mean_age ~ plot_age_dia$weighted_mean_diam)
plot(plot_age_dia$mean_age[1:10000] ~ plot_age_dia$qmd[1:10000],
     ylim = c(0, 250),
     xlab = "QMD (inches)",
     ylab = "Mean age (weighted by biomass)")

# some test values from the crosswalk table
# not too far off from crosswalk
abline(v = 6, col = "dark gray")
abline(h = 30, col = "dark gray")
abline(v = 11, col = "dark gray")
abline(h = 70, col = "dark gray")
abline(v = 24, col = "dark gray")
abline(h = 120, col = "dark gray")

text(x = 4, y = 20, labels = "Class 3 start")
text(x = 10, y = 90, labels = "Class 4 start")
text(x = 26, y = 100, labels = "Class 5 start")

plot(plot_age_dia$qmd ~ plot_age_dia$weighted_mean_diam)
plot(plot_age_dia$max_dia ~ plot_age_dia$qmd)

dia_breaks <- c(0,1,6,11,24,1000)
plot_age_dia$seral_class <- cut(plot_age_dia$qmd, dia_breaks, labels = FALSE)

seral_raster <- terra::classify(comm_raster, rcl = plot_age_dia[, c("MapCode", "seral_class")]) %>%
  terra::clamp(upper = 5, values = FALSE) 

plot(seral_raster)

writeRaster(seral_raster, "test_output_seral_class.tif", overwrite = TRUE)

#*******************************************************************************
# Stand canopy cover
#*******************************************************************************

comm <- read.csv(comm_table) %>%
  rename(TOTAGE = CohortAge)
comm[comm$SpeciesName %in% c("FX_R_SEED", "NOFX_R_SEED", "NOFX_NOR_SEED"), "SpeciesName"] <- "Shrub"

#get forest type
forest_types <- classify_forest_type_comm_output(comm_table) 

# get stand age
# TODO
plot_age_dia$mean_age

can_model <- readRDS("canopy_cover_with_CWR_type_lm.RDS")
can_model_no_cwhr <- readRDS("canopy_cover_without_CWR_type_lm.RDS")

plot_biomass <- comm %>%
  group_by(MapCode) %>%
  summarise(total_biomass = sum(CohortBiomass)/100) %>% #convert to Mg/ha
  left_join(forest_types, by = "MapCode") %>%
  left_join(plot_age_dia, by = "MapCode")

cwhr_types_with_mod <- unique(can_model$model$CWHR_type)
all_cwhr_types <- unique(plot_biomass$CWHR_type)

plot_biomass$cc <- NA

for(i in 1:length(all_cwhr_types)){
  cwhr <- all_cwhr_types[i]
  if(cwhr %in% cwhr_types_with_mod){
    plot_biomass[plot_biomass$CWHR_type == cwhr, "cc"] <- invlogit(predict(can_model, 
                                                                  newdata = plot_biomass[plot_biomass$CWHR_type == cwhr, ]))
  } else{
    plot_biomass[plot_biomass$CWHR_type == cwhr, "cc"] <- invlogit(predict(can_model_no_cwhr, 
                                                                           newdata = plot_biomass[plot_biomass$CWHR_type == cwhr, ]))
  }
}

cc_breaks <- c(0,0.1,0.25,0.4,0.6,1)
plot_biomass$cc_class <- cut(plot_biomass$cc, cc_breaks, labels = FALSE)

cc_raster <- terra::classify(comm_raster, rcl = plot_biomass[, c("MapCode", "cc_class")]) %>%
  terra::clamp(lower = 0, upper = 5, values = FALSE) 

plot(cc_raster)

writeRaster(cc_raster, "test_output_can_cov_class.tif", overwrite = TRUE)






