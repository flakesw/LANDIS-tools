#process cwhr id rasters and test which are not in reference
library(terra)
library(tidyverse)

ref <- read.csv("TCSI_Spp_suitability_values_1_28_2021.csv")

cwhr_rasters <- list.files("E:/TCSI LANDIS/CWHR_outputs2/", pattern = "CWHR")

cwhr_ids <- vector("integer")

for(i in 1:length(cwhr_rasters)){
  raster <- terra::rast(paste0("E:/TCSI LANDIS/CWHR_outputs2/", cwhr_rasters[i]))
  table(terra::values(raster))
  cwhr_ids <- c(cwhr_ids, terra::values(raster))
}

tab <- table(cwhr_ids)

df <- data.frame(cwhr_id = tab,
                 in_ref = NA) %>%
  rename(cwhr_id = cwhr_id.cwhr_ids,
         freq = cwhr_id.Freq)

df$in_ref <- df$cwhr_id %in% ref$CWHR_ID

write.csv(df, "cwhr_codes_in_ref.csv")
