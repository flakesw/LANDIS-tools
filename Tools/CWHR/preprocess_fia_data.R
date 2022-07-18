# preprocess FIA data
# inputs:
#   location of FIA data
#   which states to use
#   what area to trim the included FIA data to?

library(tidyverse)
library(sf)

#the data are stored in one folder
fia_folder <- c("D:/data/fia/rFIA_downloads/")
states <- c("CA", "NV")
boundary_loc <- "C:/Users/Sam/Documents/Research/TCSI conservation finance/Models/Inputs/masks_boundaries/WIP_Capacity_V1Draft/WIP_Capacity_V1Draft.shp"


forest_ref <- read.csv("D:/Data/fia/FIADB_REFERENCE/REF_FOREST_TYPE.csv")  

boundary <-  sf::st_read(boundary_loc) %>%
  summarise(do.union = TRUE) %>% 
  st_make_valid() %>%
  smoothr::fill_holes(threshold = 0.5) 

all_fia_plot <- paste0(fia_folder, states,"_PLOT.csv") %>%
  purrr::map(read.csv) %>%
  dplyr::bind_rows() %>% #combine states into one big list
  dplyr::filter(!is.na(LAT) & !is.na(LON)) %>%
  sf::st_as_sf(coords = c("LON", "LAT"), remove = TRUE, crs = "EPSG:4269") %>% #make spatial
  sf::st_transform(crs = st_crs(boundary))

  
all_fia_cond <- paste0(fia_folder, states,"_COND.csv") %>%
  purrr::map(read.csv) %>%
  do.call(rbind, .) #bind_rows() doesn't work here because of some columns being cast to different data types
  
fia_cond_reduced <- all_fia_cond %>%
  filter(PLT_CN %in% all_fia_plot$CN) %>% #match to FIA plots selected above
  filter(OWNCD != 46) %>%
  mutate(IS_FOREST = ifelse(FORTYPCD %in%(c(1:998)), 1, 0)) %>%
  group_by(PLT_CN) %>%
  summarise(total_cond = sum(CONDPROP_UNADJ),
            natural = sum(STDORGCD, na.rm = TRUE),
            treatment = sum(TRTCD1, na.rm = TRUE),
            proportion_forest = sum(CONDPROP_UNADJ * IS_FOREST),
            cc = sum(CONDPROP_UNADJ * LIVE_CANOPY_CVR_PCT, na.rm = TRUE)) %>%
  filter(total_cond > 0.95 ) %>%
  filter(proportion_forest > 0)

plot_forest_type <- all_fia_cond %>%
  group_by(PLT_CN, FORTYPCD) %>%
  summarise(total_fortypcd = sum(CONDPROP_UNADJ)) %>%
  slice_max(total_fortypcd)

fia_cond_reduced <- left_join(fia_cond_reduced, plot_forest_type, by = "PLT_CN")
fia_cond_reduced$forest_group <- forest_ref[match(fia_cond_reduced$FORTYPCD, forest_ref$VALUE), "TYPGRPCD"]

fia_plot_reduced <- all_fia_plot %>% 
  left_join(fia_cond_reduced, by = c("CN" = "PLT_CN")) %>%  
  sf::st_join(boundary, join = st_within) %>% # subset only plots within boundary
  dplyr::filter(do.union == TRUE) %>%
  sf::st_drop_geometry() 

#----
#Trees
fia_trees <-  paste0(fia_folder, states,"_TREE.csv") %>%
  purrr::map(read.csv) %>%
  dplyr::bind_rows() %>%
  dplyr::filter(PLT_CN %in% fia_plot_reduced$CN)

breaks <- seq(0, max(fia_trees$TOTAGE, na.rm = TRUE) + (10 - max(fia_trees$TOTAGE, na.rm = TRUE) %% 10), by = 5)

fia_trees <- fia_trees%>%
  filter(STATUSCD == 1) %>%
  mutate(DRYBIO_TOTAL = CARBON_AG * 2,
         AGE_BIN = as.numeric(base::cut(TOTAGE, breaks)),
         SPCD = as.character(SPCD))

#------------------
# write outputs

write.csv(fia_plot_reduced, "fia_plot.csv")
write.csv(fia_trees, "fia_trees.csv")
