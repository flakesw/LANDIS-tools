# Estimating shrub cover from tree data
# Inputs needed
#   pre-processed FIA data
#   regressions to get age from diameter

# Outputs:
#   Regressions for CC ~ CHWR type + Age + Biomass

library("tidyverse")
source("./LANDIS_DHSVM_functions.R")
library(sf)

logit <- function(x) log(x/(1-x))
invlogit <- function(x) exp(x)/(1+exp(x))


# read in plot and cond data
ca_fia_plot <- read.csv("D:/data/fia/rFIA_downloads/CA_PLOT.csv") %>%
  dplyr::filter(!is.na(LAT) & !is.na(LON)) %>%
  sf::st_as_sf(coords = c("LON", "LAT"), remove = TRUE, crs = "EPSG:4269") %>%
  sf::st_transform(crs = "EPSG:5070")

nv_fia_plot <- read.csv("D:/data/fia/rFIA_downloads/NV_PLOT.csv")%>%
  dplyr::filter(!is.na(LAT) & !is.na(LON)) %>%
  sf::st_as_sf(coords = c("LON", "LAT"), remove = TRUE, crs = "EPSG:4269") %>%
  sf::st_transform(crs = "EPSG:5070")

all_fia_plot <- rbind(ca_fia_plot, nv_fia_plot) %>%
  sf::st_transform(crs = "EPSG:5070")

ca_fia_cond <- read.csv("D:/Data/fia/rFIA_downloads/CA_COND.csv")
nv_fia_cond <- read.csv("D:/data/fia/rFIA_downloads/NV_COND.csv")

forest_ref <- read.csv("D:/Data/fia/FIADB_REFERENCE/REF_FOREST_TYPE.csv")

all_fia_cond <- rbind(ca_fia_cond, nv_fia_cond)

sierra_shape <- sf::st_read("C:/Users/Sam/Documents/Research/TCSI conservation finance/Models/Inputs/masks_boundaries/WIP_Capacity_V1Draft/WIP_Capacity_V1Draft.shp") %>%
  summarise(do.union = TRUE) %>% 
  st_make_valid() %>%
  smoothr::fill_holes(threshold = 0.5) %>%
  st_transform("EPSG:5070")

tcsi_shape <- sf::st_read("C:/Users/Sam/Documents/Research/TCSI conservation finance/Models/Inputs/masks_boundaries/tcsi_area_shapefile/TCSI_v2.shp") %>%
  st_transform("EPSG:5070")

#reduce plots to study area
sierra_fia_plot <- sf::st_join(all_fia_plot, sierra_shape, join = st_within) %>%
  filter(do.union == TRUE)

#plot study area and FIA plots
# plot(sf::st_geometry(sierra_shape))
# plot(sf::st_geometry(sierra_fia_plot), add = TRUE)

sierra_fia_cond <- all_fia_cond %>%
  filter(PLT_CN %in% sierra_fia_plot$CN) %>%
  # filter(OWNCD != 46) %>%
  mutate(IS_FOREST = ifelse(FORTYPCD %in%(c(1:998)), 1, 0)) %>%
  group_by(PLT_CN) %>%
  summarise(total_cond = sum(CONDPROP_UNADJ),
            natural = sum(STDORGCD, na.rm = TRUE),
            treatment = sum(TRTCD1, na.rm = TRUE),
            proportion_forest = sum(CONDPROP_UNADJ * IS_FOREST),
            cc = sum(CONDPROP_UNADJ * LIVE_CANOPY_CVR_PCT, na.rm = TRUE)) %>%
  filter(total_cond > 0.95 )

# plot_forest_type <- all_fia_cond %>%
#   mutate(FORTYPCD = ifelse(FLDSZCD == 1 | COND_STATUS_CD == 2, 1, FORTYPCD)) %>%
#   group_by(PLT_CN, FORTYPCD) %>%
#   summarise(total_fortypcd = sum(CONDPROP_UNADJ)) %>%
#   slice_max(total_fortypcd)

# sierra_fia_cond <- left_join(sierra_fia_cond, plot_forest_type, by = "PLT_CN")
# sierra_fia_cond$forest_group <- forest_ref[match(sierra_fia_cond$FORTYPCD, forest_ref$VALUE), "TYPGRPCD"]

fia_plot <- all_fia_plot %>% 
  left_join(sierra_fia_cond, by = c("CN" = "PLT_CN")) 
  # dplyr::filter(INVYR > 2006)

#----
#Trees
ca_trees <- read.csv("D:/data/fia/rFIA_downloads/CA_TREE.csv") %>%
  dplyr::filter(PLT_CN %in% sierra_fia_plot$CN)

nv_trees <- read.csv("D:/data/fia/rFIA_downloads/NV_TREE.csv") %>%
  dplyr::filter(PLT_CN %in% sierra_fia_plot$CN)

sierra_trees <- rbind(ca_trees, nv_trees)

sierra_trees <- filter(sierra_trees, PLT_CN %in% sierra_fia_plot$CN) 
  
breaks <- seq(0, max(sierra_trees$TOTAGE, na.rm = TRUE) + (10 - max(sierra_trees$TOTAGE, na.rm = TRUE) %% 10), by = 5)

fia_trees <- sierra_trees%>%
  dplyr::mutate(TOTAGE = BHAGE + 10) %>%
  filter(STATUSCD == 1) %>%
  mutate(DRYBIO_TOTAL = CARBON_AG * 2 / 2204, #convert to Mg
         AGE_BIN = as.numeric(base::cut(TOTAGE, breaks)),
         SPCD = as.character(SPCD))



################################################################################
# First step: classify FIA stands into CWHR community types
#TODO update with new classification

#species used for classification
species_class <- c("AbieConc", "AbieMagn", "CaloDecu", "JuniOcci", "PinuAlbi", "PinuBalf", "PinuCont",
             "PinuLamb", "PinuJeff", "PinuMont", "PinuPond", "PinuWash", "PseuMenz", "SequGiga", "TorrCali", 
             "TsugMert", "AcerMacr", "AlnuRhom", "AlnuRubr", "ArbuMenz", "LithDens", "PopuTrem",
             "QuerChry", "QuerKell", "QuerDoug", "QuerWisl", "UmbeCali", "Shrub")


species_ref <- read.csv("D:/Data/fia/FIADB_REFERENCE/REF_SPECIES.csv")
species_ref$species_code <- paste0(
  substr(species_ref$GENUS, 1, 4),
  substr(toupper(species_ref$SPECIES), 1, 1), 
  substr(species_ref$SPECIES, 2, 4)
)

species_ref <- read.csv("D:/Data/fia/FIADB_REFERENCE/REF_PLANT_DICTIONARY.csv") %>%
  dplyr::filter(SYMBOL %in% species_ref$SPECIES_SYMBOL) %>%
  mutate(shrub = grepl("Shrub", .$GROWTH_HABIT, ignore.case = TRUE)) %>%
  dplyr::select(c(SYMBOL, shrub)) %>%
  right_join(species_ref, by = c("SYMBOL" = "SPECIES_SYMBOL")) %>%
  mutate(SPCD = as.character(SPCD))


fia_trees <- left_join(fia_trees, 
                       species_ref[, c("SPCD", "species_code", "shrub", "SFTWD_HRDWD")],
                       by = "SPCD") %>%
  dplyr::filter(species_code %in% species_class) %>%
  dplyr::mutate(species_code = ifelse(shrub, "Shrub", species_code)) #replace shrub species with "Shrub" functional type


# sum biomass for each stand
comm_matrix <- fia_trees %>%
  dplyr::group_by(PLT_CN, species_code) %>%
  summarise(biomass = sum(DRYBIO_TOTAL, na.rm = TRUE)) %>%
  tidyr::pivot_wider(id_cols = PLT_CN, names_from = species_code, values_from = biomass, values_fill = 0)
comm_matrix$total_biomass <- rowSums(comm_matrix[, -1])
#TODO fix SMC and SCN

#fix a problem with species being absent
comm_matrix$PinuMono <- 0
comm_matrix$CWHR_type <- CWHR_type_from_comm_matrix(comm_matrix)

#add CWHR types to fia_plot dataframe
fia_plot <- left_join(fia_plot, comm_matrix[, c("PLT_CN", "CWHR_type")], by = c("CN" = "PLT_CN"))



# Get tree ages ----------------------------------------------------------------

# use regressions to estimate ages and bin to LANDIS ages
regressions_age <- readRDS("linear_models_age_from_diam.RDS")

fia_trees$age_est <- NA

for(i in 1:nrow(fia_trees)){
  if(fia_trees$species_code[i] %in% regressions_age$SpeciesName){
    fia_trees$age_est[i] <- exp(predict(regressions_age$model[match(fia_trees$species_code[i], 
                                                                    regressions_age[[1]])][[1]],
                                        newdata = data.frame(DIA = fia_trees$DIA[i])))
  }
}

fia_trees <- fia_trees %>%
  mutate(TOTAGE = ifelse(!is.na(TOTAGE), TOTAGE, age_est))%>%
  filter(STATUSCD == 1)


################################################################################
# Aggregate FIA data

## MAke plot-level summary variables

tree_summary <- fia_trees %>% 
  dplyr::group_by(PLT_CN) %>%
  # dplyr::filter(n() > 10) %>%
  # dplyr::filter(STATUSCD == 1) %>%
  dplyr::summarise(total_trees = n(),
                   trees_with_age = sum(!is.na(TOTAGE)),
                   plot_bapa = sum(0.005454*(DIA^2)*TPA_UNADJ), 
                   plot_tpa = sum(TPA_UNADJ), 
                   plot_qmd = sqrt((plot_bapa/plot_tpa)/0.005454),
                   plot_sdi = plot_tpa * ((plot_qmd/10)^(-1.605)),
                   sum_sdi = sum(TPA_UNADJ * ((DIA/10)^(-1.605))),
                   biomass = sum(DRYBIO_TOTAL * TPA_UNADJ) * 2.471, #convert to megagrams per ha
                   mean_age = mean(TOTAGE, na.rm = TRUE),
                   high_age = quantile(TOTAGE, 0.9, na.rm = TRUE),
                   .groups = "keep")

# aggregate trees to plot. We need the weighted average age and the total plot
# biomass in order to regress CC ~ age*CWHR*biomass

fia_plot2 <- fia_plot %>%
  filter(cc > 0) %>%
  mutate(cc = cc/100) %>%
  left_join(tree_summary, by = c("CN" = "PLT_CN"))

dia_breaks <- c(0, 1, 6, 11, 24, 1000)
fia_plot2$seral_stage <- cut(fia_plot2$plot_qmd, dia_breaks, labels = FALSE)


#####################################################


under <- read.csv("D:/Data/fia/rFIA_downloads/CA_P2VEG_SUBP_STRUCTURE.csv") %>%
  filter(PLT_CN %in% fia_plot2$CN,
         GROWTH_HABIT_CD == "SH",
         LAYER == 5) %>%
  group_by(PLT_CN, SUBP, CONDID) %>%
  summarise(cover = sum(COVER_PCT)) %>%
  summarise(understory_cover = mean(cover)) %>%
  left_join(fia_plot2, by = c("PLT_CN" = "CN"))


under$understory_cover <- under$understory_cover/100
under$understory_cover <- ifelse(under$understory_cover < 0.001, 0.001, under$understory_cover)
under$understory_cover <- ifelse(under$understory_cover > 0.999, 0.999, under$understory_cover)
under$biomass <- under$biomass * 100
under$biomass <- ifelse(is.na(under$biomass), 0, under$biomass)
under$biomass <- ifelse(under$biomass==0, 1, under$biomass)

plot(under$understory_cover ~ under$biomass)
plot(under$understory_cover ~ under$cc)
plot(under$understory_cover ~ under$seral_stage)
boxplot(under$understory_cover ~ under$CWHR_type)

model <- betareg::betareg(understory_cover ~ biomass + mean_age + cc*CWHR_type, data = under)
summary(model)
plot(effects::allEffects(model))

saveRDS(model, "beta_model_shrub_cover.RDS")

hist(under$understory_cover)

table(under$understory_cover > 0.5, under$cc < 0.1)



table(under$understory_cover > 0.5, under$cc < 0.1)
