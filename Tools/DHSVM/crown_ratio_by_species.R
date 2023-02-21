#get average crown ratio by species
library("tidyverse")

fia_trees <- read.csv("fia_trees.csv")

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
  mutate(SPCD = as.integer(SPCD))

species_class <- c("AbieConc", "AbieMagn", "CaloDecu", "JuniOcci", "PinuAlbi", "PinuBalf", "PinuCont",
                   "PinuLamb", "PinuJeff", "PinuMont", "PinuPond", "PinuWash", "PseuMenz", "SequGiga", "TorrCali", 
                   "TsugMert", "AcerMacr", "AlnuRhom", "AlnuRubr", "ArbuMenz", "LithDens", "PopuTrem",
                   "QuerChry", "QuerKell", "QuerDoug", "QuerWisl", "UmbeCali", "Shrub")

fia_trees <- left_join(fia_trees, 
                       species_ref[, c("SPCD", "species_code", "shrub", "SFTWD_HRDWD")],
                       by = "SPCD") %>%
  dplyr::filter(species_code %in% species_class) %>%
  dplyr::mutate(species_code = ifelse(shrub, "Shrub", species_code)) #replace shrub species with "Shrub" functional type

agg <- fia_trees %>% group_by(species_code) %>%
  summarise(average_crown_ratio = mean(CR, na.rm = TRUE))

agg
