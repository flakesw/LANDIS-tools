
################################################################################
# Estimate height:age regressions

#these are needed to convert LANDIS cohorts to FIA world

library("tidyverse")

#species used for classification
species_class <- c("AbieConc", "AbieMagn", "CaloDecu", "JuniOcci", "PinuAlbi", "PinuBalf", "PinuCont",
                   "PinuLamb", "PinuJeff", "PinuMont", "PinuPond", "PinuWash", "PseuMenz", "SequGiga", "TorrCali", 
                   "TsugMert", "AcerMacr", "AlnuRhom", "AlnuRubr", "ArbuMenz", "LithDens", "PopuTrem",
                   "QuerChry", "QuerKell", "QuerDoug", "QuerWisl", "UmbeCali", "Shrub")


species_ref <- read.csv("D:/Data/fia/FIADB_REFERENCE/REF_SPECIES.csv")
species_ref$SpeciesName <- paste0(
  substr(species_ref$GENUS, 1, 4),
  substr(toupper(species_ref$SPECIES), 1, 1), 
  substr(species_ref$SPECIES, 2, 4)
)

species_ref <- read.csv("D:/Data/fia/FIADB_REFERENCE/REF_PLANT_DICTIONARY.csv") %>%
  dplyr::filter(SYMBOL %in% species_ref$SPECIES_SYMBOL) %>%
  mutate(shrub = grepl("Shrub", .$GROWTH_HABIT, ignore.case = TRUE)) %>%
  dplyr::select(c(SYMBOL, shrub)) %>%
  right_join(species_ref, by = c("SYMBOL" = "SPECIES_SYMBOL"))


#import all FIA data for CA
fia_trees_ca <- read.csv("fia_trees.csv") %>%
  left_join(species_ref[, c("SPCD", "SpeciesName", "shrub", "SFTWD_HRDWD")],
            by = "SPCD") %>%
  dplyr::filter(SpeciesName %in% species_class &
                  !is.na(TOTAGE)) %>%
  dplyr::mutate(SpeciesName = ifelse(shrub, "Shrub", SpeciesName)) #replace shrub species with "Shrub" functional type
  # mutate(TOTAGE = BHAGE + 10) #already done in preprocess

#look at the diameter:age relationship for each species
#it's nonlinear, but should be more or less linear on a log-log scale
for(i in 1:length(unique(fia_trees_ca$SpeciesName))){
  fia_sub <- fia_trees_ca[fia_trees_ca$SpeciesName == unique(fia_trees_ca$SpeciesName)[i],]
  
  if(nrow(fia_sub[!is.na(fia_sub$TOTAGE), ]) > 5){
    plot(I((HT)+rnorm(nrow(fia_sub),0,0.1)) ~ I((TOTAGE)+rnorm(nrow(fia_sub),0,0.11)), 
         data = fia_sub,
         main = unique(fia_trees_ca$SpeciesName)[i],
         log = "xy")
  }
}

#Charles just uses log-log
library(mgcv)
# cub_regressions <- fia_trees_ca %>%
#   dplyr::filter(!is.na(HT) & !is.na(TOTAGE) & !is.na(SPCD)) %>%
#   filter(SpeciesName %in% species_class) %>%
#   dplyr::group_by(SpeciesName) %>%
#   filter(n() > 40) %>%
#   filter(TOTAGE < 500) %>%
#   dplyr::do(model = lm(log(HT) ~ log(TOTAGE), data = .))

# cub_regressions <- fia_trees_ca %>%
#   dplyr::filter(!is.na(HT) & !is.na(TOTAGE) & !is.na(SPCD)) %>%
#   filter(SpeciesName %in% species_class) %>%
#   dplyr::group_by(SpeciesName) %>%
#   filter(n() > 10) %>%
#   dplyr::do(model = gam(log(HT) ~ log(TOTAGE), data = .))

library(earth)
cub_regressions <- fia_trees_ca %>%
  dplyr::filter(!is.na(HT) & !is.na(TOTAGE) & !is.na(SPCD)) %>%
  filter(SpeciesName %in% species_class) %>%
  dplyr::group_by(SpeciesName) %>%
  filter(n() > 40) %>%
  dplyr::do(model = earth(log(HT) ~ log(TOTAGE), data = .))


for(i in 1:12){
  fia_sub <- fia_trees_ca[fia_trees_ca$SpeciesName == cub_regressions$SpeciesName[i], ]
  plot(HT ~ TOTAGE, data = fia_sub,
       main = cub_regressions$SpeciesName[i],
       xlab = "Age (years)",
       ylab = "Height (feet)")
  
  newdata = data.frame(TOTAGE = seq(0, max(fia_sub$TOTAGE, na.rm = TRUE), length.out = 1000))
  preds <- exp(predict(cub_regressions$model[i][[1]], newdata = newdata)) 
            
                 # var(residuals(cub_regressions$model[i][[1]])))
  lines(preds ~ newdata$TOTAGE, lwd = 2, color = "darkgray")

}


no_sp_regression <- earth(log(HT) ~ log(TOTAGE), 
                          data = fia_trees_ca[!is.na(fia_trees_ca$TOTAGE) & 
                                                !is.na(fia_trees_ca$DIA) &
                                                !(fia_trees_ca$SpeciesName %in% cub_regressions$SpeciesName), ])

write_rds(cub_regressions, "linear_models_ht_from_age.RDS")
write_rds(no_sp_regression, "linear_models_ht_from_age_no_sp.RDS")
