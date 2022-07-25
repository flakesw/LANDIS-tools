
################################################################################
# Estimate diameter:age regressions

#these are needed to convert LANDIS cohorts to FIA world
#TODO revise this -- which trees do we want? Can we include CWHR as a predictor?

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
    plot(I((DIA)+rnorm(nrow(fia_sub),0,0.1)) ~ I((TOTAGE)+rnorm(nrow(fia_sub),0,0.11)), 
         data = fia_sub,
         main = unique(fia_trees_ca$SpeciesName)[i],
         log = "xy")
  }
}

#Charles just uses log-log
library(mgcv)
# cub_regressions <- fia_trees_ca %>% 
#   dplyr::filter(!is.na(DIA) & !is.na(TOTAGE) & !is.na(SPCD)) %>%
#   filter(SpeciesName %in% species_class) %>%
#   dplyr::group_by(SpeciesName) %>%
#   filter(n() > 10) %>%
#   filter(TOTAGE < 500) %>%
#   dplyr::do(model = lm(log(DIA) ~ log(TOTAGE), data = .))

# cub_regressions <- fia_trees_ca %>%
#   dplyr::filter(!is.na(DIA) & !is.na(TOTAGE) & !is.na(SPCD)) %>%
#   filter(SpeciesName %in% species_class) %>%
#   dplyr::group_by(SpeciesName) %>%
#   filter(n() > 10) %>%
#   dplyr::do(model = gam(log(DIA) ~ log(TOTAGE), data = .))

library(earth)
cub_regressions <- fia_trees_ca %>%
  dplyr::filter(!is.na(DIA) & !is.na(TOTAGE) & !is.na(SPCD)) %>%
  filter(SpeciesName %in% species_class) %>%
  dplyr::group_by(SpeciesName) %>%
  filter(n() > 40) %>%
  dplyr::do(model = earth(log(DIA) ~ log(TOTAGE), data = .))


for(i in 1:12){
  fia_sub <- fia_trees_ca[fia_trees_ca$SpeciesName == cub_regressions$SpeciesName[i], ]
  plot(DIA ~ TOTAGE, data = fia_sub,
       main = cub_regressions$SpeciesName[i])
  
  newdata = data.frame(TOTAGE = seq(0, max(fia_sub$TOTAGE, na.rm = TRUE), length.out = 1000))
  preds <- exp(predict(cub_regressions$model[i][[1]], newdata = newdata)) 
            
                 # var(residuals(cub_regressions$model[i][[1]])))
  lines(preds ~ newdata$TOTAGE)
  abline(h = 24)
}



# 
# nls_regressions <- fia_trees_ca %>%
#   dplyr::filter(!is.na(DIA) & !is.na(TOTAGE) & !is.na(SPCD)) %>%
#   filter(SpeciesName %in% species_class) %>%
#   dplyr::group_by(SpeciesName) %>%
#   filter(n() > 30) %>%
#   # dplyr::do(model = lm(log(DIA) ~ poly(log(TOTAGE), 1), data = .))
#   dplyr::do(model = nls(log(DIA) ~ SSlogis(log(TOTAGE), Asym, xmid, scal), data = .))
# 
# for(i in 1:22){
#   fia_sub <- fia_trees_ca[fia_trees_ca$SpeciesName == nls_regressions$SpeciesName[i], ]
#   plot(DIA ~ TOTAGE, data = fia_sub,
#        main = nls_regressions$SpeciesName[i])
# 
#   newdata = data.frame(TOTAGE = seq(0, max(fia_sub$TOTAGE, na.rm = TRUE), length.out = 1000))
#   preds <- exp(predict(nls_regressions$model[i][[1]], newdata = newdata) +
#                  var(residuals(nls_regressions$model[i][[1]])))
#   lines(preds ~ newdata$TOTAGE)
# }
# 
# 
# mixed_model <- lme4::lmer(log(DIA) ~ poly(log(TOTAGE), 2) + (1|SpeciesName), data = fia_trees_ca[!is.na(fia_trees_ca$TOTAGE),])
# 
# for(i in 1:length(unique(mixed_model@frame$SpeciesName))){
#   fia_sub <- fia_trees_ca[fia_trees_ca$SpeciesName == unique(mixed_model@frame$SpeciesName)[i], ]
#   plot((DIA) ~ TOTAGE, data = fia_sub,
#        main = unique(mixed_model@frame$SpeciesName)[i])
# 
#   newdata = data.frame(TOTAGE = seq(0, max(fia_sub$TOTAGE, na.rm = TRUE), length.out = 1000),
#                        SpeciesName = unique(mixed_model@frame$SpeciesName)[i])
#   preds <- exp(predict(mixed_model, newdata = newdata) + var(residuals(mixed_model))/2)
#   lines(preds ~ newdata$TOTAGE)
#   abline(h = 0)
#   abline(v = 0)
# }
# 
# 
# newdata <- expand.grid(TOTAGE = seq(5, 500, by = 5),
#                        SpeciesName = unique(mixed_model@frame$SpeciesName),
#                        preds_lm = NA,
#                        preds_nls = NA,
#                        preds_mixed = NA)
# for(i in 1:nrow(newdata)){
#   newdata$preds_lm[i] <- exp(predict(cub_regressions$model[match(newdata$SpeciesName[i], cub_regressions[[1]])][[1]],
#                                      newdata = data.frame(TOTAGE = newdata$TOTAGE[i]))) 
#                              # +
#                              #   var(residuals(cub_regressions$model[match(newdata$SpeciesName[i], cub_regressions[[1]])][[1]], na.rm = TRUE))/2)
#   newdata$preds_nls[i] <- exp(predict(nls_regressions$model[match(newdata$SpeciesName[i], nls_regressions[[1]])][[1]],
#                                       newdata = data.frame(TOTAGE = newdata$TOTAGE[i])) +
#                                 var(residuals(nls_regressions$model[match(newdata$SpeciesName[i], nls_regressions[[1]])][[1]], na.rm = TRUE))/2)
# 
# }
# 
# 
# newdata$preds_mixed <- exp(predict(mixed_model, newdata = newdata) +
#                              var(residuals(mixed_model))/2)
# 
# charles_data <- read.csv("pred_values_all_spp_catrees.csv") %>%
#   rename(TOTAGE = TOTAGE2,
#          preds1 = exp.pmd1.) %>%
#   left_join(species_ref[, c("COMMON_NAME", "SpeciesName")], by = c("curr_name" = "COMMON_NAME")) %>%
#   left_join(newdata, by = c("SpeciesName", "TOTAGE"))
# 
# plot(charles_data$preds_lm ~ charles_data$TOTAGE, xlim = c(0,200), ylim = c(0, 50))
# plot(charles_data$preds_nls ~ charles_data$TOTAGE, xlim = c(0,200), ylim = c(0, 50))
# plot(charles_data$preds_mixed ~ charles_data$TOTAGE, xlim = c(0,200), ylim = c(0, 50))
# abline(h = 24)
# 
# plot(charles_data$preds_lm ~ charles_data$preds1, xlim = c(0,100), ylim = c(0,100),
#      xlab = "Charles' predicted diameter",
#      ylab = "Sam's predicted diameter")
# abline(0,1)
# 
# 
# plot(charles_data[charles_data$TOTAGE < 250, ]$preds_lm ~ charles_data[charles_data$TOTAGE < 250, ]$preds1, 
#      xlim = c(0,50), ylim = c(0,50),
#      xlab = "Charles' predicted diameter (in)",
#      ylab = "Sam's predicted diameter (in)")
# abline(0,1)
# 
# plot(charles_data[charles_data$TOTAGE < 250, ]$preds_nls ~ charles_data[charles_data$TOTAGE < 250, ]$preds1, 
#      xlim = c(0,50), ylim = c(0,50),
#      xlab = "Charles' predicted diameter (in)",
#      ylab = "Sam's predicted diameter (in)")
# abline(0,1)
# 
# plot(charles_data[charles_data$TOTAGE < 250, ]$preds_mixed ~ charles_data[charles_data$TOTAGE < 250, ]$preds1, 
#      xlim = c(0,50), ylim = c(0,50),
#      xlab = "Charles' predicted diameter (in)",
#      ylab = "Sam's predicted diameter (in)")
# abline(0,1)
# 
# write.csv(charles_data, "compare_charles_sam_diameters.csv")
# 
# charles_data$diff <- charles_data$preds1 - charles_data$preds_lm
# mean(charles_data$diff, na.rm = TRUE)

# no_sp_regression <- gam(DIA ~ log(TOTAGE), data = fia_trees_ca)
no_sp_regression <- earth(log(DIA) ~ log(TOTAGE), 
                          data = fia_trees_ca[!is.na(fia_trees_ca$TOTAGE) & 
                                                !is.na(fia_trees_ca$DIA) &
                                                !(fia_trees_ca$SpeciesName %in% cub_regressions$SpeciesName), ])

write_rds(cub_regressions, "linear_models_diam_from_age.RDS")
write_rds(no_sp_regression, "linear_models_diam_from_age_no_sp.RDS")
# write_rds(mixed_model, "mixed_model_diam_from_age.RDS")

################################################################################
# Estimate how much biomass per tree to reconstruct TPA
#
#This didn't work :(
#
# plot(DRYBIO_TOTAL ~ TOTAGE, data = fia_trees_ca)
# boxplot(DRYBIO_TOTAL ~ SpeciesName, data = fia_trees_ca)
# 
# test_model <- lm(log(DRYBIO_TOTAL) ~TOTAGE*SpeciesName, data = fia_trees_ca)
# summary(test_model)
# 
# biomass_regressions <- fia_trees_ca %>%
#   dplyr::filter(!is.na(DIA) & !is.na(TOTAGE) & !is.na(SPCD)) %>%
#   filter(SpeciesName %in% species_class) %>%
#   dplyr::group_by(SpeciesName) %>%
#   filter(n() > 20) %>%
#   dplyr::do(model = earth(log(DRYBIO_TOTAL) ~ log(TOTAGE), data = .))
# 
# 
# for(i in 1:15){
#   fia_sub <- fia_trees_ca[fia_trees_ca$SpeciesName == biomass_regressions$SpeciesName[i], ]
#   plot(DRYBIO_TOTAL ~ TOTAGE, data = fia_sub,
#        main = biomass_regressions$SpeciesName[i])
#   
#   newdata = data.frame(TOTAGE = seq(0, max(fia_sub$TOTAGE, na.rm = TRUE), length.out = 1000))
#   preds <- exp(predict(biomass_regressions$model[i][[1]], newdata = newdata)) 
#   
#   # var(residuals(cub_regressions$model[i][[1]])))
#   lines(preds ~ newdata$TOTAGE)
# }
# 
# no_sp_regression_biomass <- earth(log(DRYBIO_TOTAL) ~ log(TOTAGE), 
#                           data = fia_trees_ca[!is.na(fia_trees_ca$TOTAGE) & 
#                                                 !is.na(fia_trees_ca$DRYBIO_TOTAL) &
#                                                 !(fia_trees_ca$SpeciesName %in% cub_regressions$SpeciesName), ])
# 
# write_rds(biomass_regressions, "linear_models_biomass_from_age.RDS")
# write_rds(no_sp_regression_biomass, "linear_models_biomass_from_age_no_sp.RDS")

################################################################################
# Estimate age:diameter regressions

#these are needed to convert FIA plots into LANDIS biomass-age cohorts

#import all FIA data for CA
fia_trees_ca <- read.csv("D:/Data/fia/rFIA_downloads/CA_TREE.csv") %>%
  left_join(species_ref[, c("SPCD", "SpeciesName", "shrub", "SFTWD_HRDWD")],
            by = "SPCD") %>%
  dplyr::filter(SpeciesName %in% species_class) %>%
  dplyr::mutate(SpeciesName = ifelse(shrub, "Shrub", SpeciesName)) %>%#replace shrub species with "Shrub" functional type
  mutate(TOTAGE = BHAGE+10)

#look at the diameter:age relationship for each species
#it's nonlinear, but should be more or less linear on a log-log scale
for(i in 1:length(unique(fia_trees_ca$SpeciesName))){
  fia_sub <- fia_trees_ca[fia_trees_ca$SpeciesName == unique(fia_trees_ca$SpeciesName)[i],]
  
  if(nrow(fia_sub[!is.na(fia_sub$TOTAGE), ]) > 5){
    plot(I(log(TOTAGE)+rnorm(nrow(fia_sub),0,0.1)) ~ I(log(DIA)+rnorm(nrow(fia_sub),0,0.11)), 
         data = fia_sub,
         main = unique(fia_trees_ca$SpeciesName)[i])
  }
}

cub_regressions_age <- fia_trees_ca %>% 
  dplyr::filter(!is.na(DIA) & !is.na(TOTAGE) & !is.na(SPCD)) %>%
  filter(SpeciesName %in% species_class) %>%
  dplyr::group_by(SpeciesName) %>%
  filter(n() > 10) %>%
  dplyr::do(model = lm(log(TOTAGE) ~ poly(log(DIA), 3), data = .))

for(i in 1:18){
  fia_sub <- fia_trees[fia_trees$SpeciesName == cub_regressions_age$SpeciesName[i], ]
  plot(TOTAGE ~ DIA, data = fia_sub,
       main = cub_regressions_age$SpeciesName[i],
       ylim = c(0, 400))
  
  newdata = data.frame(DIA = seq(1, max(fia_sub$DIA, na.rm = TRUE), length.out = 1000))
  preds <- exp(predict(cub_regressions_age$model[i][[1]], newdata = newdata))
  lines(preds ~ newdata$DIA)
  abline(v = 12)
}

nls_regressions_age <- fia_trees_ca %>% 
  dplyr::filter(!is.na(DIA) & !is.na(TOTAGE) & !is.na(SPCD)) %>%
  filter(SpeciesName %in% species_class) %>%
  dplyr::group_by(SpeciesName) %>%
  filter(n() > 10) %>%
  # dplyr::do(model = lm(log(DIA) ~ poly(log(TOTAGE), 1), data = .))
  dplyr::do(model = nls(log(TOTAGE) ~ SSlogis(DIA, Asym, xmid, scal), data = .))

for(i in 1:18){
  fia_sub <- fia_trees[fia_trees$SpeciesName == agedia_regressions$SpeciesName[i], ]
  plot(TOTAGE ~ DIA, data = fia_sub,
       main = agedia_regressions$SpeciesName[i],
       ylim = c(0, 400))
  
  newdata = data.frame(DIA = seq(0, max(fia_sub$DIA, na.rm = TRUE), length.out = 1000))
  preds <- exp(predict(nls_regressions_age$model[i][[1]], newdata = newdata)) 
  lines(preds ~ newdata$DIA)
  abline(v = 24)
}


mixed_model_age <- lme4::lmer(log(TOTAGE) ~ poly(log(DIA), 3) + (1|SpeciesName), data = fia_trees_ca[!is.na(fia_trees_ca$TOTAGE),])

for(i in 1:length(unique(mixed_model_age@frame$SpeciesName))){
  fia_sub <- fia_trees_ca[fia_trees_ca$SpeciesName == unique(mixed_model_age@frame$SpeciesName)[i], ]
  plot((TOTAGE) ~ DIA, data = fia_sub,
       main = unique(mixed_model_age@frame$SpeciesName)[i],
       ylim = c(0, 400))
  
  newdata = data.frame(DIA = seq(0, max(fia_sub$TOTAGE, na.rm = TRUE), length.out = 1000),
                       SpeciesName = unique(mixed_model_age@frame$SpeciesName)[i])
  preds <- exp(predict(mixed_model_age, newdata = newdata) + var(residuals(mixed_model_age))/2)
  lines(preds ~ newdata$DIA)
  abline(h = 0)
  abline(v = 0)
  abline(v = 24)
}


#compare predictions
newdata <- expand.grid(DIA = seq(0, 200, by = 5),
                       SpeciesName = unique(cub_regressions_age$SpeciesName),
                       preds = NA)
for(i in 1:nrow(newdata)){
  newdata$preds[i] <- exp(predict(cub_regressions_age$model[match(newdata$SpeciesName[i], 
                                                                  cub_regressions_age[[1]])][[1]],
                                  newdata = data.frame(DIA = newdata$DIA[i])))
}


newdata$preds_mixed <- exp(predict(mixed_model_age, newdata = newdata) + var(residuals(mixed_model_age))/2)

charles_data <- read.csv("pred_values_all_spp_catrees.csv") %>%
  rename(TOTAGE = TOTAGE2,
         DIA = exp.pmd1.) %>%
  left_join(species_ref[, c("COMMON_NAME", "SpeciesName")], by = c("curr_name" = "COMMON_NAME"))# %>%
# left_join(newdata, by = c("SpeciesName", "DIA"))

plot(charles_data$TOTAGE ~ charles_data$DIA, xlim = c(0,50), ylim = c(0, 400))
plot(newdata$preds ~ newdata$DIA, xlim = c(0,50), ylim = c(0, 400))
plot(newdata$preds_mixed ~ newdata$DIA, xlim = c(0,50), ylim = c(0, 400))
abline(h = 100)
abline(v = 24)


write_rds(cub_regressions_age, "linear_models_age_from_diam.RDS")