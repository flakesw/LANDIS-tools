# Estimating CC from biomass

library("tidyverse")

################################################################################
# First step: classify FIA stands into CWHR community types

#species used for classification
species_class <- c("AbieConc", "AbieMagn", "CaloDecu", "JuniOcci", "PinuAlbi", "PinuBalf", "PinuCont",
             "PinuJeff", "PinuMont", "PinuPond", "PinuWash", "PseuMenz", "SequGiga", "TorrCali", 
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
  right_join(species_ref, by = c("SYMBOL" = "SPECIES_SYMBOL"))


# FIA plot data which has been trimmed to the study area and somewhat pre-processed
# using data from the COND table
# See "preprocess_fia_data.R"
fia_plot <- read.csv("fia_plot.csv")
fia_trees <- read.csv("fia_trees.csv")

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

# classify

comm_matrix$SMC <- rowSums(comm_matrix[, c("AbieConc", "CaloDecu", "PinuCont", 
                                           "PinuJeff", "PinuPond", #"PinuLamb", removed PinuLamb because absent from FIA 
                                           "PinuWash", "PseuMenz", "TorrCali",
                                           "PinuBalf")])
comm_matrix$SCN <- rowSums(comm_matrix[, c("AbieMagn", "PinuAlbi", "PinuMont", "TsugMert")])
comm_matrix$coni <-  rowSums(comm_matrix[, c("AbieConc", "AbieMagn", "CaloDecu", "JuniOcci",
                                             "PinuAlbi", "PinuCont", "PinuJeff", #PinuLamb,
                                             "PinuMont", "PinuPond", "PinuWash", "PseuMenz",
                                             "TsugMert", "PinuBalf", "SequGiga")])
comm_matrix$hard <- rowSums(comm_matrix[, c("AcerMacr", "AlnuRhom", "AlnuRubr", "ArbuMenz",
                                            "LithDens", "PopuTrem", #"QuerChry", "QuerKell","QuerWisl", "UmbeCali"
                                            "QuerDoug")])
comm_matrix <- comm_matrix %>%
  mutate(cover.coni = coni > hard & coni > Shrub,
         cover.hard = hard > coni & (coni/total_biomass) < 0.25,
         cover.mixed = hard > coni & (coni/total_biomass) > 0.25,
         cover.shrub = Shrub > hard & Shrub > coni & !cover.mixed & !cover.hard & !cover.coni)

#this is the dumbest thing I've ever written
comm_matrix$CWHR_type <- 
  ifelse(comm_matrix$cover.mixed, "mhc",
  ifelse(comm_matrix$cover.hard & 
           comm_matrix$PopuTrem / comm_matrix$hard > 0.5, "asp",
  ifelse(comm_matrix$cover.hard & 
           ((comm_matrix$ArbuMenz + comm_matrix$LithDens + comm_matrix$QuerDoug) / comm_matrix$hard) > 0.5,
            "mhw",
  ifelse(comm_matrix$cover.hard & 
           ((comm_matrix$AcerMacr + comm_matrix$AlnuRhom + comm_matrix$AlnuRubr)/comm_matrix$hard) > 0.5,
         "mri",
  ifelse(comm_matrix$cover.coni & comm_matrix$AbieConc/comm_matrix$coni > 0.5, "wfr",
  ifelse(comm_matrix$cover.coni & comm_matrix$AbieMagn/comm_matrix$coni > 0.5, "rfr",
  ifelse(comm_matrix$cover.coni & comm_matrix$PinuJeff/comm_matrix$coni > 0.5, "jpn",
  ifelse(comm_matrix$cover.coni & comm_matrix$PseuMenz/comm_matrix$coni > 0.5, "dfr",
  ifelse(comm_matrix$cover.coni & 
           (comm_matrix$PinuPond + comm_matrix$PinuWash)/comm_matrix$coni > 0.5, "ppn",
  ifelse(comm_matrix$cover.coni & comm_matrix$PinuCont / comm_matrix$coni > 0.5, "lpn",
  ifelse(comm_matrix$cover.coni & comm_matrix$JuniOcci / comm_matrix$coni > 0.5, "juo",
  ifelse(comm_matrix$cover.coni & comm_matrix$SequGiga / comm_matrix$coni > 0.5, "seg",
  ifelse(comm_matrix$Shrub, "shrub", 
  ifelse(comm_matrix$cover.coni & 
           comm_matrix$AbieConc / comm_matrix$coni < 0.5 & 
           comm_matrix$AbieMagn / comm_matrix$coni < 0.5 &
           comm_matrix$CaloDecu / comm_matrix$coni < 0.5 &
           comm_matrix$JuniOcci / comm_matrix$coni < 0.5 &
           comm_matrix$PinuAlbi / comm_matrix$coni < 0.5 &
           comm_matrix$PinuCont / comm_matrix$coni < 0.5 &
           comm_matrix$PinuJeff / comm_matrix$coni < 0.5 &
           (comm_matrix$PinuPond + comm_matrix$PinuWash) / comm_matrix$coni < 0.5 &
           comm_matrix$PseuMenz / comm_matrix$coni < 0.5 &
           comm_matrix$TsugMert / comm_matrix$coni < 0.5 &
           comm_matrix$PinuBalf / comm_matrix$coni < 0.5 &
           comm_matrix$SequGiga / comm_matrix$coni < 0.5 &
           comm_matrix$SMC > comm_matrix$SCN,
           "smc",
  ifelse(comm_matrix$cover.coni & 
           comm_matrix$AbieConc / comm_matrix$coni < 0.5 & 
           comm_matrix$AbieMagn / comm_matrix$coni < 0.5 &
           comm_matrix$CaloDecu / comm_matrix$coni < 0.5 &
           comm_matrix$JuniOcci / comm_matrix$coni < 0.5 &
           comm_matrix$PinuAlbi / comm_matrix$coni < 0.5 &
           comm_matrix$PinuCont / comm_matrix$coni < 0.5 &
           comm_matrix$PinuJeff / comm_matrix$coni < 0.5 &
           (comm_matrix$PinuPond + comm_matrix$PinuWash) / comm_matrix$coni < 0.5 &
           comm_matrix$PseuMenz / comm_matrix$coni < 0.5 &
           comm_matrix$TsugMert / comm_matrix$coni < 0.5 & 
           comm_matrix$PinuBalf / comm_matrix$coni < 0.5 &
           comm_matrix$SequGiga / comm_matrix$coni < 0.5 &
           comm_matrix$SMC < comm_matrix$SCN,
         "scn",
  ifelse(comm_matrix$cover.coni & 
           (comm_matrix$PinuJeff + comm_matrix$PinuPond + comm_matrix$PinuWash)/comm_matrix$coni > 0.5,
         "yps",
  ifelse(comm_matrix$cover.coni & 
           (comm_matrix$PinuMont + comm_matrix$PinuAlbi) / comm_matrix$coni > 0.5,
         "wps",
         NA)))))))))))))))))
  
table(comm_matrix$CWHR_type)



################################################################################
# Estimate diameter:age regressions

#these are needed to convert LANDIS cohorts to FIA world

#import all FIA data for CA
fia_trees_ca <- read.csv("D:/Data/fia/rFIA_downloads/CA_TREE.csv") %>%
  left_join(species_ref[, c("SPCD", "species_code", "shrub", "SFTWD_HRDWD")],
            by = "SPCD") %>%
  dplyr::filter(species_code %in% species_class) %>%
  dplyr::mutate(species_code = ifelse(shrub, "Shrub", species_code)) %>%#replace shrub species with "Shrub" functional type
  mutate(TOTAGE = TOTAGE + 10)

#look at the diameter:age relationship for each species
#it's nonlinear, but should be more or less linear on a log-log scale
for(i in 1:length(unique(fia_trees_ca$species_code))){
  fia_sub <- fia_trees_ca[fia_trees_ca$species_code == unique(fia_trees_ca$species_code)[i],]
  
  if(nrow(fia_sub[!is.na(fia_sub$TOTAGE), ]) > 5){
     plot(I((DIA)+rnorm(nrow(fia_sub),0,0.1)) ~ I((TOTAGE)+rnorm(nrow(fia_sub),0,0.11)), 
          data = fia_sub,
          main = unique(fia_trees_ca$species_code)[i])
  }
}

cub_regressions <- fia_trees_ca %>% 
  dplyr::filter(!is.na(DIA) & !is.na(TOTAGE) & !is.na(SPCD)) %>%
  filter(species_code %in% species_class) %>%
  dplyr::group_by(species_code) %>%
  filter(n() > 10) %>%
  dplyr::do(model = lm(log(DIA) ~ poly(log(TOTAGE), 2), data = .))

for(i in 1:18){
  fia_sub <- fia_trees[fia_trees$species_code == cub_regressions$species_code[i], ]
  plot(DIA ~ TOTAGE, data = fia_sub,
       main = cub_regressions$species_code[i])
  
  newdata = data.frame(TOTAGE = seq(0, max(fia_sub$TOTAGE, na.rm = TRUE), length.out = 1000))
  preds <- exp(predict(cub_regressions$model[i][[1]], newdata = newdata) + 
                 var(residuals(cub_regressions$model[i][[1]])))
  lines(preds ~ newdata$TOTAGE)
}
  
nls_regressions <- fia_trees_ca %>% 
  dplyr::filter(!is.na(DIA) & !is.na(TOTAGE) & !is.na(SPCD)) %>%
  filter(species_code %in% species_class) %>%
  dplyr::group_by(species_code) %>%
  filter(n() > 10) %>%
  # dplyr::do(model = lm(log(DIA) ~ poly(log(TOTAGE), 1), data = .))
  dplyr::do(model = nls(log(DIA) ~ SSlogis(log(TOTAGE), Asym, xmid, scal), data = .))

for(i in 1:18){
  fia_sub <- fia_trees[fia_trees$species_code == agedia_regressions$species_code[i], ]
  plot(DIA ~ TOTAGE, data = fia_sub,
       main = agedia_regressions$species_code[i])
  
  newdata = data.frame(TOTAGE = seq(0, max(fia_sub$TOTAGE, na.rm = TRUE), length.out = 1000))
  preds <- exp(predict(agedia_regressions$model[i][[1]], newdata = newdata) + 
                 var(residuals(agedia_regressions$model[i][[1]]))) 
  lines(preds ~ newdata$TOTAGE)
}


mixed_model <- lme4::lmer(log(DIA) ~ poly(log(TOTAGE), 2) + (1|species_code), data = fia_trees_ca[!is.na(fia_trees_ca$TOTAGE),])

for(i in 1:length(unique(mixed_model@frame$species_code))){
  fia_sub <- fia_trees_ca[fia_trees_ca$species_code == unique(mixed_model@frame$species_code)[i], ]
  plot((DIA) ~ TOTAGE, data = fia_sub,
       main = unique(mixed_model@frame$species_code)[i])
  
  newdata = data.frame(TOTAGE = seq(0, max(fia_sub$TOTAGE, na.rm = TRUE), length.out = 1000),
                       species_code = unique(mixed_model@frame$species_code)[i])
  preds <- exp(predict(mixed_model, newdata = newdata) + var(residuals(mixed_model))/2)
  lines(preds ~ newdata$TOTAGE)
  abline(h = 0)
  abline(v = 0)
}


newdata <- expand.grid(TOTAGE = seq(5, 500, by = 5),
                       species_code = unique(agedia_regressions$species_code),
                       preds_lm = NA,
                       preds_nls = NA,
                       preds_mixed = NA)
for(i in 1:nrow(newdata)){
  newdata$preds_lm[i] <- exp(predict(cub_regressions$model[match(newdata$species_code[i], cub_regressions[[1]])][[1]],
                         newdata = data.frame(TOTAGE = newdata$TOTAGE[i])) +
                           var(residuals(cub_regressions$model[match(newdata$species_code[i], cub_regressions[[1]])][[1]], na.rm = TRUE))/2)
  newdata$preds_nls[i] <- exp(predict(nls_regressions$model[match(newdata$species_code[i], nls_regressions[[1]])][[1]],
                                        newdata = data.frame(TOTAGE = newdata$TOTAGE[i])) +
                                  var(residuals(nls_regressions$model[match(newdata$species_code[i], nls_regressions[[1]])][[1]], na.rm = TRUE))/2)
  
}


newdata$preds_mixed <- exp(predict(mixed_model, newdata = newdata) + 
                             var(residuals(mixed_model))/2)

charles_data <- read.csv("pred_values_all_spp_catrees.csv") %>%
  rename(TOTAGE = TOTAGE2,
         preds1 = exp.pmd1.) %>%
  left_join(species_ref[, c("COMMON_NAME", "species_code")], by = c("curr_name" = "COMMON_NAME")) %>%
  left_join(newdata, by = c("species_code", "TOTAGE"))

plot(charles_data$preds_lm ~ charles_data$TOTAGE, xlim = c(0,200), ylim = c(0, 50))
plot(charles_data$preds_nls ~ charles_data$TOTAGE, xlim = c(0,200), ylim = c(0, 50))
plot(charles_data$preds_mixed ~ charles_data$TOTAGE, xlim = c(0,200), ylim = c(0, 50))
abline(h = 24)

plot(charles_data$preds_nls ~ charles_data$preds1, xlim = c(0,100), ylim = c(0,100),
     xlab = "Charles' predicted diameter",
     ylab = "Sam's predicted diameter")
abline(0,1)


plot(charles_data[charles_data$TOTAGE < 250, ]$preds_lm ~ charles_data[charles_data$TOTAGE < 250, ]$preds1, 
     xlim = c(0,50), ylim = c(0,50),
     xlab = "Charles' predicted diameter (in)",
     ylab = "Sam's predicted diameter (in)")
abline(0,1)

plot(charles_data[charles_data$TOTAGE < 250, ]$preds_nls ~ charles_data[charles_data$TOTAGE < 250, ]$preds1, 
     xlim = c(0,50), ylim = c(0,50),
     xlab = "Charles' predicted diameter (in)",
     ylab = "Sam's predicted diameter (in)")
abline(0,1)

plot(charles_data[charles_data$TOTAGE < 250, ]$preds_mixed ~ charles_data[charles_data$TOTAGE < 250, ]$preds1, 
     xlim = c(0,50), ylim = c(0,50),
     xlab = "Charles' predicted diameter (in)",
     ylab = "Sam's predicted diameter (in)")
abline(0,1)

write.csv(charles_data, "compare_charles_sam_diameters.csv")

charles_data$diff <- charles_data$preds1 - charles_data$preds_mixed
mean(charles_data$diff, na.rm = TRUE)




################################################################################
# Estimate age:diameter regressions

#these are needed to convert FIA plots into LANDIS biomass-age cohorts

#import all FIA data for CA
fia_trees_ca <- read.csv("D:/Data/fia/rFIA_downloads/CA_TREE.csv") %>%
  left_join(species_ref[, c("SPCD", "species_code", "shrub", "SFTWD_HRDWD")],
            by = "SPCD") %>%
  dplyr::filter(species_code %in% species_class) %>%
  dplyr::mutate(species_code = ifelse(shrub, "Shrub", species_code)) %>%#replace shrub species with "Shrub" functional type
  mutate(TOTAGE = TOTAGE + 10)

#look at the diameter:age relationship for each species
#it's nonlinear, but should be more or less linear on a log-log scale
for(i in 1:length(unique(fia_trees_ca$species_code))){
  fia_sub <- fia_trees_ca[fia_trees_ca$species_code == unique(fia_trees_ca$species_code)[i],]
  
  if(nrow(fia_sub[!is.na(fia_sub$TOTAGE), ]) > 5){
    plot(I(log(TOTAGE)+rnorm(nrow(fia_sub),0,0.1)) ~ I(log(DIA)+rnorm(nrow(fia_sub),0,0.11)), 
         data = fia_sub,
         main = unique(fia_trees_ca$species_code)[i])
  }
}

cub_regressions <- fia_trees_ca %>% 
  dplyr::filter(!is.na(DIA) & !is.na(TOTAGE) & !is.na(SPCD)) %>%
  filter(species_code %in% species_class) %>%
  dplyr::group_by(species_code) %>%
  filter(n() > 10) %>%
  dplyr::do(model = lm(log(TOTAGE) ~ poly(log(DIA), 3), data = .))

for(i in 1:18){
  fia_sub <- fia_trees[fia_trees$species_code == cub_regressions$species_code[i], ]
  plot(TOTAGE ~ DIA, data = fia_sub,
       main = cub_regressions$species_code[i],
       ylim = c(0, 400))
  
  newdata = data.frame(DIA = seq(1, max(fia_sub$DIA, na.rm = TRUE), length.out = 1000))
  preds <- exp(predict(cub_regressions$model[i][[1]], newdata = newdata) +
                   var(residuals(cub_regressions$model[i][[1]]))/2)
  lines(preds ~ newdata$DIA)
  abline(v = 12)
}

nls_regressions <- fia_trees_ca %>% 
  dplyr::filter(!is.na(DIA) & !is.na(TOTAGE) & !is.na(SPCD)) %>%
  filter(species_code %in% species_class) %>%
  dplyr::group_by(species_code) %>%
  filter(n() > 10) %>%
  # dplyr::do(model = lm(log(DIA) ~ poly(log(TOTAGE), 1), data = .))
  dplyr::do(model = nls(log(TOTAGE) ~ SSlogis(DIA, Asym, xmid, scal), data = .))

for(i in 1:18){
  fia_sub <- fia_trees[fia_trees$species_code == agedia_regressions$species_code[i], ]
  plot(TOTAGE ~ DIA, data = fia_sub,
       main = agedia_regressions$species_code[i],
       ylim = c(0, 400))
  
  newdata = data.frame(DIA = seq(0, max(fia_sub$DIA, na.rm = TRUE), length.out = 1000))
  preds <- exp(predict(nls_regressions$model[i][[1]], newdata = newdata) + 
                 var(residuals(nls_regressions$model[i][[1]]))) 
  lines(preds ~ newdata$DIA)
  abline(v = 24)
}


mixed_model <- lme4::lmer(log(TOTAGE) ~ poly(log(DIA), 3) + (1|species_code), data = fia_trees_ca[!is.na(fia_trees_ca$TOTAGE),])

for(i in 1:length(unique(mixed_model@frame$species_code))){
  fia_sub <- fia_trees_ca[fia_trees_ca$species_code == unique(mixed_model@frame$species_code)[i], ]
  plot((TOTAGE) ~ DIA, data = fia_sub,
       main = unique(mixed_model@frame$species_code)[i],
       ylim = c(0, 400))
  
  newdata = data.frame(DIA = seq(0, max(fia_sub$TOTAGE, na.rm = TRUE), length.out = 1000),
                       species_code = unique(mixed_model@frame$species_code)[i])
  preds <- exp(predict(mixed_model, newdata = newdata) + var(residuals(mixed_model))/2)
  lines(preds ~ newdata$DIA)
  abline(h = 0)
  abline(v = 0)
  abline(v = 24)
}


newdata <- expand.grid(DIA = seq(0, 200, by = 5),
                       species_code = unique(cub_regressions$species_code),
                       preds = NA)
for(i in 1:nrow(newdata)){
  newdata$preds[i] <- exp(predict(cub_regressions$model[match(newdata$species_code[i], cub_regressions[[1]])][[1]],
                                  newdata = data.frame(DIA = newdata$DIA[i])) +
                            var(residuals(cub_regressions$model[match(newdata$species_code[i], cub_regressions[[1]])][[1]], na.rm = TRUE))/2)
}


newdata$preds_mixed <- exp(predict(mixed_model, newdata = newdata) + var(residuals(mixed_model))/2)

charles_data <- read.csv("pred_values_all_spp_catrees.csv") %>%
  rename(TOTAGE = TOTAGE2,
         DIA = exp.pmd1.) %>%
  left_join(species_ref[, c("COMMON_NAME", "species_code")], by = c("curr_name" = "COMMON_NAME"))# %>%
  # left_join(newdata, by = c("species_code", "DIA"))

plot(charles_data$TOTAGE ~ charles_data$DIA, xlim = c(0,50), ylim = c(0, 400))
plot(newdata$preds ~ newdata$DIA, xlim = c(0,50), ylim = c(0, 400))
plot(newdata$preds_mixed ~ newdata$DIA, xlim = c(0,50), ylim = c(0, 400))
abline(h = 100)
abline(v = 24)


################################################################################
# Aggregate FIA data

#calculate canopy cover from COND table

#convert diameter to age
#calculate stand biomass
#



################################################################################
# Fit models

# cc ~ age*forest type*biomass




################################################################################
# Export data

