# Estimating CC from biomass
# Inputs needed
#   pre-processed FIA data
#   regressions to get age from diameter

# Outputs:
#   Regressions for CC ~ CHWR type + Age + Biomass

library("tidyverse")


logit <- function(x) log(x/(1-x))
invlogit <- function(x) exp(x)/(1+exp(x))


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
#TODO redo this with new code for TCSI

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
  
#add CWHR types to fia_plot dataframe
fia_plot <- left_join(fia_plot, comm_matrix[, c("PLT_CN", "CWHR_type")], by = c("CN" = "PLT_CN"))

################################################################################
# Aggregate FIA data

# use regressions to estimate ages and bin to LANDIS ages
regressions_age <- readRDS("linear_models_age_from_diam.RDS")

fia_trees$age_est <- NA

for(i in 1:nrow(fia_trees)){
  if(fia_trees$species_code[i] %in% regressions_age$species_code){
    fia_trees$age_est[i] <- exp(predict(regressions_age$model[match(fia_trees$species_code[i], 
                                                                    regressions_age[[1]])][[1]],
                                    newdata = data.frame(DIA = fia_trees$DIA[i])))
  }
}

fia_trees <- fia_trees %>%
  mutate(age_all = ifelse(!is.na(TOTAGE), TOTAGE, age_est))

# aggregate trees to plot. We need the weighted average age and the total plot
# biomass in order to regress CC ~ age*CWHR*biomass

tree_summary <- fia_trees %>%
  filter(DIA > 5) %>%
  group_by(PLT_CN) %>%
  summarise(total_biomass = sum(DRYBIO_TOTAL * TPA_UNADJ) / 892, #convert to megagrams per ha,
            mean_age = weighted.mean(age_est, DRYBIO_TOTAL, na.rm = TRUE))

fia_plot2 <- left_join(fia_plot, tree_summary, by = c("CN" = "PLT_CN")) %>%
  filter(fia_plot$cc > 0) %>%
  mutate(cc = cc/100)


#exploratory plots
plot(cc ~ total_biomass, data = fia_plot2)
plot(cc ~ sqrt(total_biomass), data = fia_plot2) #strong relationship
plot(logit(cc) ~ sqrt(total_biomass), data = fia_plot2) #strong relationship
plot(cc ~ mean_age, data = fia_plot2) #little relationship
boxplot(cc ~ CWHR_type, data = fia_plot2)

#remove CWHR types without enough data to fit a model
fia_plot_remove_rare <- fia_plot2 %>%
  group_by(CWHR_type) %>%
  filter(mean_age > 0 & cc > 0 & total_biomass > 0) %>%
  filter(total_biomass < 500) %>%
  filter(n() > 20) 




################################################################################
# Fit models

# cc ~ age*forest type*biomass

#fit models
cc_lm <- lm(logit(cc) ~ sqrt(total_biomass)*mean_age*CWHR_type, data = fia_plot_remove_rare)
summary(cc_lm)

cc_no_cwhr_lm <- lm(logit(cc) ~ sqrt(total_biomass)*mean_age, data = fia_plot2)
summary(cc_no_cwhr_lm)
# plot(effects::allEffects(cc_no_cwhr_lm))

cwhr_types <- unique(fia_plot_remove_rare$CWHR_type)[!is.na(unique(fia_plot_remove_rare$CWHR_type))]

newdata = expand.grid(total_biomass = seq(0, 1000, length.out = 1000),
                      mean_age = c(0, 10, 30, 70, 120, 240),
                      CWHR_type = cwhr_types)

for(i in 1:length(cwhr_types)){
  fia_plot_sub <- fia_plot_remove_rare[fia_plot_remove_rare$CWHR_type == cwhr_types[i], ]
  fia_tree_sub <- fia_trees[fia_trees$PLT_CN %in% fia_plot_sub$CN, ]
  
  plot(fia_plot_sub$cc ~ fia_plot_sub$total_biomass,
       main = cwhr_types[i], 
       xlim = c(0, 1000),
       ylim = c(0, 1))
  abline(h = c(0.4, 0.6))
  for(j in unique(newdata$mean_age)){
    newdat <- newdata[newdata$CWHR_type == cwhr_types[i] &
                                  newdata$mean_age == j, ]
    preds <- invlogit(predict(cc_lm, newdat))
    lines( preds ~ newdat$total_biomass)
    text(y = preds[750], x = 300, 
         labels = j)
    
    }
}

newdata = expand.grid(total_biomass = seq(0, 400, length.out = 1000),
                      mean_age = c(0, 10, 30, 70, 120, 240))

fia_plot_sub <- fia_plot2
fia_tree_sub <- fia_trees[fia_trees$PLT_CN %in% fia_plot_sub$CN, ]

plot(fia_plot_sub$cc ~ fia_plot_sub$total_biomass,
     main = "all", 
     xlim = c(0, 400),
     ylim = c(0, 1))
abline(h = c(0.4, 0.6))
for(j in unique(newdata$mean_age)){
  newdat <- newdata[newdata$mean_age == j, ]
  preds <- invlogit(predict(cc_no_cwhr_lm, newdat))
  lines( preds ~ newdat$total_biomass)
  text(y = preds[750], x = 300, 
       labels = j)
  
}



################################################################################
# Export models
write_rds(cc_lm, "canopy_cover_with_CWR_type_lm.RDS")
write_rds(cc_no_cwhr_lm, "canopy_cover_without_CWR_type_lm.RDS")


