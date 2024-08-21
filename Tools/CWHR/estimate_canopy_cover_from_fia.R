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
#TODO fix SMC and SCN

#fix a problem with species being absent
comm_matrix$PinuMono <- 0
comm_matrix$CWHR_type <- CWHR_type_from_comm_matrix(comm_matrix)

#add CWHR types to fia_plot dataframe
fia_plot <- left_join(fia_plot, comm_matrix[, c("PLT_CN", "CWHR_type")], by = c("CN" = "PLT_CN"))

################################################################################
# Aggregate FIA data

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
  mutate(age_all = ifelse(!is.na(TOTAGE), TOTAGE, age_est))

# aggregate trees to plot. We need the weighted average age and the total plot
# biomass in order to regress CC ~ age*CWHR*biomass

tree_summary <- fia_trees %>%
  filter(DIA > 0) %>%
  group_by(PLT_CN) %>%
  summarise(total_biomass = sum(DRYBIO_TOTAL * TPA_UNADJ) * 2.47, #convert Mg per acre to Mg per ha,
            mean_age = weighted.mean(age_est, DRYBIO_TOTAL, na.rm = TRUE),
            total_n_trees = sum(TPA_UNADJ, na.rm = TRUE),
            ba = sum(DIA^2 * 0.005454, na.rm = TRUE), #in square feet
            qmd = sqrt(ba/((0.005454*total_n_trees)))) #in inches

fia_plot2 <- left_join(fia_plot, tree_summary, by = c("CN" = "PLT_CN")) %>%
  filter(cc > 0) %>%
  mutate(cc = cc/100)

dia_breaks <- c(0, 1, 6, 11, 24, 1000)
fia_plot2$seral_stage <- cut(fia_plot2$qmd, dia_breaks, labels = FALSE)


#exploratory plots
plot(cc ~ total_biomass, data = fia_plot2)
plot(cc ~ sqrt(total_biomass), data = fia_plot2) #strong relationship
plot(cc ~ log(total_biomass), data = fia_plot2)
plot(logit(cc) ~ sqrt(total_biomass), data = fia_plot2) #strong and linear relationship
plot(log(cc) ~ log(total_biomass), data = fia_plot2)
plot(cc ~ mean_age, data = fia_plot2) #little relationship
boxplot(cc ~ CWHR_type, data = fia_plot2)
boxplot(cc ~ seral_stage, data = fia_plot2)
boxplot(total_biomass ~ seral_stage, data = fia_plot2)

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
# cc_lm <- lm(logit(cc) ~ sqrt(total_biomass)*mean_age*CWHR_type, data = fia_plot_remove_rare)
# summary(cc_lm)

# cc_lm <- earth::earth(logit(cc) ~ sqrt(total_biomass)*mean_age*CWHR_type, 
#                       data = fia_plot_remove_rare[!is.na(fia_plot_remove_rare$total_biomass) &
#                                                     !is.na(fia_plot_remove_rare$mean_age) & 
#                                                     !is.na(fia_plot_remove_rare$CWHR_type), ])
cc_lm <- earth::earth(logit(cc) ~ sqrt(total_biomass)*seral_stage*CWHR_type,
                      data = fia_plot_remove_rare[!is.na(fia_plot_remove_rare$total_biomass) &
                                                    !is.na(fia_plot_remove_rare$mean_age) &
                                                    !is.na(fia_plot_remove_rare$CWHR_type), ])


summary(cc_lm)

cc_no_cwhr_lm <- earth::earth(logit(cc) ~ sqrt(total_biomass)*seral_stage, 
                    data = fia_plot2[!is.na(fia_plot2$total_biomass) &
                                     !is.na(fia_plot2$mean_age) & 
                                     !is.na(fia_plot2$CWHR_type), ])
summary(cc_no_cwhr_lm)
# plot(effects::allEffects(cc_no_cwhr_lm))

cwhr_types <- unique(fia_plot_remove_rare$CWHR_type)[!is.na(unique(fia_plot_remove_rare$CWHR_type))]


# Make plots of cc ~ biomass ---------------------------------------------------

# layout(mat = matrix(c(1:12), nrow = 4, ncol = 3, byrow = TRUE), height = 1, width = 1)

svg(filename = "cc_regressions.svg",
    width = 7, height = 10, pointsize = 9)

col_ramp <- c('#feb24c','#a1dab4','#41b6c4','#2c7fb8','#253494')
titles <- c("Chaparral", "Lodgepole pine", "Ponderosa pine", "Western juniper",
            "Sierra mixed conifer", "Jeffrey pine", "White fir", "Subalpine mixed conifer",
            "Red fir", "Douglas-fir", "Mixed hardwood", "All vegetation types")

par(mfrow = c(4,3),
    mar = c(2,2,1,1),
    oma = c(3,3,0,0))

newdata = expand.grid(total_biomass = seq(0, 1000, length.out = 1000),
                      # mean_age = c(0, 10, 30, 70, 120, 240),
                      seral_stage = c(1,2,3,4,5),
                      CWHR_type = cwhr_types)

for(i in 1:length(cwhr_types)){
  fia_plot_sub <- fia_plot_remove_rare[fia_plot_remove_rare$CWHR_type == cwhr_types[i], ]
  fia_tree_sub <- fia_trees[fia_trees$PLT_CN %in% fia_plot_sub$CN, ]
  
  plot(fia_plot_sub$cc ~ fia_plot_sub$total_biomass,
       col = col_ramp[fia_plot_sub$seral_stage],
       bg = col_ramp[fia_plot_sub$seral_stage],
       pch = 21,
       main = titles[i], 
       xlim = c(0, 400),
       ylim = c(0, 1),
       xlab = "",
       ylab = "",
       xaxt = "n",
       yaxt = "n")
  abline(h = c(0.4, 0.6))
  
  axis(1, at = c(0,100,200,300,400), labels = FALSE)
  axis(2, at = c(0,.2,.4,.6,.8,1), labels = FALSE)
  if(i %in% c(1,4,7,10)) axis(2)
  if(i %in% c(10,11)) axis(1)
  
  for(j in unique(newdata$seral_stage)){
    newdat <- newdata[newdata$CWHR_type == cwhr_types[i] &
                                  newdata$seral_stage == j, ]
    preds <- invlogit(predict(cc_lm, newdat))
    lines( preds ~ newdat$total_biomass, col = col_ramp[j])
    text(y = preds[300], x = 300, 
         labels = j)
    
    }
}


newdata = expand.grid(total_biomass = seq(0, 400, length.out = 1000),
                      seral_stage = c(1,2,3,4,5))

fia_plot_sub <- fia_plot2[!is.na(fia_plot2$total_biomass) &
                            !is.na(fia_plot2$mean_age) & 
                            !is.na(fia_plot2$CWHR_type), ]

plot(fia_plot_sub$cc ~ fia_plot_sub$total_biomass,
     col = col_ramp[fia_plot_sub$seral_stage],
     bg = col_ramp[fia_plot_sub$seral_stage],
     pch = 21,
     main = "All vegetation types", 
     xlim = c(0, 400),
     ylim = c(0, 1),
     mar = c(0,0,0,0),
     oma = c(0,0,0,0),
     xlab = "",
     ylab = "",
     xaxt = "n",
     yaxt = "n")
abline(h = c(0.4, 0.6))

axis(1, at = c(0,100,200,300,400), labels = FALSE)
axis(2, at = c(0,.2,.4,.6,.8,1), labels = FALSE)
axis(1)

for(j in unique(newdata$seral_stage)){
  newdat <- newdata[newdata$seral_stage == j, ]
  preds <- invlogit(predict(cc_no_cwhr_lm, newdat))
  lines(preds ~ newdat$total_biomass,
         col = col_ramp[j])
  text(y = preds[750], x = 300, 
       labels = j)
}

mtext(text = expression("Site biomass (Mg ha"^-1*")"),
      side = 1,
      line = 1,
      outer = TRUE)
mtext(text = expression("Canopy cover (proportion)"),
      side = 2,
      line = 1,
      outer = TRUE)




dev.off()

################################################################################
# Export models
write_rds(cc_lm, "canopy_cover_with_CWR_type_lm.RDS")
write_rds(cc_no_cwhr_lm, "canopy_cover_without_CWR_type_lm.RDS")


