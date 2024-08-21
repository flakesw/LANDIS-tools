
################################################################################
# Estimate diameter:age regressions
addTrans <- function(color,trans)
{
  # This function adds transparancy to a color.
  # Define transparancy with an integer between 0 and 255
  # 0 being fully transparant and 255 being fully visable
  # Works with either color and trans a vector of equal length,
  # or one of the two of length 1.
  
  if (length(color)!=length(trans)&!any(c(length(color),length(trans))==1)) stop("Vector lengths not correct")
  if (length(color)==1 & length(trans)>1) color <- rep(color,length(trans))
  if (length(trans)==1 & length(color)>1) trans <- rep(trans,length(color))
  
  num2hex <- function(x)
  {
    hex <- unlist(strsplit("0123456789ABCDEF",split=""))
    return(paste(hex[(x-x%%16)/16+1],hex[x%%16+1],sep=""))
  }
  rgb <- rbind(col2rgb(color),trans)
  res <- paste("#",apply(apply(rgb,2,num2hex),2,paste,collapse=""),sep="")
  return(res)
}



#these are needed to convert LANDIS cohorts to FIA world

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

#fit models and plot output

library(earth)
earth_regressions <- fia_trees_ca %>%
  dplyr::filter(!is.na(DIA) & !is.na(TOTAGE) & !is.na(SPCD)) %>%
  filter(SpeciesName %in% species_class) %>%
  dplyr::group_by(SpeciesName) %>%
  filter(n() > 40) %>%
  dplyr::do(model = earth(log(DIA) ~ log(TOTAGE), data = .))

no_sp_regression <- earth(log(DIA) ~ log(TOTAGE), 
                          data = fia_trees_ca[!is.na(fia_trees_ca$TOTAGE) & 
                                                !is.na(fia_trees_ca$DIA), ])


#make figure
png(filename = "age_diam_figures.png",
    width = 7, height = 10,
    units = "in",
    res = 300,
    pointsize = 10)

par(mfrow = c(5, 5),
    mar = c(2,2,1,1),
    oma = c(4,4,0,0))

for(i in 1:nrow(earth_regressions)){
  fia_sub <- fia_trees_ca[fia_trees_ca$SpeciesName == earth_regressions$SpeciesName[i], ]

  plot(DIA ~ TOTAGE, 
       data = fia_sub,
       main = earth_regressions$SpeciesName[i],
       col = addTrans("blue", 20),
       bg = addTrans("blue", 20),
       pch = 21,
       xlab = "",
       ylab = "")
  
  newdata = data.frame(TOTAGE = seq(0, max(fia_sub$TOTAGE, na.rm = TRUE), length.out = 1000))
  preds <- exp(predict(earth_regressions$model[i][[1]], newdata = newdata)) 
            
                 # var(residuals(earth_regressions$model[i][[1]])))
  lines(preds ~ newdata$TOTAGE, lwd = 2)
  abline(h = 24)
}
plot(DIA ~ TOTAGE, 
     data = fia_trees_ca,
     main = "All trees",
     col = addTrans("blue", 10),
     bg = addTrans("blue", 10),
     pch = 21,
     xlab = "",
     ylab = "")

newdata = data.frame(TOTAGE = seq(0, max(fia_trees_ca$TOTAGE, na.rm = TRUE), length.out = 1000))
preds <- exp(predict(no_sp_regression, newdata = newdata)) 

# var(residuals(earth_regressions$model[i][[1]])))
lines(preds ~ newdata$TOTAGE, lwd = 2)
abline(h = 24)

mtext(text = "Tree age (yr)", side = 1, line = 1.5, outer = TRUE, cex = 1.5)
mtext(text = "Diameter (cm)", side = 2, line = 1.5, outer = TRUE, cex = 1.5)


dev.off()




write_rds(earth_regressions, "linear_models_diam_from_age.RDS")
write_rds(no_sp_regression, "linear_models_diam_from_age_no_sp.RDS")