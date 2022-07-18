# Functions to convert LANDIS -> CWHR

#*******************************************************************************
# helper functions


logit <- function(x) log(x/(1-x))
invlogit <- function(x) exp(x)/(1+exp(x))



#*******************************************************************************

classify_forest_type_rasters <- function(biomass_folder){
  # For models runs without the community output extension results, but that do have
  # the biomass by species outputs. 
  
}


#*******************************************************************************
#*

classify_forest_type_comm_output <- function(comm_table){
  #This function takes the location of the community output table from LANDIS ]
  # and returns a dataframe that maps the MapCode from the community output table
  # to a CWHR vegetation type code
  
  #The vegetation types and species have been reduced to match TCSI. This code
  #would need to be revised to accept new species and types. This would be a 
  #good place to make improvements. SF
  
  comm <- read.csv(comm_table)
  
  #replace shrub FTs with shrub
  #TODO make this flexible
  comm[comm$SpeciesName %in% c("FX_R_SEED", "NOFX_R_SEED", "NOFX_NOR_SEED"), "SpeciesName"] <- "Shrub"
  
  species <- unique(comm$SpeciesName)
  
  
  # sum biomass for each stand
  comm_matrix <- comm %>%
    dplyr::group_by(MapCode, SpeciesName) %>%
    summarise(biomass = sum(CohortBiomass, na.rm = TRUE)) %>%
    tidyr::pivot_wider(id_cols = MapCode, names_from = SpeciesName, 
                       values_from = biomass, values_fill = 0) %>%
    dplyr::ungroup()
  comm_matrix$total_biomass <- rowSums(comm_matrix[, -1])

  
  
  # classify
  
  comm_matrix$SMC <- rowSums(comm_matrix[, c("AbieConc", "CaloDecu", "PinuCont", 
                                             "PinuJeff", "PinuPond", "PinuLamb",
                                             "PseuMenz", "TorrCali")]) #"PinuWash"  "PinuBalf"
                                            
  comm_matrix$SCN <- rowSums(comm_matrix[, c("AbieMagn", "PinuAlbi", "PinuMont", "TsugMert")])
  comm_matrix$coni <-  rowSums(comm_matrix[, c("AbieConc", "AbieMagn", "CaloDecu", "JuniOcci",
                                               "PinuAlbi", "PinuCont", "PinuJeff", "PinuLamb",
                                               "PinuMont", "PinuPond", "PseuMenz", "PinuSabi",
                                               "TsugMert", "PinuMono", "TaxuBrev", "PinuAtte")]) #, "PinuBalf", "SequGiga", "PinuWash"
  comm_matrix$hard <- rowSums(comm_matrix[, c("AcerMacr", "AlnuRhom",  "ArbuMenz", #"AlnuRubr",
                                              "LithDens", "PopuTrem", "QuerChry", "QuerKell","QuerWisl",
                                              "QuerDoug", "QuerGarr", "CornNutt", "AescCali",
                                              "UmbeCali", "TorrCali")])
  comm_matrix <- comm_matrix %>%  # this is revised somewhat: hardwood has to make
                                  # up > 25% of the biomass; I added some equal signs
                                  # to prevent ties
    mutate(cover.coni = coni > 0 & coni > hard & coni > Shrub,
           cover.hard = hard > 0 & hard >= coni & (coni/total_biomass) < 0.25 & 
             (hard/total_biomass) > 0.25,
           cover.mixed = hard > 0 & coni > 0 & 
             hard >= coni & (coni/total_biomass) >= 0.25,
           cover.shrub = Shrub >= hard & Shrub >= coni)

#this is the dumbest thing I've ever written
comm_matrix$CWHR_type <- 
  ifelse(comm_matrix$cover.mixed, "mhc",
  ifelse(comm_matrix$cover.shrub, "mch", #chapparal/shrubland
  ifelse(comm_matrix$cover.hard & 
    comm_matrix$PopuTrem / comm_matrix$hard > 0.5, "asp",
  ifelse(comm_matrix$cover.hard & 
    ((comm_matrix$ArbuMenz + comm_matrix$LithDens + comm_matrix$QuerDoug + comm_matrix$QuerChry + 
        comm_matrix$QuerKell + comm_matrix$QuerWisl + comm_matrix$QuerGarr + comm_matrix$UmbeCali +
        comm_matrix$AescCali + comm_matrix$TorrCali + comm_matrix$PinuAtte) / 
       comm_matrix$hard) > 0.5,
    "mhw",
  ifelse(comm_matrix$cover.hard & 
    ((comm_matrix$AcerMacr + comm_matrix$AlnuRhom + comm_matrix$CornNutt )/comm_matrix$hard) > 0.5,
    "mri",
  ifelse(comm_matrix$cover.coni & comm_matrix$AbieConc/comm_matrix$coni > 0.5, "wfr",
  ifelse(comm_matrix$cover.coni & comm_matrix$AbieMagn/comm_matrix$coni > 0.5, "rfr",
  ifelse(comm_matrix$cover.coni & comm_matrix$PinuJeff/comm_matrix$coni > 0.5, "jpn",
  ifelse(comm_matrix$cover.coni & comm_matrix$PseuMenz/comm_matrix$coni > 0.5, "dfr",
  ifelse(comm_matrix$cover.coni & 
    comm_matrix$PinuPond/comm_matrix$coni > 0.5, "ppn",
  ifelse(comm_matrix$cover.coni & comm_matrix$PinuCont / comm_matrix$coni > 0.5, "lpn",
  ifelse(comm_matrix$cover.coni & comm_matrix$JuniOcci / comm_matrix$coni > 0.5, "juo",
  # ifelse(comm_matrix$PinuMono / comm_matrix$coni > 0.5, "pjn",       #added for TCSI
  ifelse(comm_matrix$PinuMono / comm_matrix$coni > 0.5, "juo",
  # ifelse(comm_matrix$cover.coni & comm_matrix$SequGiga / comm_matrix$coni > 0.5, "seg",
  ifelse(comm_matrix$cover.coni & 
    comm_matrix$AbieConc / comm_matrix$coni <= 0.5 & 
    comm_matrix$AbieMagn / comm_matrix$coni <= 0.5 &
    # comm_matrix$CaloDecu / comm_matrix$coni <= 0.5 &
    comm_matrix$JuniOcci / comm_matrix$coni <= 0.5 &
    comm_matrix$PinuAlbi / comm_matrix$coni <= 0.5 &
    comm_matrix$PinuCont / comm_matrix$coni <= 0.5 &
    comm_matrix$PinuJeff / comm_matrix$coni <= 0.5 &
    comm_matrix$PinuPond / comm_matrix$coni <= 0.5 &
    comm_matrix$PseuMenz / comm_matrix$coni <= 0.5 &
    # comm_matrix$TsugMert / comm_matrix$coni < 0.5 &
    # comm_matrix$PinuBalf / comm_matrix$coni < 0.5 &
    comm_matrix$SMC >= comm_matrix$SCN,
    "smc",
  ifelse(comm_matrix$cover.coni & 
    comm_matrix$AbieConc / comm_matrix$coni <= 0.5 & 
    comm_matrix$AbieMagn / comm_matrix$coni <= 0.5 &
    # comm_matrix$CaloDecu / comm_matrix$coni <= 0.5 &
    comm_matrix$JuniOcci / comm_matrix$coni <= 0.5 &
    comm_matrix$PinuAlbi / comm_matrix$coni <= 0.5 &
    comm_matrix$PinuCont / comm_matrix$coni <= 0.5 &
    comm_matrix$PinuJeff / comm_matrix$coni <= 0.5 &
    comm_matrix$PinuPond / comm_matrix$coni <= 0.5 &
    comm_matrix$PseuMenz / comm_matrix$coni <= 0.5 &
    # comm_matrix$TsugMert / comm_matrix$coni < 0.5 & 
    # comm_matrix$PinuBalf / comm_matrix$coni < 0.5 &
    # comm_matrix$SequGiga / comm_matrix$coni < 0.5 &
    comm_matrix$SMC < comm_matrix$SCN,
    "scn",
  ifelse(comm_matrix$cover.coni &
    (comm_matrix$PinuJeff + comm_matrix$PinuPond)/comm_matrix$coni > 0.5,
    "smc", #was "yps"
  ifelse(comm_matrix$cover.coni &
    (comm_matrix$PinuMont + comm_matrix$PinuAlbi) / comm_matrix$coni > 0.5,
    "smc", #was "wps"
  NA)))))))))))))))))           

test <- comm_matrix[is.na(comm_matrix$CWHR_type), ]

if(length(which(is.na(comm_matrix$CWHR_type))) > 0){
  message("Some MapCodes were not assigned CWHR types (", 
           length(which(is.na(comm_matrix$CWHR_type))),
           " sites)")
}

return(comm_matrix[, c("MapCode", "CWHR_type")])
}
