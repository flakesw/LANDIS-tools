# Functions to convert LANDIS -> CWHR

#*******************************************************************************
# helper functions


logit <- function(x) log(x/(1-x))

invlogit <- function(x) exp(x)/(1+exp(x))

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

project_to_template <- function(input_raster, template){
  #function to project landis input rasters with no CRS to an input file
  
    #replace values of template with values from input raster
    out_raster <- template %>%
      `values<-`(values(input_raster))


  return(out_raster)
}

#*******************************************************************************
#* Classify community matrix into CWHR types
#*******************************************************************************

CWHR_type_from_comm_matrix <- function(comm_matrix){  

  #match possible species to those in the comm_matrix
  SMC_names_all <- c("AbieConc", "CaloDecu", "PinuCont", 
                     "PinuJeff", "PinuPond", "PinuLamb",
                     "PseuMenz", "TorrCali", "PinuWash", "PinuBalf") 
  SMC_names_reduced <- SMC_names_all[SMC_names_all %in% names(comm_matrix)]
  
  SCN_names_all <- c("AbieMagn", "PinuAlbi", "PinuMont", "TsugMert")
  SCN_names_reduced <- SCN_names_all[SCN_names_all %in% names(comm_matrix)]
  
  coni_names_all <- c("AbieConc", "AbieMagn", "CaloDecu", "JuniOcci",
                      "PinuAlbi", "PinuCont", "PinuJeff", "PinuLamb",
                      "PinuMont", "PinuPond", "PseuMenz", "PinuSabi",
                      "TsugMert", "PinuMono", "TaxuBrev", "PinuAtte",
                      "PinuBalf", "SequGiga", "PinuWash")
  coni_names_reduced <- coni_names_all[coni_names_all %in% names(comm_matrix)]
  
  hard_names_all <- c("AcerMacr", "AlnuRhom",  "ArbuMenz", #"AlnuRubr",
                      "LithDens", "PopuTrem", "QuerChry", "QuerKell","QuerWisl",
                      "QuerDoug", "QuerGarr", "CornNutt", "AescCali",
                      "UmbeCali", "TorrCali")
  hard_names_reduced <- hard_names_all[hard_names_all %in% names(comm_matrix)]
  
  mhw_names_all <- c("ArbuMenz", "LithDens", "QuerDoug", "QuerChry", 
                   "QuerKell", "QuerWisl", "QuerGarr"," UmbeCali",
                   "AescCali", "TorrCali", "PinuAtte")
  mhw_names_reduced <- mhw_names_all[mhw_names_all %in% names(comm_matrix)]
  
  mri_names_all <- c("AcerMacr", "AlnuRhom", "CornNutt")
  mri_names_reduced <- mri_names_all[mri_names_all %in% names(comm_matrix)]

  comm_matrix$SMC <- rowSums(comm_matrix[, names(comm_matrix) %in% SMC_names_reduced])
  comm_matrix$SCN <- rowSums(comm_matrix[, names(comm_matrix) %in% SCN_names_reduced])
  comm_matrix$coni <-  rowSums(comm_matrix[,  names(comm_matrix) %in% coni_names_reduced])
  comm_matrix$hard <- rowSums(comm_matrix[, names(comm_matrix) %in% hard_names_reduced])
  

  
  comm_matrix <- comm_matrix %>%  # this is revised somewhat: hardwood has to make
                                  # up > 25% of the biomass; I added some equal signs
                                  # to prevent ties
    mutate(cover.coni = coni > 0 & coni > hard & coni > Shrub,
           cover.hard = hard > 0 & hard >= coni & (coni/total_biomass) < 0.25 & 
             (hard/total_biomass) > 0.25,
           cover.mixed = hard > 0 & coni > 0 & 
             hard >= coni & (coni/total_biomass) >= 0.25,
           cover.shrub = Shrub >= hard & Shrub >= coni)
  
  class <- #set defualt values, to be replaced by more specific types. Helps prevent NAs
  ifelse(comm_matrix$cover.mixed, "mhc",
  ifelse(comm_matrix$cover.shrub, "mch", #chapparal/shrubland
  ifelse(comm_matrix$cover.coni, "smc",  #set default -- sierra mix con
  ifelse(comm_matrix$cover.hard, "mhw", NA)))) #set default for hardwood -- mixed hardwood       
         
  class <-
  ifelse(comm_matrix$cover.hard & 
    (comm_matrix$PopuTrem / comm_matrix$hard) >= 0.5, "asp", #if hardwood and majority aspen
  ifelse(comm_matrix$cover.hard & 
     (rowSums(comm_matrix[,  names(comm_matrix) %in% mri_names_reduced])/comm_matrix$hard) >= 0.5,
   "mri", #if hardwood and majority MRI species
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
   "scn", #if no dominant species with their own types (see below), then sites that are not SMC are SNC
  #specific types for sites dominated by certain species
  ifelse(comm_matrix$cover.coni & comm_matrix$AbieConc/comm_matrix$coni >= 0.5, "wfr",
  ifelse(comm_matrix$cover.coni & comm_matrix$AbieMagn/comm_matrix$coni >= 0.5, "rfr",
  ifelse(comm_matrix$cover.coni & comm_matrix$PinuJeff/comm_matrix$coni >= 0.5, "jpn",
  ifelse(comm_matrix$cover.coni & comm_matrix$PseuMenz/comm_matrix$coni >= 0.5, "dfr",
  ifelse(comm_matrix$cover.coni & comm_matrix$PinuPond/comm_matrix$coni >= 0.5, "ppn",
  ifelse(comm_matrix$cover.coni & comm_matrix$PinuCont / comm_matrix$coni >= 0.5, "lpn",
  ifelse(comm_matrix$cover.coni & comm_matrix$JuniOcci / comm_matrix$coni >= 0.5, "juo",
  # ifelse(comm_matrix$PinuMono / comm_matrix$coni > 0.5, "pjn",       #pinyon-juniper woodlands, if desired
  ifelse(comm_matrix$cover.coni & 
             (comm_matrix$PinuMono / comm_matrix$coni) >= 0.5, "juo",
  class))))))))))) #write original default values if no conditions apply

  return(class)
}


#*******************************************************************************
#* Stand seral stage
#* *****************************************************************************

calculate_age_ht_dia <- function(comm = comm,
                                      dia_regression_rds_loc = "linear_models_diam_from_age.RDS",
                                      dia_regression_rds_no_sp_loc = "linear_models_diam_from_age_no_sp.RDS",
                                      ht_regression_rds_loc = "linear_models_ht_from_age.RDS",
                                      ht_regression_rds_no_sp_loc = "linear_models_ht_from_age_no_sp.RDS"
                                      ){
  
  #import regression equations fit from FIA data, see regression_age_diam.R
  age_to_dia_reg <- readRDS(dia_regression_rds_loc)
  age_to_dia_reg_no_sp <- readRDS(dia_regression_rds_no_sp_loc)
 
  age_to_ht_reg <- readRDS(ht_regression_rds_loc)
  age_to_ht_reg_no_sp <- readRDS(ht_regression_rds_no_sp_loc)
  
  sp_in_comm <- unique(comm$SpeciesName)
  
  comm$biomass <- comm$CohortBiomass * 180*180 / (10^6) #convert from g/m2 to Mg (per site)
  
  #fill in 
  #this takes about a minute
  for(i in 1:length(sp_in_comm)){
    sp = sp_in_comm[i]
    if(sp %in% age_to_dia_reg$SpeciesName){
     
      comm[comm$SpeciesName == sp, "dia"] <- exp(predict(age_to_dia_reg$model[age_to_dia_reg$SpeciesName == sp][[1]],
                                                      newdata = comm[comm$SpeciesName == sp, ]))
      comm[comm$SpeciesName == sp, "ht"] <- exp(predict(age_to_ht_reg$model[age_to_ht_reg$SpeciesName == sp][[1]],
                                                         newdata = comm[comm$SpeciesName == sp, ]))
      
      
    } else{
      comm[comm$SpeciesName == sp, "dia"] <- exp(predict(age_to_dia_reg_no_sp, newdata = comm[comm$SpeciesName == sp, ]))
      comm[comm$SpeciesName == sp, "ht"] <- exp(predict(age_to_ht_reg_no_sp, newdata = comm[comm$SpeciesName == sp, ]))
     
    }
  }
  
  
  #aggregate to plot
  plot_age_dia <- comm %>%
    # filter(dia > 5) %>%
    dplyr::group_by(MapCode) %>%
    dplyr::summarise(mean_age = weighted.mean(TOTAGE, CohortBiomass), #age calculated here doesn't actually affect the seral stages;
                                                                      #that's all determined by the regression equations and how
                                                                      #diameter is aggregated. It does affect CC calculoations.
                     mean_age_unweighted = mean(TOTAGE),
                     weighted_mean_diam = weighted.mean(dia, CohortBiomass),
                     mean_diam = mean(dia),
                     root_average = sqrt(weighted.mean((dia)^2, CohortBiomass)),
                     max_dia = max(dia),
                     max_ht = max(ht),
                     ht_95 = quantile(ht, .95) / 3.28 #convert to meters
                     )
  
  return(plot_age_dia)
  }


#*******************************************************************************
#* Canopy class
#*******************************************************************************
calculate_canopy_cover <- function(comm_matrix, 
                                   can_regresison_rds_loc = "canopy_cover_with_CWR_type_lm.RDS",
                                   can_regression_rds_no_sp_loc = "canopy_cover_without_CWR_type_lm.RDS"){

 
  can_model <- readRDS(can_regresison_rds_loc)
  can_model_no_cwhr <- readRDS("canopy_cover_without_CWR_type_lm.RDS")

  # cwhr_types_with_mod <- unique(can_model$model$CWHR_type)
  cwhr_types_with_mod <- can_model$xlevels$CWHR_type
  all_cwhr_types <- unique(comm_matrix$CWHR_type)
  all_cwhr_types <- all_cwhr_types[!is.na(all_cwhr_types)]
  
  comm_matrix$cc <- NA
  
  comm_matrix$total_biomass <- comm_matrix$total_biomass/100 #convert g/m2 to tonnes/ha
  comm_matrix$seral_stage <- cut(comm_matrix$weighted_mean_diam, breaks = c(0,1,6,11,24,1000), labels = FALSE)
  
  for(i in 1:length(all_cwhr_types)){
    
    # this has some dumb workarounds to avoid problems that arise when there are 
    # NAs for CWHR codes. It causes there to be different lengths for comm_matrix
    # and the predictions. I'd love to make this less dumb.
    
    cwhr <- all_cwhr_types[i]
    newdata <- subset(comm_matrix, CWHR_type == cwhr)
    newdata$cc <- NA
    if(cwhr %in% cwhr_types_with_mod){
      newdata$cc <- invlogit(predict(can_model, newdata = newdata))[, 1]
    } else{
      newdata$cc <- invlogit(predict(can_model_no_cwhr, newdata = newdata))[, 1]
    }
    
    comm_matrix[match(newdata$MapCode, comm_matrix$MapCode), ]$cc <- newdata$cc
  }
  
  return(comm_matrix$cc)
}

#*******************************************************************************
create_canopy_cover_raster <- function(output_folder,
                                      comm_matrix,
                                      comm_raster,
                                      file_prefix = "",
                                      timestep,
                                      class = TRUE){

  cc_breaks <- c(0,0.1,0.25,0.4,0.6,1)
  comm_matrix$cc_class <- cut(comm_matrix$cc, cc_breaks, labels = FALSE) - 1

  if(class == TRUE){
    message("   Writing canopy cover class raster to ", paste0(output_folder, file_prefix, "cc_class-", timestep, ".tif"))
    cc_raster <- terra::classify(comm_raster, rcl = comm_matrix[, c("MapCode", "cc_class")]) %>%
      terra::clamp(lower = 0, upper = 4, values = FALSE)
    writeRaster(cc_raster, paste0(output_folder, file_prefix, "cc_class-", timestep, ".tif"), overwrite = TRUE)
    message("   Done writing canopy cover class raster")
  } else if(class == FALSE){
    message("   Writing canopy cover continuous values raster to ", paste0(output_folder, file_prefix, "cc_continuous-", timestep, ".tif"))
    cc_raster <- terra::classify(comm_raster, rcl = comm_matrix[, c("MapCode", "cc")]) %>%
      terra::clamp(lower = 0, upper = 1, values = FALSE)
    writeRaster(cc_raster, paste0(output_folder, file_prefix, "cc_continuous-", timestep, ".tif"), overwrite = TRUE)
    message("   Done writing canopy cover continuous values raster")
  }
}

#*******************************************************************************
create_lai_raster <- function(output_folder,
                              comm_matrix,
                              comm_raster,
                              file_prefix = "",
                              timestep,
                              lai_name
                              ){
  message("   Writing LAI raster for ", lai_name, " to ", 
          paste0(output_folder, file_prefix, lai_name, "-", timestep, ".tif"))
  lai_raster <- terra::classify(comm_raster, rcl = comm_matrix[, c("MapCode", lai_name)]) %>%
    terra::clamp(lower = min(comm_matrix[[lai_name]], na.rm = TRUE), upper = max(comm_matrix[[lai_name]], na.rm = TRUE), values = FALSE)
  writeRaster(lai_raster, paste0(output_folder, file_prefix, lai_name, "-", timestep, ".tif"), overwrite = TRUE)
  message("   Done writing LAI raster")
}


#*******************************************************************************
#* Process community input files and create raster outputs
#*******************************************************************************
process_DHSVM_layers_and_write_rasters <- function(landis_folder,
                                           template,
                                           timestep,
                                           output_folder,
                                           prefix,
                                           dia_regression_rds_loc,
                                           dia_regression_rds_no_sp_loc,
                                           can_regresison_rds_loc,
                                           can_regression_rds_no_sp_loc,
                                           shrub_regression_loc,
                                           class = FALSE){
  
  message(paste0("Beginning processing LANDIS run ", landis_folder))
  
  comm_table <- paste0(landis_folder, "community-input-file-", timestep, ".csv")
  
  #otherwise the warning is converted to error
  options(warn = 1)
  comm_raster <- terra::rast(paste0(landis_folder, "output-community-",timestep, ".img")) %>%
    project_to_template(template)
  NAflag(comm_raster) <- NaN
  values(comm_raster)[values(comm_raster) == 0] <- NaN

  #make the community matrix and wrangle some
  
  comm <- read.csv(comm_table) %>%
    filter(complete.cases(.))
  
  #replace shrub FTs with shrub
  #TODO make this flexible
  comm[comm$SpeciesName %in% c("FX_R_SEED", "NOFX_R_SEED", "NOFX_NOR_SEED"), "SpeciesName"] <- "Shrub"
  species <- unique(comm$SpeciesName)
  comm <- rename(comm, TOTAGE = CohortAge)
  
  
  # sum biomass for each stand
  comm_matrix <- comm %>%
    dplyr::group_by(MapCode, SpeciesName) %>%
    summarise(biomass = sum(CohortBiomass, na.rm = TRUE)) %>%
    tidyr::pivot_wider(id_cols = MapCode, names_from = SpeciesName, 
                       values_from = biomass, values_fill = 0) %>%
    dplyr::ungroup()
  comm_matrix$total_biomass <- rowSums(comm_matrix[, -1]) # in g m-2
  comm_matrix$shrub_proportion <- comm_matrix$Shrub / comm_matrix$total_biomass
  comm_matrix$overstory_proportion = 1-comm_matrix$shrub_proportion
  
  #*******************************************************************************
  #Forest type classification
  #*******************************************************************************
  
  # this takes about 5 minutes to run. It might take a lot of RAM if you're doing
  # a big landscape. A place to improve might be to break this into chunks and do 
  # each chunk individually, then mosaic back together.
  
  
  #TODO update these
  class_table <- as.data.frame(
    rbind(c("asp", 1),
          c("mhw", 2),
          c("mri", 3),
          c("wfr", 4),
          c("rfr", 5),
          c("jpn", 6),
          c("ppn", 7),
          c("dfr", 8),
          c("mhc", 9),
          c("lpn", 10),
          c("smc", 11),
          c("scn", 12),
          c("juo", 13),
          c("mch", 14)
    )) %>%
    dplyr::rename(type = V1,
                  num = V2)
  
  comm_matrix$CWHR_type <- CWHR_type_from_comm_matrix(comm_matrix)
  
  if(length(which(is.na(comm_matrix$CWHR_type))) > 0){
    message("Some MapCodes were not assigned CWHR veg types (", 
            length(which(is.na(comm_matrix$CWHR_type))),
            " sites)")
  }
  
  #*******************************************************************************
  # Stand seral stage
  #*******************************************************************************

  #used to improve estimates of canopy cover

  plot_age_dia <- calculate_age_ht_dia(comm = comm,
                                       dia_regression_rds_loc = "linear_models_diam_from_age.RDS",
                                       dia_regression_rds_no_sp_loc = "linear_models_diam_from_age_no_sp.RDS",
                                       ht_regression_rds_loc = "linear_models_ht_from_age.RDS",
                                       ht_regression_rds_no_sp_loc = "linear_models_ht_from_age_no_sp.RDS")

  comm_matrix <- left_join(comm_matrix, plot_age_dia, by = "MapCode")
  
  #*******************************************************************************
  # Stand canopy cover
  # #*******************************************************************************
  comm_matrix$cc <- calculate_canopy_cover(comm_matrix = comm_matrix,
                                           can_regresison_rds_loc = "canopy_cover_with_CWR_type_lm.RDS",
                                           can_regression_rds_no_sp_loc = "canopy_cover_without_CWR_type_lm.RDS")

  create_canopy_cover_raster(output_folder = output_folder,
                             comm_matrix = comm_matrix,
                             comm_raster = comm_raster,
                             file_prefix = prefix,
                             timestep = timestep,
                             class = FALSE
  )
  
  #*******************************************************************************
  #* Estimate shrub cover
  #*******************************************************************************
  shrub_model <- readRDS(shrub_regression_loc)
  cwhr_update <- ifelse(comm_matrix$CWHR_type == "mri", "mhc", comm_matrix$CWHR_type)
  
  comm_matrix$shrub_cover <- predict(shrub_model, newdata = data.frame(biomass = comm_matrix$total_biomass,
                                                                 cc = comm_matrix$cc,
                                                                 mean_age = comm_matrix$mean_age_unweighted,
                                                                 CWHR_type = cwhr_update))

  
  #*******************************************************************************
  # LAI -- split by overstory and understory
  #*******************************************************************************

  necn_folder <- paste0(landis_folder, "NECN/")
  lai_files <- list.files(necn_folder, pattern = "LAI")
  lai_files <- gsub("-", "", lai_files)
  lai_years <- readr::parse_number(lai_files)

  timestep2 <- timestep
  if(timestep == 0){
    timestep2 <- min(lai_years)
    message(paste0("Timestep is lower than first NECN output; setting timestep = ", timestep2))
  } else if(!(timestep %in% lai_years)){
    error_flag <- TRUE
    message("Timestep has no LAI layer -- check that timesteps align for community
            output file and NECN")
  }
  if(error_flag){
    next()
  }

  lai_raster <- terra::rast(paste0(necn_folder, "LAI-", timestep2, ".img"))
  # plot(lai_raster)
  lai_vals <- data.frame(cbind(terra::values(comm_raster), terra::values(lai_raster)))
  names(lai_vals) <- c("MapCode", "LAI")

  comm_matrix <- left_join(comm_matrix, lai_vals,
                           by = "MapCode")
  comm_matrix$LAI_shrub <- comm_matrix$LAI*comm_matrix$shrub_proportion #assign LAI weighted by biomass
  comm_matrix$LAI_tree <- comm_matrix$LAI*comm_matrix$overstory_proportion

  create_lai_raster(output_folder = output_folder,
                    comm_matrix = comm_matrix,
                    comm_raster = comm_raster,
                    file_prefix = scenario_name,
                    timestep = timestep,
                    lai_name = "LAI_shrub")
  create_lai_raster(output_folder = output_folder,
                    comm_matrix = comm_matrix,
                    comm_raster = comm_raster,
                    file_prefix = scenario_name,
                    timestep = timestep,
                    lai_name = "LAI_tree")

  #*******************************************************************************
  #* Create CWHR DHSVM Codes
  #* **********************************************************************
  # if(exists("HS_reference")){
  #   # Check if all the CWHR codes match the allowed types in reference
  #   # and replace with valid codes
  # 
  #   comm_matrix <- comm_matrix %>%
  #     mutate(in_ref = CWHR_code %in% HS_reference) %>%
  #     mutate(CWHR_code = case_when(CWHR_type == "mch" ~ 1400,
  #                                  !in_ref & seral_class == 1 ~ as.numeric(paste0(num, seral_class, 0)),
  #                                  !in_ref & seral_class > 1 & cc_class == 0 ~ as.numeric(paste0(num, seral_class, 1)),
  #                                  TRUE ~ CWHR_code))
  #   comm_matrix$in_ref <- comm_matrix$CWHR_code %in% HS_reference
  # 
  #   #are they all fixed now?
  #   if(sum(!comm_matrix$in_ref) > 0){
  #     message("Some CWHR codes not found in reference (", sum(!comm_matrix$in_ref), " sites)")
  #   }
  # }
  
  reference_types <- paste0(1:14, rep(c("01","02"), each = 14)) #allowable veg type codes
  # comm_bak <- comm_matrix #remove
  # comm_matrix <- comm_bak #remove
  
  # print(comm_matrix[is.na(comm_matrix$CWHR_type), ], width = Inf)
  
  
  comm_matrix <- comm_matrix %>% 
    dplyr::inner_join(class_table, by = c("CWHR_type" = "type")) %>%
    mutate(DHSVM_code = as.numeric(paste0(as.character(num), ifelse(shrub_proportion > 0.01 | shrub_cover > 0.2, "02", "01"))))%>%
    # mutate(in_ref = DHSVM_code %in% reference_types) %>%
    mutate(DHSVM_code = case_when(shrub_cover > 0.2 & cc < 0.1 ~ 1402,
                                  TRUE ~ DHSVM_code)) %>%
    mutate(DHSVM_code = case_when(CWHR_type == "mch" ~ 1402,
                                  # !in_ref ~ NA,
                                  TRUE ~ DHSVM_code)) 
  # table(comm_matrix$CWHR_type)
  # table(comm_matrix$DHSVM_code)
  
  message("   Writing veg ID raster to ", paste0(output_folder, prefix, "veg-ID-", timestep, ".tif"))
  
  cwhr_raster <- terra::classify(comm_raster, rcl = comm_matrix[, c("MapCode", "DHSVM_code")], others = NA)
  
  writeRaster(cwhr_raster, paste0(output_folder, prefix, "veg-ID-", timestep, ".tif"), overwrite = TRUE)
  message("   Done writing veg ID raster")
  
  
  #*****************************************************************************
  #* Write Canopy height raster
  #*****************************************************************************
  
  message("   Writing canopy height raster to ", paste0(output_folder, prefix, "can_ht-", timestep, ".tif"))
  
  #set height for chaparral
  comm_matrix <-  mutate(comm_matrix, ht_95 = case_when(DHSVM_code == 1402 ~ 3,
                             TRUE ~ ht_95))
  
  ht_raster <- terra::classify(comm_raster, rcl = comm_matrix[, c("MapCode", "ht_95")]) %>%
    terra::clamp(lower = min(unique(comm_matrix$ht_95)),
                 upper = max(unique(comm_matrix$ht_95)),
                 values = FALSE)
  writeRaster(ht_raster, paste0(output_folder, prefix, "can-ht-", timestep, ".tif"), overwrite = TRUE)
  message("   Done writing height raster")
  
  gc()
  
  message(paste0("Done processing LANDIS run ", landis_folder))
  
}
