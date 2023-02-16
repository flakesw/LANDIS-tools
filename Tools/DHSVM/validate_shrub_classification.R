#compare calculated chaparral to observed
library(terra)
library(tidyverse)


evt <- rast("./LF2016_EVT_200_CONUS/LC16_EVT_200.tif")
plot(evt)

dhsvm <- rast("E:/TCSI LANDIS/test/projected/Scenario1 - cnrm - Run 1 -  veg-ID- 0 .tif.tif") %>%
  terra::project(evt, method = "near", mask = TRUE)

plot(dhsvm)

evt_subset <- terra::mask(evt, dhsvm) %>% crop(dhsvm)
plot(evt_subset)
table(values(evt_subset))


woodland_types <- c(7008:7014, 7019, 7029,
                    7031, 7033, 7034, 7044,
                     7062, 7112, 7114, 7151, 7264,
                    7265, 9003, 9022)

open_types <- c(7067, 7071, 7079, 7080:7099, 
                7103, 7105, 7125:7137, 7145,
                7146, 7152, 7153, 7191, 7192, 7195,
                7196, 7197, 7198, 7199, 7234, 7292:7300, 7662,
                7903, 7904, 7913, 7914,  
                7923, 7924, 7928, 7929, 7943, 7944, 
                7960:7998, 9005, 9006, 9008, 9010, 9011,
                9017, 9062, 9125:9134, 9213, 9272, 9301, 
                9629, 9633, 9827, 9828, 9829)

understory_types_dhsvm <- c(102, 202, 302, 402, 502, 602, 702, 802, 902, 1002, 1102, 1202, 1302, 1402)

evt_at_dhsvm_shrub <- values(terra::mask(evt_subset, dhsvm, maskvalues = 1402, inverse = TRUE)) %>%
  `[`(!is.na(.))
table(evt_at_dhsvm_shrub %in% open_types) #>50% open at sites classified as shrubs
table((values(evt_subset) %>% `[`(!is.na(.))) %in% open_types) #only 23% open in whole area



dhsvm_at_evt_shrub <- values(terra::mask(dhsvm, evt_subset, maskvalues = open_types, inverse = TRUE)) %>%
  `[`(!is.na(.))
table((dhsvm_at_evt_shrub))


