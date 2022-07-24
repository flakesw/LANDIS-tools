

##########
# Evaluate canopy cover
library(terra)
library(sf)
sf::sf_use_s2(FALSE)


###########

mask <- terra::rast("C:/Users/Sam/Documents/Research/TCSI conservation finance/Models/Inputs/masks_boundaries/mask.tif") 

boundary <- sf::st_read("C:/Users/Sam/Documents/Research/TCSI conservation finance/Models/Inputs/masks_boundaries/tcsi_area_shapefile/TCSI_v2.shp") %>%
  sf::st_make_valid()
sf::st_write(boundary, "TCSI_boundary_valid.shp")

nlcd <- terra::rast("D:/Data/nlcd_2016_treecanopy_2019_08_31/nlcd_2016_treecanopy_2019_08_31.img")
crs(nlcd) <- "EPSG:5070"
# cc_raster <- terra::rast("initial_canopy_cover_modeled.tif")
cc_raster <- project_to_template(cc_raster_cont, mask)

backcast <- terra::rast("cancov_MIROC5_85_Scenario1_1_v2_0.tif")

mean(values(backcast), na.rm = TRUE)
mean(values(cc_raster), na.rm = TRUE)
mean(values(nlcd2), na.rm = TRUE)


mask <- terra::project(mask, "EPSG:5070")
cc_raster <- terra::project(cc_raster, mask)

# terra::writeRaster(cc_raster, "new_initial_cc.tif")

# boundary <- boundary %>%
#   sf::st_transform(crs = st_crs(nlcd))
nlcd2 <- nlcd  %>%
  crop(cc_raster) %>%
  resample(cc_raster) %>%
  mask(cc_raster)
values(nlcd2) <- as.numeric(values(nlcd2))/100
plot(nlcd2)
vals_nlcd <- values(nlcd2)
vals_nlcd <- vals_nlcd[vals_nlcd> 0 & !is.na(vals_nlcd)]

vals_cc <- values(cc_raster)
vals_cc <- vals_cc[vals_cc> 0 & !is.na(vals_cc)]

plot(density(vals_nlcd), col = "#66C2A5", lwd = 2,
     xlab = "Canopy Cover",
     ylab = "Density",
     main = "Modeled canopy cover vs NLCD canopy cover")
text(x = 0.8, y = 3.3, col = "#66C2A5", labels = c("NLCD\nCanopy Cover"))
lines(density(vals_cc), col = "#FC8D62", lwd = 2)
text(x = 0.45, y = 3, col = "#FC8D62", labels = c("LANDIS initial\ncanopy cover"))

hist(vals_cc, breaks = c(0, .1, .25, .4, .6, 1), col = addTrans("#FC8D62", 120), lwd = 2)
hist(vals_nlcd, breaks = c(0, .1, .25, .4, .6, 1), col = addTrans("#66C2A5", 120), lwd = 2, add = TRUE)
abline(v = c(0.4, 0.6))

mean(values(nlcd2), na.rm = TRUE)
mean(values(cc_raster), na.rm = TRUE)
