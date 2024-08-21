

library(basemaps)
library(terra)

Yiel_total <- template


set_defaults(Yiel_total, map_service = "esri", map_type = "world_hillshade")
x <- basemap_raster(Yiel_total, map_service = "esri", map_type = "world_hillshade")
x_terr <- rast(x)


library(ggplot2)
library(tidyterra)

ggplot() +
  geom_spatraster_rgb(data = x_terr) +
  geom_spatraster(data = Yiel_total)