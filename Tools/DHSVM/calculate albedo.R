#create albedo timeseries for each veg type

dhsvm <- rast("E:/TCSI LANDIS/test/projected/Scenario1 - cnrm - Run 1 -  veg-ID- 0 .tif.tif") 
dhsvm


albedo <- rast("./average_monthly_albedo.tif")%>%
  terra::project(dhsvm, method = "bilinear", mask = TRUE)
albedo
plot(albedo)
hist(values(albedo))


veg_categories <- unique(values(dhsvm)[!is.na(values(dhsvm))])
veg_table <- data.frame(veg_cat = numeric(0),
                        albedos = numeric(0))
for(i in 1:length(veg_categories)){

veg_categories[i]
albedos <- unlist(lapply(albedo, FUN = function(x) {mean(values(x)[values(dhsvm) == veg_categories[i]], na.rm = TRUE)}))
veg_table <- rbind(veg_table, cbind(veg_categories[i], albedos))

}

veg_table$month = rep(1:12, times = length(veg_categories))
veg_table$albedos <- veg_table$albedos/1000
veg_table <- veg_table[order(veg_table$V1), ]

write.csv(veg_table, "albedo_by_month_with_snow.csv")



