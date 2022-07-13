# Classification of stands based on species abundance
# Based on work by Charles Maxwell

library(raster)
library(rasterVis)
library(RColorBrewer)

mycolors = c(brewer.pal(name="Dark2", n = 8), brewer.pal(name="Paired", n = 10))

TCSI_template <- raster("./Models/Inputs/masks_boundaries/mask.tif")
TCSI_template
plot(TCSI_template)

####sierraN###################################

tcsi_dir <- "E:/TCSI LANDIS/Scenario1 - historical - Run 1/biomass/"
TCSI_template <- raster("./Models/Inputs/masks_boundaries/mask.tif")

rac <- raster(paste0(tcsi_dir, "AbieConc-",i,".img"))
ram <- raster(paste0(tcsi_dir, "AbieMagn-",i,".img"))
rcd <- raster(paste0(tcsi_dir, "CaloDecu-",i,".img"))
rjo <- raster(paste0(tcsi_dir, "JuniOcci-",i,".img"))
rpa <- raster(paste0(tcsi_dir, "PinuAlbi-",i,".img"))
# rpb <- raster(paste0(tcsi_dir, "PinuBalf-",i,".img"))
rpc <- raster(paste0(tcsi_dir, "PinuCont-",i,".img"))
rpj <- raster(paste0(tcsi_dir, "PinuJeff-",i,".img"))
rpl <- raster(paste0(tcsi_dir, "PinuLamb-",i,".img"))
rpm <- raster(paste0(tcsi_dir, "PinuMont-",i,".img"))
rpp <- raster(paste0(tcsi_dir, "PinuPond-",i,".img"))
rpw <- raster(paste0(tcsi_dir, "PinuWash-",i,".img"))
rpz <- raster(paste0(tcsi_dir, "PseuMenz-",i,".img"))
# rsg <- raster(paste0(tcsi_dir, "SequGiga-",i,".img"))
rtc <- raster(paste0(tcsi_dir, "TorrCali-",i,".img"))
rtm <- raster(paste0(tcsi_dir, "TsugMert-",i,".img"))

racer <- raster(paste0(tcsi_dir, "AcerMacr-",i,".img"))
ralnu <- raster(paste0(tcsi_dir, "AlnuRhom-",i,".img"))
# ralnr <- raster(paste0(tcsi_dir, "AlnuRubr-",i,".img"))
rarbu <- raster(paste0(tcsi_dir, "ArbuMenz-",i,".img"))
rlith <- raster(paste0(tcsi_dir, "LithDens-",i,".img"))
rpopu <- raster(paste0(tcsi_dir, "PopuTrem-",i,".img"))
rquec <- raster(paste0(tcsi_dir, "QuerChry-",i,".img"))
rquek <- raster(paste0(tcsi_dir, "QuerKell-",i,".img"))
rqudo <- raster(paste0(tcsi_dir, "QuerDoug-",i,".img"))
rquew <- raster(paste0(tcsi_dir, "QuerWisl-",i,".img"))
rumbe <- raster(paste0(tcsi_dir, "UmbeCali-",i,".img"))

rfrs <- raster(paste0(tcsi_dir, "FX_R_SEED-",i,".img"))
rnrs <- raster(paste0(tcsi_dir, "NOFX_R_SEED-",i,".img"))
rnns <- raster(paste0(tcsi_dir, "NOFX_NOR_SEED-",i,".img"))

total.biomass <- raster(paste0(tcsi_dir, "TotalBiomass-",i,".img"))

biomass.SMC = rac + rcd + rpc + rpj + rpl + rpp + rpw + rpz + rtc # + rpb
biomass.SCN = ram + rpa + rpm + rtm

biom.coni <- rac + ram + rcd + rjo + rpa +  rpc + rpj + rpl + 
  rpm + rpp + rpw + rpz  + rtc + rtm #rpb + rsg
biom.hard <- racer + ralnu  + rarbu + rlith  + rumbe + rpopu + rquec + rquek + rquew + rqudo #+ ralnr
biom.shrub <- rfrs + rnrs + rnns

cover.coni <- biom.coni > biom.hard & biom.coni > biom.shrub
cover.hard <- biom.hard > biom.coni & (biom.coni/total.biomass) < 0.25
cover.mixed <- biom.hard > biom.coni & (biom.coni/total.biomass) > 0.25
cover.shrub <- biom.shrub > biom.hard & biom.shrub > biom.coni & cover.mixed == FALSE & cover.hard == FALSE & cover.coni == FALSE

plot(cover.coni)
plot(cover.hard)
plot(cover.mixed)
plot(cover.shrub)

mhc <- cover.mixed
asp <- cover.hard == 1 & (rpopu / biom.hard) > 0.5
mhw <- cover.hard == 1 & ((rarbu + rlith + rquec + rquek + rqudo + rquew + rumbe) / biom.hard) > 0.5 #this one
# mri <- cover.hard == 1 & ((racer + ralnu + ralnr) / biom.hard) > 0.5
mri <- cover.hard == 1 & ((racer + ralnu) / biom.hard) > 0.5
wfr <- cover.coni == 1 & (rac/biom.coni) > 0.5
rfr <- cover.coni == 1 & (ram/biom.coni) > 0.5
jpn <- cover.coni == 1 & (rpj /biom.coni) > 0.5
dfr <- cover.coni == 1 & (rpz/biom.coni) >0.5
ppn <- cover.coni == 1 & ((rpp + rpw)/biom.coni) > 0.5
lpn <- cover.coni == 1 & (rpc/biom.coni) > 0.5
juo <- cover.coni == 1 & (rjo/biom.coni) > 0.5
# seg <- cover.coni == 1 & (rsg/biom.coni) > 0.5
mch <- cover.shrub == 1 #this one

smc <- cover.coni == 1 & ((rac/biom.coni) < 0.5) & ((ram/biom.coni) < 0.5) & ((rcd/biom.coni) < 0.5) & ((rjo /biom.coni) < 0.5) &
  ((rpa/biom.coni) < 0.5) & ((rpc/biom.coni) < 0.5) & ((rpj/biom.coni) < 0.5) & ((rpl/biom.coni) < 0.5) &
  (((rpp + rpw)/biom.coni) < 0.5) & ((rpz/biom.coni) < 0.5)  &
  ((rtc/biom.coni) < 0.5) & ((rtm/biom.coni) < 0.5) & biomass.SMC > biomass.SCN # & ((rpb/biom.coni) < 0.5) & ((rsg/biom.coni) < 0.5)

scn <- cover.coni == 1 & ((rac/biom.coni) < 0.5) & ((ram/biom.coni) < 0.5) & ((rcd/biom.coni) < 0.5) & ((rjo /biom.coni) < 0.5) &
  ((rpa/biom.coni) < 0.5)  & ((rpc/biom.coni) < 0.5) & ((rpj/biom.coni) < 0.5) & ((rpl/biom.coni) < 0.5) &
  (((rpp + rpw)/biom.coni) < 0.5) & ((rpz/biom.coni) < 0.5)  &
  ((rtc/biom.coni) < 0.5) & ((rtm/biom.coni) < 0.5) & biomass.SMC < biomass.SCN #& ((rpb/biom.coni) < 0.5) & ((rsg/biom.coni) < 0.5)

yps <- cover.coni == 1 & ((rpj + rpp + rpw)/biom.coni) > 0.5
wps <- cover.coni == 1 & ((rpl + rpm + rpa)/biom.coni) > 0.5

mhc[mhc[] == 1] <- 1
mhw[mhw[] == 1] <- 2
mri[mri[] == 1] <- 3
wfr[wfr[] == 1] <- 4
dfr[dfr[] == 1] <- 5
rfr[rfr[] == 1] <- 6
ppn[ppn[] == 1] <- 7
jpn[jpn[] == 1] <- 8
lpn[lpn[] == 1] <- 9
# seg[seg[] == 1] <- 10

asp[asp[] == 1] <- 12
juo[juo[] == 1] <- 13
smc[smc[] == 1] <- 14
scn[scn[] == 1] <- 15


mch[mch[] == 1] <- 18



r.fortype <- mhc + mhw + mri + wfr + dfr + rfr + ppn + jpn + lpn  + smc + scn + asp + juo + mch #+ seg
plot(r.fortype)
freq(r.fortype)
 
# error_cell <- r.fortype
# values(error_cell) <- 0
# error_cell[r.fortype[] == 2] <- 1
# 
# mch[which(error_cell[] == 1)] #mch and mhw are the problem

r.fortype@crs <- TCSI_template@crs
r.fortype@extent <- TCSI_template@extent

writeRaster(r.fortype, paste0("./Analysis/DHSVM/Forest_type", i, ".tif"), overwrite = T)

####Combine maps####

CA_cwhr2 <- as.factor(raster(paste0("./Analysis/DHSVM/Forest_type", i, ".tif")))
tar <- levels(CA_cwhr2)[[1]]
tar[["cwhr_class"]] <- c("inactive", "Mixed hardwood conifer",	"Mixed hardwood",	"Mixed riparian",	"white fir",	"doug-fir",	
                         "red fir",	"ponderosa pine",	"Jeffrey pine",	"lodgepole pine",	"aspen",
                         "western juniper",	"Sierra mixed conifer",	"High elevation Sierra mixed conifer", "shrub")

levels(CA_cwhr2) <- tar

levelplot(CA_cwhr2, col.regions=mycolors)