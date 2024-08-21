#LAI example figure

lai <- function(biomass, maxLAI, KLAI){
  return(maxLAI*(biomass/(biomass+KLAI)))
}

maxLAI = 7
KLAI = 12000
minLAI = 0.1
biomass <- seq(0, 20000, length.out = 200)


lai_est <- lai(biomass, maxLAI, KLAI)
lai_est <- ifelse(lai_est < minLAI, minLAI, lai_est)

plot(lai_est ~ biomass, 
     type = "l",
     lwd = 2,
     ylim = c(0, 10),
     xlab = "Cohort biomass",
     ylab = "Cohort LAI")
abline(h = maxLAI/2)
abline(v = KLAI)


maxLAI = 5
KLAI = 3500
minLAI = 0.1
biomass <- seq(0, 20000, length.out = 200)

lai_est <- lai(biomass, maxLAI, KLAI)
lai_est <- ifelse(lai_est < minLAI, minLAI, lai_est)
lines(lai_est ~ biomass, type = "l", lwd = 2)

maxLAI = 3
KLAI = 2000
minLAI = 0.1
biomass <- seq(0, 20000, length.out = 200)

lai_est <- lai(biomass, maxLAI, KLAI)
lai_est <- ifelse(lai_est < minLAI, minLAI, lai_est)
lines(lai_est ~ biomass, type = "l")





