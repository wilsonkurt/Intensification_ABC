#Add code to ignore if the model PrDens goes to NA at any point (which happens if we drop to a population of 0)

createmodelspd <- function(model) {
  
  #if the simulated model had a population crash to 0 (and therefore PrDens is NA after the pop = 0 point) record it as an NA
  if (any(is.na(model$PrDens)|is.nan(model$PrDens)|model$PrDens==Inf))
  {
    return(data.frame(
      calBP = seq(5000,0,-1),
    #  ucs.norm.calBP = seq(5000,0,-1),
      ucs.norm.prdens = NA,
     # ucs.nnorm.calBP = seq(5000,0,-1),
      ucs.nnorm.prdens = NA,
      #cs.norm.calBP = seq(5000,0,-1),
      cs.norm.prdens = NA,
      #cs.nnorm.calBP = seq(5000,0,-1),
      cs.nnorm.prdens = NA,
      param.combo = param.combo
    ))
  }  

#Now create spds from our model simulated data
d <- uncalibrate(model, calCurves = 'intcal20', verbose = F)#uncalibrate the dates from the simulated growth model
n.samples <- length(index)#we will sample the uncalibtrated dates equal to the unique bins
errors <- sample(dat_cal.Plat.norm$metadata$Error, size = n.samples, replace =TRUE) #generate date errors by sampling from the observed errors of the actual dates
uncal.samples <- sample(d$CRA, size = n.samples, prob = d$PrDens, replace = TRUE)#sample the n.samples (so # of dates) we had using that calcurve from the uncalibrated growth model output - weight the probability of selection by the PrDens of the growth model output
cal.samples <- sample(d$CRA, size = n.samples, prob = d$Raw, replace = TRUE)#do the same thing, but instead of weighting hte probability by the PrDens, weight by the Raw output from the uncalibrate process
ids <- idCounter + (1:n.samples)#assign an ID (this will just be a unique number from 1 to the total # samples we generated) 
uncalsample.norm <- calibrate(uncal.samples, errors, calCurves = 'intcal20', ids = ids, normalised = TRUE, verbose = FALSE) #calibrate a new SPD using the artificial 14c dates sampled from the growthModel output based on sampling using PrDens while  normalising the SPD
uncalsample.nnorm <- calibrate(uncal.samples, errors, calCurves = 'intcal20', ids = ids, normalised = FALSE, verbose = FALSE) #calibrate a new SPD using the artificial 14c dates sampled from the growthModel output based on sampling using Raw while not normalising the SPD
calsample.norm <- calibrate(cal.samples, errors, calCurves = 'intcal20', ids = ids, normalised = TRUE, verbose = FALSE)#calibrate a new SPD using the artificial 14c dates sampled from the growthModel output based on sampling using Raw while normalising the SPD
calsample.nnorm <- calibrate(cal.samples, errors, calCurves = 'intcal20', ids = ids, normalised = FALSE, verbose = FALSE) #calibrate a new SPD using the artificial 14c dates sampled from the growthModel output based on sampling using Raw while not normalising the SPD

#Now make an spd using the resampled dates for each of normalised, not normalised, and each version of back calculating 14c dates
spd.uncalsample.norm <- spd(uncalsample.norm, timeRange = c(5000,0), spdnormalised = TRUE, verbose = FALSE)
spd.uncalsample.norm$metadata$combo <- i #add the param combo signifier to the metadata
spd.uncalsample.nnorm <- spd(uncalsample.nnorm, timeRange = c(5000,0), spdnormalised = TRUE, verbose = FALSE)
spd.uncalsample.nnorm$metadata$combo <- i
spd.calsample.norm <- spd(calsample.norm, timeRange = c(5000,0), spdnormalised = TRUE, verbose = FALSE)
spd.calsample.norm$metadata$combo <- i
spd.calsample.nnorm <- spd(calsample.nnorm, timeRange = c(5000,0), spdnormalised = TRUE, verbose = FALSE)
spd.calsample.nnorm$metadata$combo <- i

return(data.frame(
  calBP = seq(5000,0,-1),
  #ucs.norm.calBP = spd.uncalsample.norm$grid$calBP,
  ucs.norm.prdens = spd.uncalsample.norm$grid$PrDens,
  #ucs.nnorm.calBP = spd.uncalsample.nnorm$grid$calBP,
  ucs.nnorm.prdens = spd.uncalsample.nnorm$grid$PrDens,
  #cs.norm.calBP = spd.calsample.norm$grid$calBP,
  cs.norm.prdens = spd.calsample.norm$grid$PrDens,
  #cs.nnorm.calBP = spd.calsample.nnorm$grid$calBP,
  cs.nnorm.prdens = spd.calsample.nnorm$grid$PrDens,
  param.combo = param.combo
))

}