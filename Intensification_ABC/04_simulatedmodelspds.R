#Add code to ignore if the model PrDens goes to NA at any point (which happens if we drop to a population of 0)

createmodelspd <- function(model) {
  
  #if the simulated model throws an NA in the population size, record the evaluation values as NAs as well
  if (any(is.na(model$PrDens)|is.nan(model$PrDens)|model$PrDens==Inf))
    
  {
    
    mod.results.df <- data.frame(param.combo = param.combo,
                                 euc.uncal.norm = NA,
                                 euc.uncal.nnorm = NA,
                                 euc.cal.norm = NA,
                                 euc.cal.nnorm = NA,
                                 nrmse.uncal.norm = NA,
                                 nrmse.uncal.nnorm = NA,
                                 nrmse.cal.norm = NA,
                                 nrmse.cal.nnorm = NA)
    mod.results.df[,14:16017] <- NA
    
    return(mod.results.df) 
    
  }  
  
  #Now create spds from our model simulated data
  d <- uncalibrate(model, calCurves = 'intcal20', verbose = F)#uncalibrate the spd from the simulated growth model
  n.samples <- bins.use #line 34 is Dinappoli code - for use in CHPC we just use the set # not the length of allbins object
  #n.samples <- length(allbins)#we will sample the uncalibtrated dates equal to the unique bins
  
 # errors <- sample(dat_cal.Plat.norm$metadata$Error, size = n.samples, replace =TRUE) #generate date errors by sampling 
  errors <- sample(region.errors, size = n.samples, replace =TRUE) #generate date errors by sampling 
  
  #from the observed errors of the actual dates
  uncal.samples <- sample(d$CRA, size = n.samples, prob = d$PrDens, replace = TRUE)#sample the n.samples (so # of dates) 
  #we had using the calcurve from the uncalibrated growth model output - weight the probability of selection by the 
  #PrDens of the growth model output. d$CRA is sampling a year to serve as the year of the date
  cal.samples <- sample(d$CRA, size = n.samples, prob = d$Raw, replace = TRUE)#do the same thing, but instead of weighting 
  #hte probability by the PrDens, weight by the Raw output from the uncalibrate process
  ids <- idCounter + (1:n.samples)#assign an ID (this will just be a unique number from 1 to the total # samples we generated) 
  uncalsample.norm <- calibrate(uncal.samples, errors, calCurves = 'intcal20', ids = ids, normalised = TRUE, verbose = FALSE) #calibrate 
  #a new SPD using the artificial 14c dates sampled from the growthModel output based on sampling using PrDens while  normalising the SPD
  uncalsample.nnorm <- calibrate(uncal.samples, errors, calCurves = 'intcal20', ids = ids, normalised = FALSE, verbose = FALSE) #calibrate 
  #a new SPD using the artificial 14c dates sampled from the growthModel output based on sampling using Raw while not normalising the SPD
  calsample.norm <- calibrate(cal.samples, errors, calCurves = 'intcal20', ids = ids, normalised = TRUE, verbose = FALSE)#calibrate a new 
  #SPD using the artificial 14c dates sampled from the growthModel output based on sampling using Raw while normalising the SPD
  calsample.nnorm <- calibrate(cal.samples, errors, calCurves = 'intcal20', ids = ids, normalised = FALSE, verbose = FALSE) #calibrate a 
  #new SPD using the artificial 14c dates sampled from the growthModel output based on sampling using Raw while not normalising the SPD
  
  #Now make an spd using the resampled dates for each of normalised, not normalised, and each version of back calculating 14c dates
  spd.uncalsample.norm <- spd(uncalsample.norm, timeRange = c(4000,0), spdnormalised = TRUE, verbose = FALSE)
  spd.uncalsample.norm$metadata$combo <- i #add the param combo signifier to the metadata
  spd.uncalsample.nnorm <- spd(uncalsample.nnorm, timeRange = c(4000,0), spdnormalised = TRUE, verbose = FALSE)
  spd.uncalsample.nnorm$metadata$combo <- i
  spd.calsample.norm <- spd(calsample.norm, timeRange = c(4000,0), spdnormalised = TRUE, verbose = FALSE)
  spd.calsample.norm$metadata$combo <- i
  spd.calsample.nnorm <- spd(calsample.nnorm, timeRange = c(4000,0), spdnormalised = TRUE, verbose = FALSE)
  spd.calsample.nnorm$metadata$combo <- i
  
  #observed.spds.norm
  #observed.spds.nnorm
  
  #Compute the average euclidean distance measure for simulated spd versus all 1000 observed spd possibilities. 
  #To do so subtract the simulated spd from each of the 1000 observed spds to get the distance at each year for 
  #each of the 1000 possibilites. Then square the sum those distances (get the cumulative distance for each of 
  #the 1000 observations). Take the square root of that value.
  euc.uncal.norm <- data.frame()
  for (i in 2:length(observed.spds.norm)) {
    euc <- sqrt(sum((observed.spds.norm[yrs.match.range,i] - spd.uncalsample.norm$grid$PrDens[yrs.match.range])^2))
    euc.uncal.norm <- rbind(euc.uncal.norm, euc)
    
  }
  euc.uncal.norm <- colMeans(euc.uncal.norm)  
  
  euc.uncal.nnorm <- data.frame()
  for (i in 2:length(observed.spds.nnorm)) {
    euc <- sqrt(sum((observed.spds.nnorm[yrs.match.range,i] - spd.uncalsample.nnorm$grid$PrDens[yrs.match.range])^2))
    euc.uncal.nnorm <- rbind(euc.uncal.nnorm, euc)
    
  }
  euc.uncal.nnorm <- colMeans(euc.uncal.nnorm)  
  
  euc.cal.norm <- data.frame()
  for (i in 2:length(observed.spds.norm)) {
    euc <- sqrt(sum((observed.spds.norm[yrs.match.range,i] - spd.calsample.norm$grid$PrDens[yrs.match.range])^2))
    euc.cal.norm <- rbind(euc.cal.norm, euc)
    
  }
  euc.cal.norm <- colMeans(euc.cal.norm)  
  
  euc.cal.nnorm <- data.frame()
  for (i in 2:length(observed.spds.nnorm)) {
    euc <- sqrt(sum((observed.spds.nnorm[yrs.match.range,i] - spd.calsample.nnorm$grid$PrDens[yrs.match.range])^2))
    euc.cal.nnorm <- rbind(euc.cal.nnorm, euc)
    
  }
  euc.cal.nnorm <- colMeans(euc.cal.nnorm) 
  
  
  #Compute the average euclidean distance adapted with delta for simulated spd vs all 1000 observed spd possibilities.
  #We are performing matching only on the simulated years 4000 to 1000 yBP
  euc.delta.uncal.norm <- data.frame()
  for (i in 2:length(observed.spds.norm)) {
    euc.delta <- sqrt(sum(((observed.spds.norm[yrs.match.range,i] - spd.uncalsample.norm$grid$PrDens[yrs.match.range]) / 
                         (spd.uncalsample.norm$grid$PrDens[yrs.match.range] + delta.val))^2))
    euc.delta.uncal.norm <- rbind(euc.delta.uncal.norm, euc.delta)
  }
  euc.delta.uncal.norm <- colMeans(euc.delta.uncal.norm)
  
  euc.delta.uncal.nnorm <- data.frame()
  for (i in 2:length(observed.spds.nnorm)) {
    euc.delta <- sqrt(sum(((observed.spds.nnorm[yrs.match.range,i] - spd.uncalsample.nnorm$grid$PrDens[yrs.match.range]) / 
                             (spd.uncalsample.nnorm$grid$PrDens[yrs.match.range] + delta.val))^2))
    euc.delta.uncal.nnorm <- rbind(euc.delta.uncal.nnorm, euc.delta)
  }
  euc.delta.uncal.nnorm <- colMeans(euc.delta.uncal.nnorm)
  
  euc.delta.cal.norm <- data.frame()
  for (i in 2:length(observed.spds.norm)) {
    euc.delta <- sqrt(sum(((observed.spds.norm[yrs.match.range,i] - spd.calsample.norm$grid$PrDens[yrs.match.range]) / 
                             (spd.calsample.norm$grid$PrDens[yrs.match.range] + delta.val))^2))
    euc.delta.cal.norm <- rbind(euc.delta.cal.norm, euc.delta)
  }
  euc.delta.cal.norm <- colMeans(euc.delta.cal.norm)
  
  euc.delta.cal.nnorm <- data.frame()
  for (i in 2:length(observed.spds.nnorm)) {
    euc.delta <- sqrt(sum(((observed.spds.nnorm[yrs.match.range,i] - spd.calsample.nnorm$grid$PrDens[yrs.match.range]) / 
                             (spd.calsample.nnorm$grid$PrDens[yrs.match.range] + delta.val))^2))
    euc.delta.cal.nnorm <- rbind(euc.delta.cal.nnorm, euc.delta)
  }
  euc.delta.cal.nnorm <- colMeans(euc.delta.cal.nnorm)
 
  
  #Compute the average euclidean distance  weighted towards key years simulated spd vs all 1000 observed spd possibilities.
  #We are performing matching only on the simulated years 4000 to 1000 yBP
  euc.weighted.uncal.norm <- data.frame()
  for (i in 2:length(observed.spds.norm)) {
   euc.weighted <-  sqrt((1/sum(weights.vec) *
                       (sum(
                         weights.vec[i] * (observed.spds.norm[yrs.match.range,i] - spd.uncalsample.norm$grid$PrDens[yrs.match.range])^2
                        ))
                      )
                    )
   euc.weighted.uncal.norm <- rbind(euc.weighted.uncal.norm, euc.weighted)
  }
  euc.weighted.uncal.norm <- colMeans(euc.weighted.uncal.norm)
  
  euc.weighted.uncal.nnorm <- data.frame()
  for (i in 2:length(observed.spds.nnorm)) {
    euc.weighted <-  sqrt((1/sum(weights.vec) *
                             (sum(
                               weights.vec[i] * (observed.spds.nnorm[yrs.match.range,i] - spd.uncalsample.nnorm$grid$PrDens[yrs.match.range])^2
                             ))
    )
    )
    euc.weighted.uncal.nnorm <- rbind(euc.weighted.uncal.nnorm, euc.weighted)
  }
  euc.weighted.uncal.nnorm <- colMeans(euc.weighted.uncal.nnorm)
  
  euc.weighted.cal.norm <- data.frame()
  for (i in 2:length(observed.spds.norm)) {
    euc.weighted <-  sqrt((1/sum(weights.vec) *
                             (sum(
                               weights.vec[i] * (observed.spds.norm[yrs.match.range,i] - spd.calsample.norm$grid$PrDens[yrs.match.range])^2
                             ))
    )
    )
    euc.weighted.cal.norm <- rbind(euc.weighted.cal.norm, euc.weighted)
  }
  euc.weighted.cal.norm <- colMeans(euc.weighted.cal.norm)
  
  euc.weighted.cal.nnorm <- data.frame()
  for (i in 2:length(observed.spds.nnorm)) {
    euc.weighted <-  sqrt((1/sum(weights.vec) *
                             (sum(
                               weights.vec[i] * (observed.spds.nnorm[yrs.match.range,i] - spd.calsample.nnorm$grid$PrDens[yrs.match.range])^2
                             ))
    )
    )
    euc.weighted.cal.nnorm <- rbind(euc.weighted.cal.nnorm, euc.weighted)
  }
  euc.weighted.cal.nnorm <- colMeans(euc.weighted.cal.nnorm)
  
  
  mod.results.df <- data.frame(param.combo = param.combo,
    euc.uncal.norm = euc.uncal.norm,
    euc.uncal.nnorm = euc.uncal.nnorm,
    euc.cal.norm = euc.cal.norm,
    euc.cal.nnorm = euc.cal.nnorm,
    euc.delta.uncal.norm = euc.delta.uncal.norm,
    euc.delta.uncal.nnorm = euc.delta.uncal.nnorm,
    euc.delta.cal.norm = euc.delta.cal.norm,
    euc.delta.cal.nnorm = euc.delta.cal.nnorm,
    euc.weighted.uncal.norm = euc.weighted.uncal.norm,
    euc.weighted.uncal.nnorm = euc.weighted.uncal.nnorm,
    euc.weighted.cal.norm = euc.weighted.cal.norm,
    euc.weighted.cal.nnorm = euc.weighted.cal.nnorm)
  mod.results.df[,14:4014] <- spd.uncalsample.norm$grid$PrDens
  mod.results.df[,4015:8015] <- spd.uncalsample.nnorm$grid$PrDens
  mod.results.df[,8016:12016] <- spd.calsample.norm$grid$PrDens
  mod.results.df[,12017:16017] <- spd.calsample.nnorm$grid$PrDens
 
  return(mod.results.df) 
  
  #return a dataframe with the evaluation metrics included
#  return(data.frame(
#    param.combo = param.combo,
#    euc.uncal.norm = euc.uncal.norm,
#    euc.uncal.nnorm = euc.uncal.nnorm,
#    euc.cal.norm = euc.cal.norm,
#    euc.cal.nnorm = euc.cal.nnorm,
#    nrmse.uncal.norm = nrmse.uncal.norm,
#    nrmse.uncal.nnorm = nrmse.uncal.nnorm,
#    nrmse.cal.norm = nrmse.cal.norm,
#    nrmse.cal.nnorm = nrmse.cal.nnorm
    
    #code below would make an output containing all the spds, not just the distance evaluation metrics
    #calBP = seq(4000,0,-1),
    #ucs.norm.calBP = spd.uncalsample.norm$grid$calBP,
    #ucs.norm.prdens = spd.uncalsample.norm$grid$PrDens,
    #ucs.nnorm.calBP = spd.uncalsample.nnorm$grid$calBP,
    #ucs.nnorm.prdens = spd.uncalsample.nnorm$grid$PrDens,
    #cs.norm.calBP = spd.calsample.norm$grid$calBP,
    #cs.norm.prdens = spd.calsample.norm$grid$PrDens,
    #cs.nnorm.calBP = spd.calsample.nnorm$grid$calBP,
    #cs.nnorm.prdens = spd.calsample.nnorm$grid$PrDens,
  ##))
  
}

####CODE BELOW WOULD ENABLE COMPUTATION OF NRMSE FOR PATTERN MATCHING
# #Compute the average normalised root mean square error for simulated spd vs all 1000 observed spd possibilities.
# #We are performing matching only on the simulated years 4000 to 1000 yBP
# nrmse.uncal.norm <- data.frame()
# for (i in 2:length(observed.spds.norm)) {
#   nrmse <- sqrt(sum((observed.spds.norm[yrs.match.range,i] - spd.uncalsample.norm$grid$PrDens[yrs.match.range])^2) / length(model$CalBP[yrs.match.range])) / sd(observed.spds.norm[yrs.match.range,i])
#   nrmse.uncal.norm <- rbind(nrmse.uncal.norm, nrmse)
# }
# nrmse.uncal.norm <- colMeans(nrmse.uncal.norm)
# 
# nrmse.uncal.nnorm <- data.frame()
# for (i in 2:length(observed.spds.nnorm)) {
#   nrmse <- sqrt(sum((observed.spds.nnorm[yrs.match.range,i] - spd.uncalsample.nnorm$grid$PrDens[yrs.match.range])^2) / length(model$CalBP[yrs.match.range])) / sd(observed.spds.nnorm[yrs.match.range,i])
#   nrmse.uncal.nnorm <- rbind(nrmse.uncal.nnorm, nrmse)
# }
# nrmse.uncal.nnorm <- colMeans(nrmse.uncal.nnorm)
# 
# nrmse.cal.norm <- data.frame()
# for (i in 2:length(observed.spds.norm)) {
#   nrmse <- sqrt(sum((observed.spds.norm[yrs.match.range,i] - spd.calsample.norm$grid$PrDens[yrs.match.range])^2) / length(model$CalBP[yrs.match.range])) / sd(observed.spds.norm[yrs.match.range,i])
#   nrmse.cal.norm <- rbind(nrmse.cal.norm, nrmse)
# }
# nrmse.cal.norm <- colMeans(nrmse.cal.norm)
# 
# nrmse.cal.nnorm <- data.frame()
# for (i in 2:length(observed.spds.nnorm)) {
#   nrmse <- sqrt(sum((observed.spds.nnorm[yrs.match.range,i] - spd.calsample.nnorm$grid$PrDens[yrs.match.range])^2) / length(model$CalBP[yrs.match.range])) / sd(observed.spds.nnorm[yrs.match.range,i])
#   nrmse.cal.nnorm <- rbind(nrmse.cal.nnorm, nrmse)
# }
# nrmse.cal.nnorm <- colMeans(nrmse.cal.nnorm)
