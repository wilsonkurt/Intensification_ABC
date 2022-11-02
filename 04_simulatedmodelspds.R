#Add code to ignore if the model PrDens goes to NA at any point (which happens if we drop to a population of 0)

createmodelspd <- function(model) {
  
  #if the simulated model had a population crash to 0 (and therefore PrDens is NA after the pop = 0 point) record it as an NA
  if (any(is.na(model$PrDens)|is.nan(model$PrDens)|model$PrDens==Inf))
  {
    return(data.frame(
      #   calBP = seq(5000,0,-1),
      # #  ucs.norm.calBP = seq(5000,0,-1),
      #   ucs.norm.prdens = NA,
      #  # ucs.nnorm.calBP = seq(5000,0,-1),
      #   ucs.nnorm.prdens = NA,
      #   #cs.norm.calBP = seq(5000,0,-1),
      #   cs.norm.prdens = NA,
      #   #cs.nnorm.calBP = seq(5000,0,-1),
      #   cs.nnorm.prdens = NA,
      param.combo = param.combo,
      euc.uncal.norm = NA,
      euc.uncal.nnorm = NA,
      euc.cal.norm = NA,
      euc.cal.nnorm = NA,
      nrmse.uncal.norm = NA,
      nrmse.uncal.nnorm = NA,
      nrmse.cal.norm = NA,
      nrmse.cal.nnorm = NA
    ))
  }  
  
  #Now create spds from our model simulated data
  d <- uncalibrate(model, calCurves = 'intcal20', verbose = F)#uncalibrate the spd from the simulated growth model
  n.samples <- length(allbins)#we will sample the uncalibtrated dates equal to the unique bins
  errors <- sample(dat_cal.Plat.norm$metadata$Error, size = n.samples, replace =TRUE) #generate date errors by sampling 
  #from the observed errors of the actual dates
  uncal.samples <- sample(d$CRA, size = n.samples, prob = d$PrDens, replace = TRUE)#sample the n.samples (so # of dates) 
  #we had using that calcurve from the uncalibrated growth model output - weight the probability of selection by the 
  #PrDens of the growth model output
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
  spd.uncalsample.norm <- spd(uncalsample.norm, timeRange = c(5000,0), spdnormalised = TRUE, verbose = FALSE)
  spd.uncalsample.norm$metadata$combo <- i #add the param combo signifier to the metadata
  spd.uncalsample.nnorm <- spd(uncalsample.nnorm, timeRange = c(5000,0), spdnormalised = TRUE, verbose = FALSE)
  spd.uncalsample.nnorm$metadata$combo <- i
  spd.calsample.norm <- spd(calsample.norm, timeRange = c(5000,0), spdnormalised = TRUE, verbose = FALSE)
  spd.calsample.norm$metadata$combo <- i
  spd.calsample.nnorm <- spd(calsample.nnorm, timeRange = c(5000,0), spdnormalised = TRUE, verbose = FALSE)
  spd.calsample.nnorm$metadata$combo <- i
  
  #observed.spds.norm
  #observed.spds.nnorm
  
  #Compute the average euclidean distance measure for simulated spd versus all 1000 observed spd possibilities. 
  #To do so subtract the simulated spd from each of the 1000 observed spds to get the distance at each year for 
  #each of the 1000 possibilites. Then square the sum those distances (get the cumulative distance for each of 
  #the 1000 observations). Take the square root of that value.
  euc.uncal.norm <- data.frame()
  for (i in 2:length(observed.spds.norm)) {
    euc <- sqrt(sum((observed.spds.norm[,i] - spd.uncalsample.norm$grid$PrDens)^2))
    euc.uncal.norm <- rbind(euc.uncal.norm, euc)
    
  }
  euc.uncal.norm <- colMeans(euc.uncal.norm)  
  
  euc.uncal.nnorm <- data.frame()
  for (i in 2:length(observed.spds.nnorm)) {
    euc <- sqrt(sum((observed.spds.nnorm[,i] - spd.uncalsample.nnorm$grid$PrDens)^2))
    euc.uncal.nnorm <- rbind(euc.uncal.nnorm, euc)
    
  }
  euc.uncal.nnorm <- colMeans(euc.uncal.nnorm)  
  
  euc.cal.norm <- data.frame()
  for (i in 2:length(observed.spds.norm)) {
    euc <- sqrt(sum((observed.spds.norm[,i] - spd.calsample.norm$grid$PrDens)^2))
    euc.cal.norm <- rbind(euc.cal.norm, euc)
    
  }
  euc.cal.norm <- colMeans(euc.cal.norm)  
  
  euc.cal.nnorm <- data.frame()
  for (i in 2:length(observed.spds.nnorm)) {
    euc <- sqrt(sum((observed.spds.nnorm[,i] - spd.calsample.nnorm$grid$PrDens)^2))
    euc.cal.nnorm <- rbind(euc.cal.nnorm, euc)
    
  }
  euc.cal.nnorm <- colMeans(euc.cal.nnorm) 
  
  
  #Compute the average normalised root mean square error for simulated spd vs all 1000 observed spd possibilities.
  nrmse.uncal.norm <- data.frame()
  for (i in 2:length(observed.spds.norm)) {
    nrmse <- sqrt(sum((observed.spds.norm[,i] - spd.uncalsample.norm$grid$PrDens)^2) / length(model$CalBP)) / sd(observed.spds.norm[,i])
    nrmse.uncal.norm <- rbind(nrmse.uncal.norm, nrmse)
  }
  nrmse.uncal.norm <- colMeans(nrmse.uncal.norm)
  
  nrmse.uncal.nnorm <- data.frame()
  for (i in 2:length(observed.spds.nnorm)) {
    nrmse <- sqrt(sum((observed.spds.nnorm[,i] - spd.uncalsample.nnorm$grid$PrDens)^2) / length(model$CalBP)) / sd(observed.spds.nnorm[,i])
    nrmse.uncal.nnorm <- rbind(nrmse.uncal.nnorm, nrmse)
  }
  nrmse.uncal.nnorm <- colMeans(nrmse.uncal.nnorm)
  
  nrmse.cal.norm <- data.frame()
  for (i in 2:length(observed.spds.norm)) {
    nrmse <- sqrt(sum((observed.spds.norm[,i] - spd.calsample.norm$grid$PrDens)^2) / length(model$CalBP)) / sd(observed.spds.norm[,i])
    nrmse.cal.norm <- rbind(nrmse.cal.norm, nrmse)
  }
  nrmse.cal.norm <- colMeans(nrmse.cal.norm)
  
  nrmse.cal.nnorm <- data.frame()
  for (i in 2:length(observed.spds.nnorm)) {
    nrmse <- sqrt(sum((observed.spds.nnorm[,i] - spd.calsample.nnorm$grid$PrDens)^2) / length(model$CalBP)) / sd(observed.spds.nnorm[,i])
    nrmse.cal.nnorm <- rbind(nrmse.cal.nnorm, nrmse)
  }
  nrmse.cal.nnorm <- colMeans(nrmse.cal.nnorm)
  
  #return a dataframe with the evaluation metrics included
  return(data.frame(
    param.combo = param.combo,
    euc.uncal.norm = euc.uncal.norm,
    euc.uncal.nnorm = euc.uncal.nnorm,
    euc.cal.norm = euc.cal.norm,
    euc.cal.nnorm = euc.cal.nnorm,
    nrmse.uncal.norm = nrmse.uncal.norm,
    nrmse.uncal.nnorm = nrmse.uncal.nnorm,
    nrmse.cal.norm = nrmse.cal.norm,
    nrmse.cal.nnorm = nrmse.cal.nnorm
    
    #code below would make an output containing all the spds, not just the distance evaluation metrics
    #calBP = seq(5000,0,-1),
    #ucs.norm.calBP = spd.uncalsample.norm$grid$calBP,
    #ucs.norm.prdens = spd.uncalsample.norm$grid$PrDens,
    #ucs.nnorm.calBP = spd.uncalsample.nnorm$grid$calBP,
    #ucs.nnorm.prdens = spd.uncalsample.nnorm$grid$PrDens,
    #cs.norm.calBP = spd.calsample.norm$grid$calBP,
    #cs.norm.prdens = spd.calsample.norm$grid$PrDens,
    #cs.nnorm.calBP = spd.calsample.nnorm$grid$calBP,
    #cs.nnorm.prdens = spd.calsample.nnorm$grid$PrDens,
  ))
  
}
