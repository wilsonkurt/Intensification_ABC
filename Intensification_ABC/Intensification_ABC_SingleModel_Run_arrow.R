## ----setup, include=FALSE-----------------------------------
#knitr::opts_chunk$set(echo = TRUE)

start_time <- Sys.time()
## -----------------------------------------------------------
#C:/Users/u1132502/Box Sync/ABC_Work/CHPC/Intensification_ABC
library(spatstat.utils)
#library(doParallel)#parallel processing
#library(doSNOW)#parallel processing
library(rcarbon)#14c date calibration and spd recreation
#library(foreach)
library(truncnorm)#used in making prior probability distribution of r
library(arrow)#used to speed up spd process through use of an arrow table
library(duckdb)#used to work with arrow data table format in windows
library(dplyr)#used in data table processing to speed up computation
library(tidyr)#used in data table processing to speed up computation

#source("/uufs/chpc.utah.edu/common/home/u1132502/Desktop/Intensification_ABC/03_loggrowthmodel.R")#creates and enables the log growth function for simulation
#source("/uufs/chpc.utah.edu/common/home/u1132502/Desktop/Intensification_ABC/04_simulatedmodelspds.R")#generates SPDs from the simulated population dynamic

source("C:/Users/u1132502/Box Sync/ABC_Work/CHPC/Intensification_ABC/03_loggrowthmodel.R")#creates and enables the log growth function for simulation
source("C:/Users/u1132502/Documents/GitHub/Intensification_ABC/Intensification_ABC/04_simulatedmodelspds_singlemodel_run_arrow.R")#generates SPDs from the simulated population dynamic


## -----------------------------------------------------------
#read in the climate data for the Colorado Plateau
load("C:/Users/u1132502/Box Sync/ABC_Work/CHPC/Intensification_ABC/data/pla.clim.df.Rdata")#read the 4ka climate data with npp calculated using Miami model

#pla.clim.df <- read.csv("/uufs/chpc.utah.edu/common/home/u1132502/Desktop/Intensification_ABC/data/pla.clim.df.csv")#read the 4ka climate data with npp calculated using Miami model
pla.clim.df <- approx(pla.clim.df$year, pla.clim.df$npp, xout = seq(4000,0,-1))#interpolate npp between each of the decadal reconstructions for annual resolution
pla.clim.df <- data.frame(age = pla.clim.df$x, npp = pla.clim.df$y)

#repeat for Great Basin
load("C:/Users/u1132502/Box Sync/ABC_Work/CHPC/Intensification_ABC/data/bas.clim.df.Rdata")#read the 4ka climate data with npp calculated using Miami model

#bas.clim.df <- read.csv("/uufs/chpc.utah.edu/common/home/u1132502/Desktop/Intensification_ABC/data/bas.clim.df.csv")#read the 4ka climate data with npp calculated using Miami model
bas.clim.df <- approx(bas.clim.df$year, bas.clim.df$npp, xout = seq(4000,0,-1))#interpolate npp between each of the decadal reconstructions for annual resolution
bas.clim.df <- data.frame(age = bas.clim.df$x, npp = bas.clim.df$y)

#read in the saved spds
load("C:/Users/u1132502/Box Sync/ABC_Work/CHPC/Intensification_ABC/data/Observed_spds_Plateau.Rdata")

#observed.spds.plat <- read.csv("/uufs/chpc.utah.edu/common/home/u1132502/Desktop/Intensification_ABC/data/Observed_spds_Plateau.csv")
observed.spds.plat.norm <- Observed_spds_Plateau[,c(1,seq(2,2000,2))]#grab the normalized spds and put in own object
observed.spds.plat.nnorm <- Observed_spds_Plateau[,c(1,seq(3,2001,2))]#make a df for the not normalized spds

load("C:/Users/u1132502/Box Sync/ABC_Work/CHPC/Intensification_ABC/data/Observed_spds_Basin.Rdata")

#observed.spds.basin <- read.csv("/uufs/chpc.utah.edu/common/home/u1132502/Desktop/Intensification_ABC/data/Observed_spds_Basin.csv")
observed.spds.basin.norm <- Observed_spds_Basin[,c(1,seq(2,2000,2))]#grab the normalized spds and put in own object
observed.spds.basin.nnorm <- Observed_spds_Basin[,c(1,seq(3,2001,2))]#make a df for the not normalized spds

#read in the errors for use in simulated spds
#Plateau.errors <- read.csv("/uufs/chpc.utah.edu/common/home/u1132502/Desktop/Intensification_ABC/data/Plateau_errors.csv")
#Basin.errors <- read.csv("/uufs/chpc.utah.edu/common/home/u1132502/Desktop/Intensification_ABC/data/Basin_errors.csv")

load("C:/Users/u1132502/Box Sync/ABC_Work/CHPC/Intensification_ABC/data/Plateau_errors.Rdata")
load("C:/Users/u1132502/Box Sync/ABC_Work/CHPC/Intensification_ABC/data/Basin_errors.Rdata")

#set the number of bins per region
allbins <- 3192 #this is the # unique sites/bins for Plateau dates 4000 ybp and younger
allbins.basin <- 1364 #this is the # unique sites/bins for Basin dates 4000 ybp and younger

#set the number of years (equating to time period) for model matching (here using 4000 to 1000 ybp - see below)
yrs.match.range <- 1:3001



## -----------------------------------------------------------
plat.spd.avgs <- data.frame(year = observed.spds.plat.norm$calBP, spd.avg = rowMeans(observed.spds.plat.norm[,2:1001]))
cp.delta <- mean(plat.spd.avgs[,2]) * 0.05

bas.spd.avgs <- data.frame(year = observed.spds.basin.norm$calBP, spd.avg = rowMeans(observed.spds.basin.norm[,2:1001]))
gb.delta <- mean(bas.spd.avgs[,2]) * 0.05


## -----------------------------------------------------------
weights.vec <- c(rep(1, times = 1500), rep(10, times = 1500))


## -----------------------------------------------------------

#Globals
nsim <- 1 #1000000 models total, but 1 param set each model #set the number of parameter combinations to use
tstart <- 4000 #set the beginning of the simulation to be 4000 yBP

#Landscape variables
npp_lo <- vector()#a vector to hold the npp min threshold for any population to be present
npp_hi <- vector()#a vector to hold the npp max threshold over which no increase in k occurs (unless by intensification)
for (i in 1:nsim) {
  npp_limits <- sort(runif(2, min = 0.17, max = 1.02))#grab 2 npp values to be low and high and sort from low to high
  #the min npp during our study period is 0.4763764 while the max is 0.7187248 We therefore randomly sample from values below the min to above the max. Min npp represents the min npp needed to support at least some population. Max represents the npp value at which the maximum k (without intensification) is reached. Any npp above that is essentially unusable by the population and doesn't result in larger k.
  npp_lo <- append(npp_lo, npp_limits[1])  
  npp_hi <- append(npp_hi, npp_limits[2])
}

#Intensification variables
tseistart <- round(runif(nsim, min = 1800, max = 3500))#select the year in which intensification will begin
SEI_max <- (rexp(nsim, 0.10)) #distribution produces a right skewed distribution with a mean ~10 and a tail up to ~60 (comparable to the ethnographic cross-cultural estimate of magnitude of population density difference from WNAI for agr vs non-agr groups)
K_m <- 1 #this is the intensification constant applied to K_t on each time step. It is 1 b/c prior to SEI beginning, there is no SEI effect on K_t. But once SEI begins, it is a multiplier, so 1 + SEI at that moment will add to k.

#Population variables
cK0 <- rexp(nsim, rate = 20)#set a starting population as a proportion of Kmax
r <- rtruncnorm(n=nsim,a=0.0001,b=0.5,mean=0.04,sd=0.1)#set a population growth rate. Mean 0.04 from Zahid et al 2016 and Page et al. 2016.
#Range enables sampling across other estimates in published lit and from DiNappoli and Porcic ABC papers.
Kmax <- 1 #hold Kmax constant at a value of 1, with only SEI allowed to push population beyond that amount
#Kmax <- exp(runif(nsim, min = log(100), max = log(10000)))#set a maximum carrying capacity

#make a dataframe to hold nsim parameter combinations
params <- data.frame(tstart = tstart, r = r, cK0 = cK0, Kmax = Kmax, npp_lo = npp_lo, npp_hi = npp_hi, tseistart = tseistart, SEI_max = SEI_max)

params[,9:4009] <- NA
colnames(params) <- c(colnames(params[,1:8]), seq(4000,0,-1))

idCounter <- 0



#set objects to be used in pattern matching to plateau data
region.errors <- Plateau_errors$Errors
bins.use <- allbins
observed.spds.norm <- observed.spds.plat.norm
observed.spds.nnorm <- observed.spds.plat.nnorm
delta.val <- cp.delta


    #establish npp boundary conditions
    npp.conds <- data.frame(npp = c(params$npp_lo, params$npp_hi, max(pla.clim.df$npp)),
                      K = c(0, params$Kmax, params$Kmax)) #dataframe where k at npp_lo = 0 (min threshold for population survival - anything at or below this will have a k of 0), at npp_hi = max k, and at highest observed npp = maxk - this is all used below for the npp to k relationship for the model run
    
    ## Derive time-dependent K
    ## this will give us k at each time step (or at least k without any impact of SEI - so k based on npp)
    ## First filter climate to tstart
    clim.df.run <- pla.clim.df[pla.clim.df$age <= tstart & pla.clim.df$age >= 0,] #get climate data from 4000 to 0 yBP
    time <- clim.df.run$age
    #time <- max(clim.df$age) - clim.df$age## Need to invert time
    K_t <- approx(npp.conds$npp, npp.conds$K, clim.df.run$npp)$y #get a list of interpolated points connecting npp and k at each observation of npp from our climate reconstruction. This gives us the carrying capacity at each time step w/out any SEI occurring.
    K_t[is.na(K_t)] <- 0 #Make sure any periods where popluation should fall to 0 are 0s and not NAs
    
    #set.seed(i)
    sim.growthModel <- data.frame(log_growth_t(N0 = params$cK0, r = params$r, K = K_t, K_m = K_m, time = time, tseistart = params$tseistart[i], SEI_max = params$SEI_max[i])) #run the log growth model function

    sim.growthModel.spd <- sim.growthModel[c(1,4)] #select year and simulated PrDens for the calgrid object
    class(sim.growthModel.spd) <- 'CalGrid' #change to a calgrid object
    
    results <- createmodelspd(model = sim.growthModel.spd)
    params[,9:4009] <- sim.growthModel$SEI
    
    model.result <- cbind(params, results)
   
    cp.file <- paste("C:/Users/u1132502/Box Sync/ABC_Work/CHPC/Intensification_ABC/Output/CP/CP_Model_results.arrowtest", Sys.getpid(), ".Rdata", sep="")
    save(model.result, file = cp.file)

# Great Basin simulations

#set objects to be used in pattern matching to basin data
region.errors <- Basin_errors$Errors
bins.use <- allbins.basin
observed.spds.norm <- observed.spds.basin.norm
observed.spds.nnorm <- observed.spds.basin.nnorm
delta.val <- gb.delta


    #establish npp boundary conditions
    npp.conds <- data.frame(npp = c(params$npp_lo, params$npp_hi, max(bas.clim.df$npp)),
                      K = c(0, params$Kmax, params$Kmax)) #dataframe where k at npp_lo = 0 (min threshold for population survival - anything at or below this will have a k of 0), at npp_hi = max k, and at highest observed npp = maxk - this is all used below for the npp to k relationship for the model run
    
    ## Derive time-dependent K
    ## this will give us k at each time step (or at least k without any impact of SEI - so k based on npp)
    ## First filter climate to tstart
    clim.df.run <- bas.clim.df[bas.clim.df$age <= tstart & bas.clim.df$age >= 0,] #get climate data from 4000 to 0 yBP
    time <- clim.df.run$age
    #time <- max(clim.df$age) - clim.df$age## Need to invert time
    K_t <- approx(npp.conds$npp, npp.conds$K, clim.df.run$npp)$y #get a list of interpolated points connecting npp and k at each observation of npp from our climate reconstruction. This gives us the carrying capacity at each time step w/out any SEI occurring.
    K_t[is.na(K_t)] <- 0 #Make sure any periods where popluation should fall to 0 are 0s and not NAs
    
    #set.seed(i)
    sim.growthModel <- data.frame(log_growth_t(N0 = params$cK0, r = params$r, K = K_t, K_m = K_m, time = time, tseistart = params$tseistart[i], SEI_max = params$SEI_max[i])) #run the log growth model function

    sim.growthModel.spd <- sim.growthModel[c(1,4)] #select year and simulated PrDens for the calgrid object
    class(sim.growthModel.spd) <- 'CalGrid' #change to a calgrid object
    

    gb.results <- createmodelspd(model = sim.growthModel.spd)
    params[,9:4009] <- sim.growthModel$SEI
    
    gb.model.result <- cbind(params, gb.results)
    
    gb.file <- paste("C:/Users/u1132502/Box Sync/ABC_Work/CHPC/Intensification_ABC/Output/GB/GB_Model_results.arrowtest", Sys.getpid(), ".Rdata", sep="")
    save(gb.model.result, file = gb.file)
    
end_time <- Sys.time()
run_time <- end_time - start_time
run_time

