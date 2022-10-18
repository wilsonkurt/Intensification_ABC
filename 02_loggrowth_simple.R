## Single random rin of log growth model
## Sample the following
## tstart (initial settlement)
## npp_lo: minimum NPP for growth
## npp_hi: maximum NPP for growth
## cK0: initial number of individuals
## r: growth rate
## Kmax: max carrying capacity
## NPPlo: NPP value corresponding to K=0
## NPPhi: NPP value corresponding to K=Kmax
set.seed(42)

library(ggplot2)

log_growth <- function(N0, r, K, time) {
  N_t <- rep(NA, length(time))
  N_t <- N0
  for (i in 2:length(time)) {
    # N_t[i] = K / (1 + ((K - N0) / N0) * exp(-r * time[i]))
    N_t[i] = N_t[i-1] + r * N_t[i-1] * (1 - (N_t[i-1] / K))
  }
  return(list(time = time,
              N = N_t))
}

## Simple example
myN0 <- 2
myr <- 0.05
myK <- 500
mytime <- seq(1, 500)
results <- log_growth(N0 = myN0, 
                      r = myr, 
                      K = myK, 
                      time = mytime)

plot(results$time, results$N, 
     type = 'l')

## Try with random parameters
clim.df <- read.csv("./data/clim.df.csv")

tstart <- round(runif(1, min = 3000, max = 5000))
npp_limits <- sort(runif(2, min = 0.1, max = 0.8))
npp_lo <- npp_limits[1]
npp_hi <- npp_limits[2]
cK0 <- round(runif(1, min = 2, max = 10))
r <- exp(runif(1, min = log(0.01), max = log(0.5)))
Kmax <- exp(runif(1, min = log(100), max = log(10000)))
## Kmax: max carrying capacity

## NPP - K relationship
plot.df <- data.frame(npp = c(min(clim.df$npp), npp_lo, npp_hi, max(clim.df$npp)),
                      K = c(0, 0, Kmax, Kmax))
ggplot(plot.df, aes(x = npp, y = K)) +
  geom_line()

## Derive time-dependent K
## First filter climate to tstart
clim.df <- clim.df[clim.df$age < tstart, ]
time <- max(clim.df$age) - clim.df$age## Need to invert time somehow
K_t <- approx(plot.df$npp, plot.df$K, clim.df$npp)$y
plot(clim.df$age, K_t)

## New function with time dependent K
log_growth_t <- function(N0, r, K, time) {
  if (length(K) != length(time)) {stop("K and time have different legnths")}
  
  N_t <- rep(NA, length(time))
  N_t[1] <- N0
  for (i in 2:length(time)) {
    # N_t[i] = K[i] / (1 + ((K[i] - N0) / N0) * exp(-r * time[i]))
    N_t[i] = N_t[i-1] + r * N_t[i-1] * (1 - (N_t[i-1] / K[i]))
    
  }
  return(list(time = time,
              K = K,
              N = N_t))
}

results <- log_growth_t(N0 = cK0, 
                        r = r, 
                        K = K_t, 
                        time = time)
plot(results$time, results$N, 
     type = 'l')
