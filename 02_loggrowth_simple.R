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

library(ggplot2)

log_growth <- function(N0, r, K, time) {
  N_t <- rep(NA, length(time))
  N_t <- N0
  for (i in 2:length(time)) {
    N_t[i] = K / (1 + ((K - N0) / N0) * exp(-r * time[i]))
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
