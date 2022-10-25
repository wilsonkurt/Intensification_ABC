##This rscript creates generates the 



## New function with time dependent K
## this provides log growth model that will output based upon the chosen parameter values and reading in k at each time step
log_growth_t <- function(N0, r, K, time) {
  if (length(K) != length(time)) {stop("K and time have different legnths")}
  
  N_t <- rep(NA, length(time))
  N_t[1] <- N0
  for (j in 2:length(time)) {
    # N_t[i] = K[i] / (1 + ((K[i] - N0) / N0) * exp(-r * time[i]))
    N_t[j] = N_t[j-1] + r * N_t[j-1] * (1 - (N_t[j-1] / K[j]))
    
  }
  
  return(list(CalBP = rev(time),
              K = K,
              N = N_t,
              PrDens = N_t / sum(N_t)))
}


#Code below would directly emulate Dinapoli et al if inserted in place of the current return() code
#d <- data.frame(CalBP = rev(time), PrDens = (N_t / sum(N_t))) #, K = K, N = N_t
#class(d) <- 'CalGrid'

#return(d)