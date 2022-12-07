##This rscript creates generates the simulated population growth based on sampled parameters.

## New function with time dependent K
## this provides log growth model that will output based upon the chosen parameter values and reading in k at each time step. This first
## set of code does not include intensification.

# log_growth_t <- function(N0, r, K, time) {
#   if (length(K) != length(time)) {stop("K and time have different legnths")}
#   
#   N_t <- rep(NA, length(time))
#   N_t[1] <- N0
#   for (j in 2:length(time)) {
#     # N_t[i] = K[i] / (1 + ((K[i] - N0) / N0) * exp(-r * time[i]))
#     N_t[j] = N_t[j-1] + r * N_t[j-1] * (1 - (N_t[j-1] / K[j]))
#     
#   }
#   
#   return(list(CalBP = rev(time),
#               K = K,
#               N = N_t,
#               PrDens = N_t / sum(N_t)))
# }


#Code below would directly emulate Dinapoli et al if inserted in place of the current return() code
#d <- data.frame(CalBP = rev(time), PrDens = (N_t / sum(N_t))) #, K = K, N = N_t
#class(d) <- 'CalGrid'
#return(d)

log_growth_t <- function(N0, r, K, time, tseistart, K_m, SEI_max) {
  if (length(K) != length(time)) {stop("K and time have different legnths")}
  
  N_t <- rep(NA, length(time)) #make a vector of population size for each year that is NA to begin
  N_t[1] <- N0 #set initial population to starting population
  SEI_now <- 0 #start SEI as no SEI
  SEI_year <- rep(NA, length(time)) #make a vector that records the amount of SEI that happened in each year
  SEI_year[1] <- 0 #start vector at 0
  
  for (j in 2:length(time)) {
    
    if (time[j] <= tseistart & N_t[j-1] >= (K[j-1] * 0.9) & SEI_now < (SEI_max - 1)) { #if the simulation is at or beyond the point where SEI may
      #begin, and there is a population crunch/crisis happening b/c population is near k, and if the amount of SEI that has occurred
      #is not yet the maximum, then set a new SEI amount to occur on this turn
      SEI_now <- runif(1, min = SEI_now, (SEI_max - 1)) #We do SEI_max - 1 b/c the K_m (SEI constant) is 1 and we increase from there.
    }
    
    #simulate population, each time saving the output population size
    N_t[j] <- N_t[j-1] + (r * N_t[j-1] * (1 - (N_t[j-1] / (K[j] * (K_m + SEI_now))))) 
    #within the k code (K[j] * (K_m + SEI_now[j])) this is saying set the carrying capacity for this year/turn to be the 
    #carrying capacity from climatically driven NPP times the bonus from the amount of intensification that has occurred
    
    #If however, population has grown beyond the max carrying capacity even with SEI, drop pop to the max K with SEI factored in
    if (N_t[j] > (K[j] * (K_m + SEI_now))) {N_t[j] <- (K[j] * (K_m + SEI_now))}
    
    SEI_year[j] <- if (SEI_now == 0) {SEI_now} else {(SEI_now + 1)} #record the SEI that happened this year adding + 1 if SEI happened
    #to record if the multiplier is greater than 1.
  }
  
  return(list(CalBP = time,
              K = K,
              N = N_t,
              PrDens = N_t / sum(N_t),
              SEI = SEI_year))
}

## this set of code incorporates SEI into the popluation simulation


