#####################################
# Ichihara - Iterative SD algorithm #
#####################################

#############
# Functions #
#############

# Resulting correction factor
ichihara.correct <- function(k){
  c <- rnorm(100000)
  c <- c[c > - k & c <= k]
  return(1/sd(c))
}

# Recursively defined function 
ichihara.iterate <- function(data, mean = NA, sd = NA, step = 0, k = k, c.k = c.k){
  
  # Case 1: First step -> find mean and sd!
  if(is.na(sd) | is.na(mean)){
    
    return(ichihara.iterate(data = data, mean = mean(data), sd = sd(data), 
                            step = step + 1, k = k, c.k = c.k))
  }
  
  # Filter
  sd.prev <- sd(data)
  data <- data[data > mean - k * sd & data <= mean + k * sd]
  sd.new <- sd(data)
  
  # if(plot){print(blot(step = step, mean.it = mean(data), 
  #                     sd.it = c.k * sd(data)))}
  
  # Case 2: Last step (equilibrium)
  if( (sd.prev-sd.new)/sd.prev < 1.0E-5){
    #cat("Ichihara algo, final step: -", step ,"-\n")
    return(c(mean(data), c.k*sd(data), step, length(data)))
  }
  
  return(ichihara.iterate(data = data, mean = mean(data), sd = c.k * sd(data),
                          step = step + 1, k = k, c.k = c.k))
  
}

#############
# Algorithm #
#############

ichihara <- function(data, refQt = 0.95, k = 1.4, CI = T, R = 500){
  
  Result = list()
  x <- NULL
  
  Result$call <- match.call()
  Result$refQt <- refQt
  Result$k <- k
  Result$CI <- CI
  Result$R <- R
  
  # Truncation factor (mean +/- k * sd)
  Result$k <- k
  Result$c.k <- ichihara.correct(k)
  Result$prop <- pnorm(k)-pnorm(-k)
  
  # Correction factor
  c.k <- ichihara.correct(k)
  
  # Run iterative function from above
  x <- data[!is.na(data) & data != 0]
  run <- ichihara.iterate(data = x, k = k, c.k = c.k)
  
  Result$mean <- run[1]
  Result$sd <- run[2]
  Result$final.iter <- run[3]
  Result$final.length <- run[4]
  Result$final.excluded <- length(x)-run[4]
  
  # Reference limits
  Result$lower.limit <- qnorm(p = (1 - refQt) / 2, mean = run[1], sd = run[2])
  Result$upper.limit <- qnorm(p = 1 - ((1 - refQt) / 2), mean = run[1], sd = run[2])
  
  if(CI){
    # https://en.wikipedia.org/wiki/Reference_range#Confidence_interval_of_limit
    # se = sqrt(((sd^2)/n) + (((z.score.RI^2)*(sd^2))/(2*n)));
    z.score.RI <- qnorm(1 - ((1 - refQt) / 2))
    z.score.CI <- qnorm(1 - ((1 - 0.9) / 2))
    st.error <- sqrt((run[2]^2/run[4]) + (z.score.RI^2*(run[2]^2))/(2*run[4]))
    
    Result$lower.limit.low <- Result$lower.limit - z.score.CI * st.error
    Result$lower.limit.high <- Result$lower.limit + z.score.CI * st.error
    Result$upper.limit.low <- Result$upper.limit - z.score.CI * st.error
    Result$upper.limit.high <- Result$upper.limit + z.score.CI * st.error
    
  }else{
    Result$lower.limit.low <- Result$lower.limit.high <- NA
    Result$upper.limit.low <- Result$upper.limit.high <- NA
  }
  
  Result$RLs <- c(Result$lower.limit.low, Result$lower.limit, Result$lower.limit.high,
                  Result$upper.limit.low, Result$upper.limit, Result$upper.limit.high)
  
  Result$RLs[Result$RLs < 0] <- 0
  
  if(is.na(run[1]) | is.na(run[1])){ Result$fail = TRUE }else{ Result$fail = F }
  
  return(Result)
}

# if (!interactive()) {
#   
#   # Testing
#   data <- sample(refineR::testcase1, size = 2000, replace = F)
#   S <- ichihara(data = data, CI = T)
#   
#   # Plotting
#   hist(data, breaks = 50, freq = F)
#   abline(v = S$RLs, col = 2)
#   
#   xx <- seq(0,20,by=0.01)
#   lines(x=xx, dnorm(xx, mean = S$mean, sd = S$sd))
#   
# }

