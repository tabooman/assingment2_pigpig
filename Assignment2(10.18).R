## Group members : 1.Yitian Zhu 2.Duoer Chen 3.Yaqi Liu
## Contributions : 1.Yitian Zhu 34% 2.Duoer Chen 34% 3.Yaqi Liu 32%
## Explanations : We roughly divided the work as follows:
## Yitian Zhu was mainly responsible for writing the code before function encapsulation.
## Duoer Chen was responsible for bootscrap part of the code writing and modification.
## Yaqi Liu was responsible for writing deconv code and function wrapping.
## But throughout the work we worked together on the project.
## After extensive discussion, debugging and polishing,
## we finally arrived at the current outcome.


# setwd("C:\\Users\\huawei\\Documents\\uoe ex stat programming")

data <- read.table("engcov.txt")
subset_data <- data[1:150,]
death_realnum <- subset_data$nhs

# Define parameters
meanlog <- 3.152
sdlog <- 0.451

# Generate days from 1 to 80
period <- 1:80

# Calculate the probability density for each day
probabilities <- dlnorm(period, meanlog, sdlog)

# Normalised probability
normalized_probabilities <- probabilities / sum(probabilities)

n=29442
duration0 <- sample(period, size = n, replace = TRUE, prob = normalized_probabilities)
death_date0 <- c(rep(62:211, times = death_realnum))
t0 <- death_date0 - duration0

# Make sure infection dates are not before day 1
t0[which(t0 < 0)] <- 1


# Function encapsulation
# deconv: Estimate infection dates based on death data
# Inputs:
# - t: Days with deaths
# - deaths: Number of deaths on each day
# - n.rep: Number of iterations
# - bs: Controls bootstrap
# - t0: Initial infection dates (optional)
# Outputs:
# - P: Best p-values history
# - inft: Infection frequency matrix
# - t0: Updated infection dates
deconv <- function(t,deaths,n.rep,bs,t0){
  
  n=sum(deaths) # Total number of deaths
  
  # Generate original t0 if t0 is NULL
  if (is.null(t0)) {
    duration0 <- sample(period, size = n, replace = TRUE, prob = normalized_probabilities)
    death_date0 <- rep(t, times = deaths)
    t0 <- death_date0 - duration0
    t0[which(t0 < 0)] <- 1
  }
  
  freq_s <- rep(0, 310)
  freq_r <- rep(0, 310)
  freq_r[t] <- deaths
  
  inft <- matrix(0, nrow=310, ncol=n.rep)
  Pearson_history <- numeric(n.rep)
  Pearson_mhistory <- numeric(n.rep)
  t0_move <- t0
  t0_update <- t0
  
  moves_50 <- c(-8, -4, -2,  -1, 1, 2, 4, 8)
  moves_75 <- c(-4, -2, -1, 1, 2, 4)
  moves_100 <- c(-2, -1, 1, 2)
  
  
  for (i in 1:n.rep) {
    freq_r <- rep(0, 310)
    freq_r[t] <- deaths
    if(bs==TRUE){
      for(m in 1:length(freq_r)){
        lamada <- freq_r[m]
        freq_r[m] <- rpois(1,lamada)
      }
      new_n <- sum(freq_r)
      if (new_n < length(t0_move)) {
        t0_move <- sample(t0_move, size = new_n, replace = F)
      } else if (new_n > length(t0_move)) {
        t0_move <- sample(t0_move, size = new_n, replace = T) 
      }
    }
    
    # Calculate new p0 at the start of each iteration
    duration_move <- sample(period, size = sum(freq_r), replace = TRUE, prob = normalized_probabilities)
    death_date <- t0_move + duration_move
    freq_s <- tabulate(death_date,nbins=310)
    Pearson0 <- sum((freq_r - freq_s)^2 / pmax(1, freq_s))
    Pearson_best <- Pearson0 
    
    # Determine move set
    if (i <= 50) {
      moves <- moves_50
    } else if (i <= 75) {
      moves <- moves_75
    } else {
      moves <- moves_100
    }
    
    t0_update <- t0_move
    
    # Update t0_move
    for (j in sample(1:sum(freq_r))) {
      move_j <- sample(moves, 1)
      t0_move[j] <- t0_move[j] + move_j
      t0_move[j] <- max(t0_move[j], 1)
      t0_move[j] <- min(t0_move[j], 310)
      
      freq_s[t0_move[j]+duration_move[j]] <- freq_s[t0_move[j]+duration_move[j]]+1
      freq_s[t0_move[j]-move_j+duration_move[j]] <- freq_s[t0_move[j]-move_j+duration_move[j]]-1
      freq_deathdate_move <- freq_s
      
      Pearson_move <- sum((freq_r - freq_deathdate_move)^2 / pmax(1, freq_deathdate_move))
      
      if (Pearson_move < Pearson_best) {
        Pearson_best <- Pearson_move
        t0_update[j] <- t0_move[j]
      } 
      
      else {
        freq_s[t0_move[j]+duration_move[j]] <- freq_s[t0_move[j]+duration_move[j]]-1
        freq_s[t0_move[j]-move_j+duration_move[j]] <- freq_s[t0_move[j]-move_j+duration_move[j]]+1
        t0_move[j] <- t0_update[j] 
        
      }
      
    }
    
    Pearson_history[i] <- Pearson_best
    Pearson_mhistory[i] <- Pearson_move
    inft[, i] <- freq_s
    print(paste("Iteration:", i, "Updated Pearson_history[i]:", Pearson_best))
    
    
    # Plotting true and predicted values for each iteration
    days <- 1:310
    true_values <- freq_r
    predicted_values <- inft[,i]
    infection_day <- tabulate(t0_update,nbins = 310)
    
    # Draw the lines of TrueValues
    par(mar=c(4, 4, 2, 2))  # Adjust margins
    plot(days, true_values, type = "l", col = "blue", lwd = 2, 
         xlab = "Days", ylab = "Values", main = "Comparison of True and Predicted Values")
    
    # Add dotted lines of PredictedValues and t0_update
    lines(days, predicted_values, col = "red", lwd = 2, lty = 2)
    lines(days, infection_day, col = "green", lwd = 2, lty = 2)
    
    # Add legend
    legend("topright", legend = c("True Values", "Predicted Values", "Infection Day"), 
           col = c("blue", "red","green"), lty = c(1, 2, 2), lwd = 2,cex = 0.5)
    
    
    abline(v=84,lwd = 2, lty = 2,col="black")
  }
  
  return(list(P = Pearson_history, inft = inft, t0 = t0_update,freq_r))
}


# Create the output datasets
dev.new()
result1 <- deconv(t=subset_data$julian, deaths=subset_data$nhs, n.rep=100, bs=FALSE, t0=NULL)
result <- deconv(t=subset_data$julian, deaths=subset_data$nhs, n.rep=100, bs=TRUE, t0=result1$t0)
#deconv(t=subset_data$julian,deaths = subset_data$nhs,n.rep = 100,bs=F, t0=NULL) 
#deconv(t=subset_data$julian,deaths = subset_data$nhs,n.rep = 100,bs=TRUE, t0=t0_update)




  
