setwd("/Users/zhuyitian/Desktop/UOE_studymaterial/statistical_programming/assignment2")
data <- read.table("engcov.txt")
subset_data <- data[1:150,]
death_realdate <- subset_data$julian
death_realnum <- subset_data$nhs

# 定义参数
meanlog <- 3.152
sdlog <- 0.451

# 生成从 1 到 80 的天数
period <- 1:80

# 计算每一天的概率密度
probabilities <- dlnorm(period, meanlog, sdlog)

# 归一化概率
normalized_probabilities <- probabilities / sum(probabilities)

n=29442
duration0 <- sample(period, size = n, replace = TRUE, prob = normalized_probabilities)
death_date0 <- c(rep(62:211, times = death_realnum))
t0 <- death_date0 - duration0

# 将负的初始化感染日期变成为第一天
t0[which(t0 < 0)] <- 1


freq_s <- rep(0, 310)
freq_r <- rep(0, 310)

freq_stimulation <- tabulate(death_date)
freq_s[1:length(freq_stimulation)] <- freq_stimulation

freq_r[62:211] <- death_realnum

n.rep <- 100
inft <- matrix(0, nrow=310, ncol=n.rep)
freq_s <- rep(0, 310)
freq_s1 <- rep(0,310)
p_history <- numeric(n.rep)
p_mhistory <- numeric(n.rep)
t0_move <- t0
t0_update <- t0

moves_50 <- c(-16, -8, -4, -2, 2, 4, 8, 16)
moves_75 <- c(-4, -2, -1, 1, 2, 4)
moves_100 <- c(-2, -1, 1, 2)

for (i in 1:n.rep) {
  # 在每次 i 循环开始前计算新的 p0
  duration_move <- sample(period, size = n, replace = TRUE, prob = normalized_probabilities)
  death_date <- t0_move + duration_move
  freq_s <- tabulate(death_date,nbin=310)
  p0 <- sum((freq_r - freq_s)^2 / pmax(1, freq_s))
  p_best <- p0 
  
  if (i <= 50) {
    moves <- moves_50
  } else if (i <= 75) {
    moves <- moves_75
  } else {
    moves <- moves_100
  }
  
  t0_update <- t0_move
  for (j in sample(1:n)) {
    move_j <- sample(moves, 1)
    t0_move[j] <- t0_move[j] + move_j
    t0_move[j] <- max(t0_move[j], 1)
    t0_move[j] <- min(t0_move[j], 310)
    
    freq_s[t0_move[j]+duration_move[j]] <- freq_s[t0_move[j]+duration_move[j]]+1
    freq_s[t0_move[j]-move_j+duration_move[j]] <- freq_s[t0_move[j]-move_j+duration_move[j]]-1
    freq_deathdate_move <- freq_s
    
    p_move <- sum((freq_r - freq_deathdate_move)^2 / pmax(1, freq_deathdate_move))
    
    if (p_move < p_best) {
      p_best <- p_move
      t0_update[j] <- t0_move[j]
    } 
    
    else {
      freq_s[t0_move[j]+duration_move[j]] <- freq_s[t0_move[j]+duration_move[j]]-1
      freq_s[t0_move[j]-move_j+duration_move[j]] <- freq_s[t0_move[j]-move_j+duration_move[j]]+1
      t0_move[j] <- t0_update[j] 
      
    }
    
  }
  
  p_history[i] <- p_best
  p_mhistory[i] <- p_move
  inft[, i] <- freq_s
  print(p_best)
}



#创建图的数据集
days <- 1:310
data <- data.frame(Days = days, 
                   TrueValues = freq_r, 
                   PredictedValues = inft[,100])

# Plot the graph
ggplot(data, aes(x = Days)) +
  geom_line(aes(y = TrueValues, color = "True Values"), size = 1) +
  geom_line(aes(y = PredictedValues, color = "Predicted Values"), size = 1, linetype = "dashed") +
  labs(title = "Comparison of True and Predicted Values",
       x = "Days",
       y = "Values") +
  scale_color_manual(values = c("True Values" = "blue", "Predicted Values" = "red")) +
  theme_minimal() +
  theme(legend.title = element_blank())

