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

n = 29442
duration0 <- sample(period, size = n, replace = TRUE, prob = normalized_probabilities)
death_date0 <- c(rep(62:211, times = death_realnum))
t0 <- death_date0 - duration0

# 将负的初始化感染日期变成为第一天
t0[which(t0 < 0)] <- 1

# 初始化真实死亡人数和模拟死亡人数
freq_s <- rep(0, 310)
freq_r <- rep(0, 310)
freq_r[62:211] <- death_realnum

# 初始化模拟死亡日期的频率
duration_new <- sample(period, size = n, replace = TRUE, prob = normalized_probabilities)
death_date <- t0 + duration_new
freq_s <- tabulate(death_date, nbins = 310)

# 参数设置
n.rep <- 100
inft <- matrix(0, nrow = 310, ncol = n.rep)
p_history <- numeric(n.rep)
t0_move <- t0
t0_update <- t0

# 定义不同阶段的移动范围
moves_50 <- c(-16, -8, -4, -2, 2, 4, 8, 16)
moves_75 <- c(-4, -2, -1, 1, 2, 4)
moves_100 <- c(-2, -1, 1, 2)

for (i in 1:n.rep) {
  # 在每次 i 循环开始前生成一次 duration_new
  duration_new <- sample(period, size = n, replace = TRUE, prob = normalized_probabilities)
  death_date <- t0_move + duration_new
  freq_s <- tabulate(death_date, nbins = 310)
  
  # 计算当前 p0
  p0 <- sum((freq_r - freq_s)^2 / pmax(1, freq_s))
  p_best <- p0 
  
  # 根据迭代次数调整移动范围
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
    old_t0 <- t0_move[j]
    new_t0 <- old_t0 + move_j
    new_t0 <- max(min(new_t0, 310), 1)  # 确保新值在1到310之间
    
    # 更新频率数组 freq_s
    freq_s[old_t0 + duration_new[j]] <- freq_s[old_t0 + duration_new[j]] - 1
    freq_s[new_t0 + duration_new[j]] <- freq_s[new_t0 + duration_new[j]] + 1
    
    # 计算新的 p 值
    p_move <- sum((freq_r - freq_s)^2 / pmax(1, freq_s))
    
    if (p_move < p_best) {
      p_best <- p_move
      t0_update[j] <- new_t0
    } else {
      # 如果 p 值不更好，回滚 freq_s
      freq_s[new_t0 + duration_new[j]] <- freq_s[new_t0 + duration_new[j]] - 1
      freq_s[old_t0 + duration_new[j]] <- freq_s[old_t0 + duration_new[j]] + 1
    }
  }
  
  # 更新 t0_move 为 t0_update
  t0_move <- t0_update
  p_history[i] <- p_best
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

