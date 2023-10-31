require(ggplot2)
require(reshape2)
require(latex2exp)

# a)

get_corr <- function(rows, cols) {
  # Input:
    # rows: Vector with values of theta_1
    # cols: Vector with values of theta_2
  # Returns: The correlation matrix, corr(theta_1, theta_2)
  m <- length(rows)
  n <- length(cols)
  
  row_matrix <- matrix(rep(rows, n), m, n)
  col_matrix <- t(matrix(rep(cols, m), n, m))
  
  abs_diff <- abs(row_matrix - col_matrix)
  
  corr <- (1 + 15*abs_diff)*exp(-15*abs_diff)
  
  return(corr)
}

get_conditionals <- function(theta_range, obs) {
  # Input:
    # theta_range: a regular grid of theta parameters
    # obs: data frame with pairs of theta and y values, (theta, y(theta))
  # Returns:
    # mean_C: The conditional mean vector
    # Sig_C: The conditional variance matrix
  n <- length(theta_range)
  
  mean_A <-rep(0.5, n)
  
  Sig_AA <- get_corr(theta_range, theta_range)* (0.5^2)

  
  Sig_BB <- get_corr(obs$theta, obs$theta) * (0.5^2)
  
  Sig_AB <- get_corr(theta_range, obs$theta) * (0.5^2)
  
  mean_C <- 0.5 + Sig_AB %*% solve(Sig_BB, obs$y - 0.5)
  
  Sig_C <- Sig_AA - Sig_AB %*% solve(Sig_BB, t(Sig_AB))
  
  return(list(mean_C, Sig_C))
}

run_sim_conditional <- function(theta_range, mean_C, Sig_C, nR=100) {
  # Input:
    # theta_range: a regular grid of theta parameters
    # mean_C: The conditional mean vector
    # Sig_C: The conditional variance matrix
    # nR: Number of realizations for each mean_C
  # Returns:
    # result: Dataframe with the nR simulations
  n = length(mean_C)
  
  # Simulate several realizations
  yMat = matrix(0, nrow = n, ncol = nR)
  Lt = t(chol(Sig_C + diag(n)*1e-8))
  for(i in 1:nR){
    yMat[,i] = Lt%*%rnorm(n) + mean_C
  }
  
  result_sim <- data.frame(theta = theta_range, y=yMat)
  
  return(result_sim)
}

get_prediction_interval <- function(theta_range, mean_C, Sig_C, p=0.90) {
  # Input:
    # theta_range: a regular grid of theta parameters
    # mean_C: The conditional mean vector
    # Sig_C: The conditional variance matrix
    # p: Probability for which to calculate a prediction interval
  # Returns:
    # prediction_interval: The prediction interval for each mean_C,
    # with the given variance from Sig_C,
    # and the probability range given by p
  q <- qnorm((1-p)/2, lower.tail = FALSE)
  prediction_interval <- data.frame(
    theta = theta_range, 
    y = mean_C,
    upper = mean_C + q*sqrt(diag(Sig_C)),
    lower = mean_C - q*sqrt(diag(Sig_C))
  )
  
  return(prediction_interval)
}


obs1 <- data.frame(
  theta = c(0.30, 0.35, 0.39, 0.41, 0.45),
  y = c(0.50, 0.32, 0.40, 0.35, 0.60)
)

theta_range <- seq(0.25, 0.50, 0.005)

conds <- get_conditionals(theta_range, obs1)

mean_C1 <- conds[[1]]
Sig_C1 <- conds[[2]]

result_sim1 <- run_sim_conditional(theta_range, mean_C1, Sig_C1)

prediction_interval1 <- get_prediction_interval(theta_range, mean_C1, Sig_C1)


# Plot
p1 <-   
  ggplot(data = melt(result_sim1, id="theta", value.name = "y")) + 
  geom_point(aes(x = theta, y = y, col="Realizations"), alpha=0.05) +
  geom_point(data = obs1, aes(x=theta, y=y, col="Observations"), alpha=1) +
  ggtitle(TeX("100 realizations of $y(\\theta)$, for a range of $\\theta$, given the 5 evaluation points")) + 
  xlab(TeX("$\\theta$")) +
  ylab(TeX("$y$")) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_color_manual(name = "", values = c("red", "black"))

p1

p2 <-   
  ggplot(data = prediction_interval1) + 
  geom_ribbon(aes(x = theta, ymin = lower, ymax=upper), fill="grey70") +
  geom_line(aes(x = theta, y = y, col="Prediction"), size=0.5) +
  geom_point(data = obs1, aes(x=theta, y=y, col="Observations")) +
  scale_color_manual(name = "Plot", values = c("red", "black")) +
  ggtitle(TeX("Prediction interval for $y(\\theta)$, p=0.90, given the 5 evaluation points")) + 
  xlab(TeX("$\\theta$")) +
  ylab(TeX("$y$")) +
  scale_x_continuous(expand = c(0, 0))


p2

# b)

get_y_lt_prob <- function(theta_range, mean_C, Sig_C, crit=0.30) {
  # Input:
    # theta_range: a regular grid of theta parameters
    # mean_C: The conditional mean vector
    # Sig_C: The conditional variance matrix
    # crit: Upperbound for y, i.e. want y(theta) < crit
  # Returns:
    # goal_result: Dataframe with pairs of theta and Pr[y(theta) < crit]
    
  n <- length(theta_range)
  goal_prob <- pnorm(rep(crit, n), mean = mean_C, sd = sqrt(diag(Sig_C)))
  
  goal_result <- data.frame(theta = theta_range, P = goal_prob)
  
  return(goal_result)
}

goal_result1 <- get_y_lt_prob(theta_range, mean_C1, Sig_C1)

i_max1 <- which.max(goal_result1$P)

# Plot
p3 <-   
  ggplot(data = goal_result1) + 
  geom_step(aes(x = theta, y = P, col="P")) +
  geom_point(data = goal_result1[i_max1,], aes(x=theta, y=P, col="Max")) +
  ggtitle(TeX("Probability that $y(\\theta) < 0.3$, given the 5 evaluation points")) +
  xlab(TeX("$\\theta$")) +
  ylab(TeX("Probability")) +
  scale_color_discrete(name="Plot", labels = c(paste("Max: (", toString(c(goal_result1$theta[i_max1], round(goal_result1$P[i_max1], 2))), ")"), TeX("$Pr[y(\\theta)<0.3]$"))) +
  scale_x_continuous(expand = c(0, 0))

p3

# c)

obs2 <- rbind(obs1, c(0.33, 0.40))

conds <- get_conditionals(theta_range, obs2)

mean_C2 <- conds[[1]]
Sig_C2 <- conds[[2]]

result2 <- run_sim_conditional(theta_range, mean_C2, Sig_C2)

prediction_interval2 <- get_prediction_interval(theta_range, mean_C2, Sig_C2)

goal_result2 <- get_y_lt_prob(theta_range, mean_C2, Sig_C2)

i_max2 <- which.max(goal_result2$P)

# Plot
p4 <-   
  ggplot(data = melt(result2, id="theta", value.name = "y")) + 
  geom_point(aes(x = theta, y = y, col="Realizations"), alpha=0.05) +
  geom_point(data = obs2, aes(x=theta, y=y, col="Observations"), alpha=1) +
  ggtitle(TeX("100 realizations of $y(\\theta)$, for a range of $\\theta$, given the 6 evaluation points")) + 
  xlab(TeX("$\\theta$")) +
  ylab(TeX("$y$")) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_color_manual(name = "", values = c("red", "black"))

p4

p5 <-   
  ggplot(data = prediction_interval2) + 
  geom_ribbon(aes(x = theta, ymin = lower, ymax=upper), fill="grey70") +
  geom_line(aes(x = theta, y = y, col="Prediction"), size=0.5) +
  geom_point(data = obs2, aes(x=theta, y=y, col="Observations")) +
  scale_color_manual(name = "Plot", values = c("red", "black")) +
  ggtitle(TeX("Prediction interval for $y(\\theta)$, p=0.90, given the 6 evaluation points")) + 
  xlab(TeX("$\\theta$")) +
  ylab(TeX("$y$")) +
  scale_x_continuous(expand = c(0, 0))

p5

p6 <-   
  ggplot(data = goal_result2) + 
  geom_step(aes(x = theta, y = P, col="P")) +
  geom_point(data = goal_result2[i_max2,], aes(x=theta, y=P, col="Max")) +
  ggtitle(TeX("Probability that $y(\\theta) < 0.3$, given the 6 evaluation points")) +
  xlab(TeX("$\\theta$")) +
  ylab(TeX("Probability")) +
  scale_color_discrete(name="Plot", labels = c(paste("Max: (", toString(c(goal_result2$theta[i_max2], round(goal_result2$P[i_max2], 2))), ")"), TeX("$Pr[y(\\theta)<0.3]$"))) +
  scale_x_continuous(expand = c(0, 0))

p6