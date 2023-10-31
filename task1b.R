
lam = 5 # birth rate / hour
mu = 6  # death rate / hour

sim = 30               # number of simulations
averageW = rep(0,sim)  # average time in system

for (j in 0:sim) {
  s = c(0)
  x = c(0)                                              # number of patients
  t = c(0)                                              # time
  while(tail(t,1) <= 50*24){                            # while time less than 50 years
    state = tail(x,1)                                   # last state
    if(state == 0){                                     # if last state 0, no death
      s = c(s, rexp(n = 1, rate = lam))                 # sojourn time in state
      t = c(t, tail(t,1) + tail(s,1))                   # add sojourn time
      x = c(x, 1)                                       # increment number of patients
    } else{                                             # last state not 0, death included
      s = c(s, rexp(n = 1, rate = lam+mu))              # sojourn time in state
      t = c(t, tail(t,1) + tail(s,1))                   # add sojourn time
      if(runif(1) < lam/(lam+mu)){                      # probability for 'birth'
        x = c(x, state+1)                               # incrementing number of patients
      } else{
        x = c(x, state-1)                               # finished treating patient
      }
    }
  }
  L = mean(x%*%s/sum(s))                                # long run mean number of patients in UCC
  averageW[j] = L/lam                                   # little's law, find W
}

# Plot for 12 hours
plot(NULL, NULL, main = 'Number of patients in the UCC as a function of time(h)', xlim = c(0, 12), ylim = c(0, max(x[which(t <= 12)])), xlab = "Time (h)", ylab = "Number of patients in UCC", cex.axis = 1.5, cex.lab = 1.5)
for(i in 1:length(which(t <= 12))){
  lines(t[i:(i+1)], rep(x[i],2), lwd = 2)
}

mean(averageW)                                   # estimated expected time in UCC
1/(mu-lam)                                       # analytic expected time in UCC

error = qt(.975, sim-1) * sd(averageW)/sqrt(sim) # confidence interval
mean(averageW) - error                           # min interval
mean(averageW) + error                           # max interval

# Histogram
timeInState = rep(0, max(x))
for(i in 1:(length(x)-1)){
  timeInState[x[i]+1] = timeInState[x[i]+1]+t[i+1]-t[i]
}
sum(timeInState/t[length(x)])                 # normalized histogram
barplot(timeInState/t[length(x)], names.arg = 0:max(x), xlab = "Number of patients in UCC", ylab = "Mean fraction of time", cex.axis = 1.5, cex.lab = 1.5, cex = 1.5)
