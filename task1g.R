
lam = 5   # birth rate / hour
mu = 6    # death rate / hour
p = 0.80  # probability that arrived patient is urgent

sim = 30                  # number of simulations
averageWU = rep(0,sim)    # average time for urgent patient in system
averageWN = rep(0,sim)    # average time for normal patient in system

for (j in 0:sim) {
  u = c(0)                                              # number of urgent patients
  n = c(0)                                              # number of normal patients
  t = c(0)                                              # time
  s = c(0)                                              # sojourn time
  
  while(tail(t,1) <= 50*24){                            # while time less than 50 days
    stateU = tail(u,1)                                  # last state urgent patients
    stateN = tail(n,1)                                  # last state normal patients
    
    if(stateU == 0 && stateN == 0){                     # if last state 0, no death
      s = c(s, rexp(n = 1, rate = lam))                 # sojourn time in state
      t = c(t, tail(t,1) + tail(s,1))                   # add sojourn time
      if (runif(1) < p){                                
        u = c(u, 1)                                     # arrived urgent patient
        n = c(n, stateN)
      } else{
        u = c(u, stateU)
        n = c(n, 1)                                     # arrived normal patient
      }
    } else{                                             # last state not 0, 'death' included
      s = c(s, rexp(n = 1, rate = lam+mu))              # sojourn time in state
      t = c(t, tail(t,1) + tail(s,1))                   # add sojourn time
      
      if(runif(1) < lam/(lam+mu)){                      # probability for 'birth'
        if(runif(1) < p){                               # the arrived patient is urgent
          u = c(u, stateU+1)                            # increment number of urgent patients
          n = c(n, stateN)
        } else{                                         # the arrived patient is normal
          u = c(u, stateU)
          n = c(n, stateN+1)                            # increment number of normal patients
        }
      } else{
        if(tail(u,1) > 0){                              # treating urgent patients, if there are any
          u = c(u, stateU-1)                            # finished treating urgent patient
          n = c(n, stateN)
        } else{
          u = c(u, stateU)
          n = c(n, stateN-1)                            # finished treating normal patient
        }
      }
    }
  }
  LU = mean(u%*%s/sum(s))                               # long run mean number of urgent patients
  LN = mean(n%*%s/sum(s))                               # long run mean number of normal patients
  
  averageWU[j] = LU/(p*lam)                             # long run expected time in UCC for urgent patients
  averageWN[j] = LN/((1-p)*lam)                         # long run expected time in UCC for normal patients
}


# Plot for 12 hours
plot(NULL, NULL, main = 'Number of urgent and normal patients in the UCC as a function of time(h)', xlim = c(0, 12), ylim = c(0, max(c(u[which(t <= 12)], n[which(t <= 12)]))), xlab = "Time (h)", ylab = "Number of patients in UCC", cex.axis = 1.5, cex.lab = 1.5)
for(i in 1:length(which(t <= 12))){
  lines(t[i:(i+1)], rep(u[i],2), lwd = 2, col = 'red')
  lines(t[i:(i+1)], rep(n[i],2), lwd = 2, col = 'blue')
}
legend(0.5, 8, legend=c('Urgent patients', 'Normal patients'), col = c('red', 'blue'), lty=1:1, cex=1)


# Calculate expected time for urgent patient in UCC
mean(averageWU)                                   # estimated expected time for urgent patient in UCC
1/(mu-p*lam)                                      # analytic expected time for urgent patient in UCC

error = qt(.975, sim-1) * sd(averageWU)/sqrt(sim) # confidence interval
mean(averageWU) - error                           # min interval
mean(averageWU) + error                           # max interval


# Calculate expected time for normal patient in UCC
mean(averageWN)                                   # estimated expected time for normal patient in UCC
mu/((mu-lam)*(mu-p*lam))                          # analytic expected time for normal patient in UCC

error = qt(.975, sim-1) * sd(averageWN)/sqrt(sim) # confidence interval
mean(averageWN) - error                           # min interval
mean(averageWN) + error                           # max interval


# Histogram for urgent patients in UCC
timeInState = rep(0, max(u))
for(i in 1:(length(u)-1)){
  timeInState[u[i]+1] = timeInState[u[i]+1]+t[i+1]-t[i]
}
timeInState
sum(timeInState/t[length(u)])                 # normalized histogram

barplot(timeInState/t[length(u)], names.arg = 0:max(u), xlab = "Number of urgent patients in UCC", ylab = "Mean fraction of time", cex.axis = 1.5, cex.lab = 1.5, cex = 1.5)
