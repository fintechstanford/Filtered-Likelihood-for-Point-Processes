##############################################################################################################
### This code computes an approximation of the log-likelihood in equation (5) of the paper "Filtered Likelihood for Point Processes" by Kay Giesecke and Gustavo Schwenkler. The copyright to these codes remains with the authors.
##############################################################################################################

# Enter the parameters
alpha = 3
beta = 3
gamma = 2
kappa = 2.5
gamma1 = -0.5
gamma2 = 1
gamma3 = 0.2
gamma4 = 5
gamma5 = -1/8
gamma6 = 0.5
par = c(alpha, beta, gamma, kappa, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6)

Y0 = 2.8
Z0 = 1
epsilon = 0.01

# Fix the time horizon
tau = 40

# True likelihood
trueloglik = 8655

# Load all relevenat algorithms
source("algorithms.R")
source("input_probabilities.R")

# These are some auxiliary time parameters 
freq_zeta = 21;
day_freq = 1;
dpy = 250;
dpw = 5;
hpd = 20;
mpd = hpd*60;
tau2 = tau*dpy
precision = day_freq/hpd;

# Load the simulated event and factors data
ep = read.table("simulated_events.txt", header = FALSE, skip = 1)
exp_f = read.table("simulated_factors.txt", header = FALSE, skip = 1)

# Only the fist factor is observed, the second is a frailty. The vector col corresponds to the factor Y in the paper. It is observed one a month
obs = seq(from = 1, to = nrow(exp_f), by = freq_zeta*hpd)
cov = exp_f[obs,1:2]

# Construct the marked process (T,D)
Tp = ep[,1];
Dp = ep[,2];


### Error analysis for likelihood filter ###

n_seq = c(10, 50, 100, 250, 500, 750, 1000, 2500, 5000, 10000, 20000, 50000, 100000)
loglik = rep(0, times = length(n_seq))
rmse = rep(0, times = length(n_seq))
rmse_rel = rep(0, times = length(n_seq))
compeffort = rep(0, times = length(n_seq))

for (ii in 1:length(n_seq)) {
  
  # Fix the parameter n for the state space discretization. The state space will be discretized based on n squared discretization points. 
  n = n_seq[ii]
  
  # Construct a time discretization and supplement it with the observation times of Y. Please enter the time discretazation width in T_diff. 
  T_diff = pmax(ceiling(60*(n^(-1/4))),1)
  Td = T_diff
  seq_times2 = seq(from = 0, to = tau2, by = T_diff)
  seq_times = c()
  j1 = 1;
  j2 = 1;
  while ((j1 <= nrow(cov)) | (j2 <= length(seq_times2))) {
    min_date = min(cov[j1,1], seq_times2[j2], na.rm = TRUE)
    if (!is.na(cov[j1,1]) & (min_date == cov[j1,1])) {
      seq_times = append(seq_times,cov[j1,1]);
      j1 = j1+1;
      j2 = if (!is.na(seq_times2[j2]) & (seq_times2[j2] < min_date + precision)) j2+1 else j2;
    } else if (!is.na((seq_times2[j2])) & ((min_date == seq_times2[j2]))) {
      seq_times = append(seq_times,seq_times2[j2]);
      j2 = j2+1;
      j1 = if (!is.na(cov[j1,1]) & (cov[j1,1] < min_date + precision)) j1+1 else j1;
    }
  }
  T_diff2 = seq_times[2:length(seq_times)] - seq_times[1:(length(seq_times)-1)]
  T_diff2 = sort(T_diff2)
  T_diff = rep(0, times = length(T_diff2));
  T_diff[1] = T_diff2[1];
  j3 = 1;
  for (i in 2:length(T_diff2)) {
    if (abs(T_diff2[i] - T_diff[j3]) > precision) {
      j3 = j3 + 1;
      T_diff[j3] = T_diff2[i]
    }
  }
  T_diff = T_diff[1:j3]
  N_T_diff = j3;
  
  # Fix the state space discretization
  p = pmin(0.25, 0.5/(sqrt(n)))
  UY = max(cov[,2], cov[1,2]*exp(sqrt(2*log(1/p))))
  LY = min(cov[,2], cov[1,2]*exp(-sqrt(2*log(1/p))))
  UZ = Z0*exp(sqrt(2*log(1/p)))
  LZ = Z0*exp(-sqrt(2*log(1/p)))
  J = pmax(ceiling((n^(1/4))*20), 2)

  
  # These matrices are needed to vectorize and speed up the computations
  zero_mat = matrix(0, ncol = J, nrow = J)
  start_F = matrix(1, ncol = J, nrow = J)
  temp_P = matrix(ncol = J, nrow = J)
  xi = matrix(0, ncol = J, nrow = J);
  lambda_nJ = matrix(ncol = J, nrow = J);
  jumps = matrix(ncol = J, nrow = J);
  lambda = matrix(ncol = J, nrow = J);
  P_Z_nJ = array(dim = c(J, J, N_T_diff))
  
  # State space discretizations for Y and Z given n and J
  IY = seq(from = LY, to = UY, length = J)
  IZ = seq(from = LZ, to = UZ, length = J)
  
  # Compute the standard errors
  time_start = Sys.time()
  loglik[ii] = filtered_likelihood_events(par, tau2)
  time_end = Sys.time()
  compeffort[ii] = as.double(difftime(time_end, time_start, units = "secs"))
  
  rmse[ii] = abs(loglik[ii] - trueloglik)
  rmse_rel[ii] = abs(loglik[ii] - trueloglik)/trueloglik
  
  print(c(n, J, Td, compeffort[ii], loglik[ii], rmse[ii], rmse_rel[ii]))
  if (ii == 1) {
    write(c(n, compeffort[ii], loglik[ii], rmse[ii], rmse_rel[ii]), "rmse_sa_asympopt_2.txt", ncolumns = 5, sep = "\t")
  } else {
    write(c(n, compeffort[ii], loglik[ii], rmse[ii], rmse_rel[ii]), "rmse_sa_asympopt_2.txt", ncolumns = 5, sep = "\t", append = TRUE)
  }
}
