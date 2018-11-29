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

# Load all relevenat algorithms
source("algorithms_higherdim.R")
source("input_probabilities_higherdim.R")
require(compiler)
enableJIT(3)
library("bigmemory")
library("bigalgebra")

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

dY = 4
dZ = 1

# True likelihood
trueloglik = 8659.09

n_seq = c(10, 20, 30, 40, 50, 100, 500, 750, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 100000, 300000, 600000, 1e6)
loglik = rep(0, times = length(n_seq))
rmse = rep(0, times = length(n_seq))
rmse_rel = rep(0, times = length(n_seq))
compeffort = rep(0, times = length(n_seq))

for (ii in 1:length(n_seq)) {
  
  # Fix the parameter n for the state space discretization. The state space will be discretized based on n squared discretization points. 
  n = n_seq[ii]
  
  # Construct a time discretization and supplement it with the observation times of Y. Please enter the time discretazation width in T_diff. 
  temp = n^(-1/4)
  T_diff = ceiling(pmin(freq_zeta*0.5, pmax(200*temp,1)))
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
  
  # Fix number of discretization points per axis
  J = pmax(ceiling((n^(1/4))*20), 2)
  JY = pmax(2, floor(3*(n^(1/(4*(dY))))))
  JZ = pmax(2, floor(5*(n^(1/(4*(dZ))))))
  
  # Fix the state space discretization
  pY = 1/(0.75*JY)
  pZ = 1/(JZ)
  UY = cov[1,2]*exp(sqrt(2*log(1/pY)))/3
  LY = cov[1,2]*exp(-sqrt(2*log(1/pY)))/3
  UZ = Z0*exp(sqrt(2*log(1/pZ)))
  LZ = Z0*exp(-sqrt(2*log(1/pZ)))
  
  # State space discretizations for Y and Z given n and J
  IY = seq(from = LY, to = UY, length = JY)
  IZ = seq(from = LZ, to = UZ, length = JZ)
  
  # These matrices are needed to vectorize and speed up the computations
  zero_mat_Y = matrix(0, ncol = JY, nrow = JY)
  start_F = matrix(1, ncol = JZ^dZ, nrow = JY^dY)
  temp_P = matrix(0, ncol = JZ^dZ, nrow = JY^dY)
  xi = matrix(0, ncol = JZ^dZ, nrow = JY^dY);
  lambda_nJ = lambda_nJ_Y = matrix(1, ncol = JZ^dZ, nrow = JY^dY);
  lambda_nJ_Z = matrix(1, nrow = JZ^dZ, ncol = JY^dY);
  jumps = matrix(0, ncol = JZ^dZ, nrow = JY^dY);
  lambda = matrix(0, ncol = JZ^dZ, nrow = JY^dY);
  P_Z_nJ = array(dim = c(JZ, JZ, N_T_diff))
  tr_prob_Y = matrix(0, nrow = JY^dY, ncol = JY^dY)
  tr_prob_Z = matrix(0, nrow = JZ^dZ, ncol = JZ^dZ)
  indices_Y = as.matrix(expand.grid(rep(list(1:JY), times = dY)), dimnames = NULL)
  colnames(indices_Y) = NULL
  indices_Z = as.matrix(expand.grid(rep(list(1:JZ), times = dZ)))
  colnames(indices_Z) = NULL
  
  # Compute filter approximation
  time_start = Sys.time()
  loglik[ii] = filtered_likelihood_events_higherdim_2(par, tau2)
  time_end = Sys.time()
  compeffort[ii] = as.double(difftime(time_end, time_start, units = "secs"))
  
  rmse[ii] = abs(loglik[ii] - trueloglik)
  rmse_rel[ii] = abs(loglik[ii] - trueloglik)/trueloglik
  
  print(round(c(n, JY, JZ, Td, compeffort[ii], loglik[ii], rmse[ii], rmse_rel[ii]), 3))
  if (ii == 1) {
    write(c(n, compeffort[ii], loglik[ii], rmse[ii], rmse_rel[ii]), paste("rmse_sa_asympopt_dY_", dY, "_dZ_", dZ, ".txt", sep = ""), ncolumns = 5, sep = "\t")
  } else {
    write(c(n, compeffort[ii], loglik[ii], rmse[ii], rmse_rel[ii]), paste("rmse_sa_asympopt_dY_", dY, "_dZ_", dZ, ".txt", sep = ""), ncolumns = 5, sep = "\t", append = TRUE)
  }
  
}
