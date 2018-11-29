##############################################################################################################
### Implementation of all algorithms of the paper "Filtered Likelihood for Point Processes" by Kay Giesecke and Gustavo Schwenkler. The copyright to these codes remains with the authors. ###

### PLEASE DO NOT CHANGE THESE CODES ###

##############################################################################################################

library("stats")
library("MASS")
library("nlme")
library("numDeriv")


# This function computes an approximation of the log-likelihood of Theorem 3.1 of the paper
approximate_likelihood = function(par2, tau2) {
	
	### Model parameters
	alpha = as.double(par2[1]);
	beta = as.double(par2[2]);
	gamma = as.double(par2[3]);
	kappa = as.double(par2[4]);
	gamma1 = as.double(par2[5]);
	gamma2 = as.double(par2[6]);
	gamma3 = as.double(par2[7]);
	gamma4 = as.double(par2[8]);
	gamma5 = as.double(par2[9]);
	gamma6 = as.double(par2[10]);
	
	LN = filtered_likelihood_events(par2, tau2)
	LX = conditional_likelihood_factors(par2, tau2)
			
	LN + LX
}


# This function implements Algorithm 5.2 of the paper for g = 1 and returns an approximation of the filtered likelihood of the event data
filtered_likelihood_events = function(par2, tau2) {
  
  ### Model parameters
  alpha = as.double(par2[1]);
  beta = as.double(par2[2]);
  gamma = as.double(par2[3]);
  kappa = as.double(par2[4]);
  gamma1 = as.double(par2[5]);
  gamma2 = as.double(par2[6]);
  gamma3 = as.double(par2[7]);
  gamma4 = as.double(par2[8]);
  gamma5 = as.double(par2[9]);
  gamma6 = as.double(par2[10]);
  
  Z0 = 1;
  Y0 = cov[1,2];
  
  # Value function for the point process likelihood
  xi[] = 0
  
  # If the parameter falls outside of the parameter space, return a log-likelihood of -infinity 
  if ( (alpha < 0) | (beta < 0) | (gamma < 0) | (kappa < 0) | (gamma2 < 0) | (gamma3 < 0) | (gamma4 < 0) | (gamma6 < 0) )  {
    xi = -Inf;
  } else {
    
    # Used sequence of events
    i2 = findInterval(tau2, Tp)
    T2 = Tp[1:i2]
    D2 = Dp[1:i2]
    N = length(T2)
    
    
    ### Initialization
    
    # Start the covariate and frailty location probabilities
    s1 = if (findInterval(Y0, IY) > 0) findInterval(Y0, IY) else findInterval(Y0, IY) + 1;
    s1 = if ((s1 < J) & (IY[s1] < Y0)) s1 + 1 else s1 
    s2 = if (findInterval(Z0, IZ) > 0) findInterval(Z0, IZ) else findInterval(Z0, IZ) + 1;
    s2 = if ((s2 < J) & (IZ[s2] < Z0)) s2 + 1 else s2
    
    # Intensity evaluated over the grid when no jump occurs
    for (i in 1:J) {
      lambda_nJ[i,] = alpha*IY[i] + beta*IZ;
    }
    
    # Value vector
    xi[s1, s2] = 1;
    
    # Number of time steps to go
    N3 = findInterval(tau2, seq_times)
    
    # Initialize Hawkes event sums
    int_jumps_lambda = 0		
    
    ### Iteration
    
    i2 = 1;
    i3 = 1;
    
    # These will be used to rescale the log-likelihood if it becomes too small or too large for numerical tractability
    number_divided = 0;
    number_divided2 = 0;
    number_divided_LX = 0;
    number_divided2_LX = 0;
    
    # Transition probabilities for covariate and frailty when no jumps are observed
    for (j in 1:N_T_diff) {
      P_Z_nJ[,,j] = trans_prob_Z(IZ, IZ, T_diff[j]/dpy, gamma5, gamma6, 0, 0, 0);
    }
    
    for (i in 2:N3) {
      
      t2 = seq_times[i];
      t1 = seq_times[i-1];
      
      diff_t = t2-t1;
      
      # Initialize temporary value matrix
      temp_F = start_F;
      
      
      ### Jumps
      i2 = if (findInterval(t1, T2) == 0) 1 else findInterval(t1, T2);
      i2 = if (T2[i2] < t1) i2 = i2+1 else i2;
      i3 = findInterval(t2, T2)
      
      jumps[] = 1;
      if (i3 >= i2) {				
        # Define variable F
        for (i5 in i2:i3) {
          temp = if (i2 == 1) 0 else exp(-kappa*(T2[i5]-T2[1:(i5-1)])/dpy);
          temp = gamma*sum(temp) 
          temp2 = if ((i5 == i2) & (i2 < i3)) exp(-(lambda_nJ + temp)*(T2[i5]-t1)/dpy) else if ((i5 == i2) & (i2 == i3)) exp(-(lambda_nJ + temp)*(t2-t1)/dpy) else if (i5 == i3) exp(-(lambda_nJ + temp)*(t2-T2[i5])/dpy) else exp(-(lambda_nJ + temp)*(T2[i5+1]-T2[i5])/dpy);
          jumps = jumps*(lambda_nJ + temp)*temp2
        }			
      } else {
        temp = if (i2 == 1) 0 else exp(-kappa*(t2-T2[1:(i2-1)])/dpy);
        temp = gamma*sum(temp) 
        temp2 = exp(-(lambda_nJ + temp)*(t2-t1)/dpy);
        jumps = jumps*temp2
      }
      temp_F = temp_F*jumps
      
      # Transition matrix for Z
      what_P = findInterval(diff_t, T_diff)
      tr_prob_Z = P_Z_nJ[,,what_P];
      
      
      # Transition matrix for Y
      
      pre_i4 = findInterval(t1, cov[,1]);
      post_i4 = findInterval(t2, cov[,1])
      post_i4 = if ((t2 - cov[post_i4,1] >= precision) | (pre_i4 == post_i4)) post_i4 = post_i4 + 1 else post_i4;
      
      i6 = findInterval(cov[post_i4,1], T2)
      i7 = findInterval(t1,T2);
      i5 = findInterval(cov[pre_i4,1], T2)
      
      if ((post_i4 <= nrow(cov)) & (abs(cov[post_i4,1] - t2) < precision)) { 
        # If value of covariate is observed, then deterministic transition to one of the states
        value = if (findInterval(cov[post_i4,2], IY) == 0) 1 else findInterval(cov[post_i4,2], IY);
        value = if ((value < J) & (IY[value] < cov[post_i4,2])) value + 1 else value;
        tr_prob_Y = zero_mat;
        tr_prob_Y[value,] = rep(1, times = J)
        
      } else if (post_i4 > nrow(cov)) {
        
        if (i7 > 0) {
          jumps_Y1 = exp(-gamma4*(t1-T2[1:i7])/dpy)*D2[1:i7]
          jumps_Y1 = sum(jumps_Y1)	
        } else {
          jumps_Y1 = 0	
        }
        
        if (i3 > 0) {
          jumps_Y2 = exp(-gamma4*(t2-T2[1:i3])/dpy)*D2[1:i3]
          jumps_Y2 = sum(jumps_Y2)	
        } else {
          jumps_Y2 = 0	
        }
        
        obs_cov_0 = if (t1 - cov[pre_i4,1] < precision) cov[pre_i4,2] else IY;
        
        tr_prob_Y = trans_prob_Z(IY, obs_cov_0, diff_t/dpy, gamma1, gamma2, gamma3, jumps_Y1, jumps_Y2)
        
      } else {
        # Otherwise, calculate the proper transition matrix considering the partial observability of the covariates
        
        if (i7 > 0) {
          jumps_Y1 = exp(-gamma4*(t1-T2[1:i7])/dpy)
          jumps_Y1 = sum(jumps_Y1)	
        } else {
          jumps_Y1 = 0	
        }
        
        if (i3 > 0) {
          jumps_Y2 = exp(-gamma4*(t2-T2[1:i3])/dpy)
          jumps_Y2 = sum(jumps_Y2)	
        } else {
          jumps_Y2 = 0	
        }
        
        if (i6 > 0) {
          jumps_Y3 = exp(-gamma4*(cov[post_i4,1]-T2[1:i6])/dpy)
          jumps_Y3 = sum(jumps_Y3)	
        } else {
          jumps_Y3 = 0	
        }
        
        obs_cov_0 = if (t1 - cov[pre_i4,1] < precision) cov[pre_i4,2] else IY;
        
        # Conditional normal transition distribution based on observations of the covariates
        tr_prob_Y = trans_prob_Y2(IY, obs_cov_0, cov[post_i4,2], t1, t2, cov[post_i4,1], gamma1, gamma2, gamma3, jumps_Y1, jumps_Y2, jumps_Y3);
      }
      
      # Update value vector
      temp_P = xi %*% t(tr_prob_Z);
      temp_P = tr_prob_Y %*% temp_P;
      xi = temp_P * temp_F;	
      
      
      if (is.na(sum(xi))) {
        xi = zero_mat;
        break;
      }
      
      if (sum(xi) == Inf) {
        break;
      }
      
      if (sum(xi) > 1e+150) {
        xi = xi/1e+150;
        number_divided = number_divided + 1;
      }
      
      if ((sum(xi) < 1e-50) & (sum(xi) > 0)) {
        xi = xi/1e-50;
        number_divided2 = number_divided2 + 1;
      }
      
      if (sum(xi) == 0) {
        break;	
      } 
      
    }
    
    ### Termination
    jumps = sum(exp(-kappa*(tau2 - T2)/dpy))
    int_jumps_lambda = gamma*(length(T2) - jumps)/kappa
    
    ### For log-likelihood
    xi = if (sum(xi) > 0 ) log(sum(xi)) + log(10)*(number_divided*150 - number_divided2*50) - int_jumps_lambda else -Inf;
    
  }
  
  ### Return end value	
  if (is.na(xi)) {
    xi = -Inf;
  }
  
  xi		
}


# This function computes the conditional log-likelihood of the observed factor data given the event data. This corresponds to L_F^* in Theorem 3.1
conditional_likelihood_factors = function(par3, tau3) {
	
	### Model parameters
  alpha = as.double(par3[1]);
  beta = as.double(par3[2]);
  gamma = as.double(par3[3]);
  kappa = as.double(par3[4]);
  gamma1 = as.double(par3[5]);
  gamma2 = as.double(par3[6]);
  gamma3 = as.double(par3[7]);
  gamma4 = as.double(par3[8]);
  gamma5 = as.double(par3[9]);
  gamma6 = as.double(par3[10]);
	
	X0 = cov[1,2];
	
	# Likelihood function for observed covariates and frailty
	LX = 0;
	
	if ((gamma2 < 0) | (gamma3 < 0) | (gamma4 < 0))  {
		LX = -Inf;
	} else {
		
		# Used sequence of events
		i2 = findInterval(tau3, Tp)
		T2 = Tp[1:i2]
		D2 = Dp[1:i2]
		N = length(T2)
	
		jumps_Y = 0

		number_divided_LX = 0;
		number_divided2_LX = 0;
		
		N3 = findInterval(tau3, cov[,1])
				
		
		for (i in 2:N3) {
					
			t2 = cov[i,1];
			t1 = cov[i-1,1];
			diff_t = t2-t1;
			
			i6 = findInterval(t2, T2)
			i5 = findInterval(t1, T2)
			
			if (i5 > 0) {
				jumps_Y0 = exp(-gamma4*(t1-T2[1:i5])/dpy)*D2[1:i5]
				jumps_Y0 = sum(jumps_Y0)	
			} else {
				jumps_Y0 = 0	
			}
				
			if (i6 > 0) {
				jumps_Y3 = exp(-gamma4*(t2-T2[1:i6])/dpy)*D2[1:i6]
				jumps_Y3 = sum(jumps_Y3)	
			} else {
				jumps_Y3 = 0	
			}

			
			tildeX1 = (cov[i-1,2] - gamma3*jumps_Y0)
			tildeX1 = if ((tildeX1 > -1e-5) & (tildeX1 < 1e-10)) 1e-10 else tildeX1;
			
			tildeX2 = (cov[i,2] - gamma3*jumps_Y3)
			tildeX2 = if ((tildeX2 > -1e-5) & (tildeX2 < 1e-10)) 1e-10 else tildeX2;
			
			
			if ((tildeX1 <= 0) | (tildeX2 <= 0)) {
				LX = -Inf 
				break;
			} else {
				mu_1 = tildeX2 / tildeX1;
				mu_1 = if (is.na(mu_1)) -Inf else mu_1

				LX = LX - (((log(mu_1) - gamma1*(t2 - t1)/dpy)^2)/(2*(gamma2^2)*(t2-t1)/dpy)) - log(tildeX2*gamma2*sqrt((t2-t1)/dpy))
			} 
			
#		print(c(i, LX, tildeX1, tildeX2, mu_1))		
		}
	}
	
	### Return end value	

	if (is.na(LX)) {
		LX = -Inf;
		print("NA LX end")
	}
	
#	print(c(LX, par3));
	LX	
}


# This function implements Algorithm 5.2 of the paper for g(V, tau) = Lambda(V; theta) and returns an approximation of the filtered intensity path of the marked point process
filtered_intensity = function(par2, tau2) {
	
	### Model parameters
	alpha = as.double(par2[1]);
	beta = as.double(par2[2]);
	gamma = as.double(par2[3]);
	kappa = as.double(par2[4]);
	gamma1 = as.double(par2[5]);
	gamma2 = as.double(par2[6]);
	gamma3 = as.double(par2[7]);
	gamma4 = as.double(par2[8]);
	gamma5 = as.double(par2[9]);
	gamma6 = as.double(par2[10]);
	
	# Value function for the point process likelihood
	xi = matrix(0, ncol = J, nrow = J);
	num_xi = xi

	# Used sequence of events
	i2 = findInterval(tau2, Tp)
	T2 = Tp[1:i2]
	D2 = Dp[1:i2]
	N = length(T2)

	### Initialization
	
	# Start the covariate and frailty location probabilities
	s1 = if (findInterval(Y0, IY) > 0) findInterval(Y0, IY) else findInterval(Y0, IY) + 1;
	s1 = if ((s1 < J) & (IY[s1] < Y0)) s1 + 1 else s1 
	s2 = if (findInterval(Z0, IZ) > 0) findInterval(Z0, IZ) else findInterval(Z0, IZ) + 1;
	s2 = if ((s2 < J) & (IZ[s2] < Z0)) s2 + 1 else s2
	
	# Intensity evaluated over the grid when no jump occurs
	for (i in 1:J) {
	  lambda_nJ[i,] = alpha*IY[i] + beta*IZ;
	}
	
	# Value vector
	xi[s1, s2] = 1;
		
	# Number of time steps to go
	N3 = findInterval(tau2, seq_times)
	fi = rep(0, times = N3)
		
	# Initialize Hawkes event sums
	int_jumps_lambda = 0		
	
	### Iteration
	
	i2 = 1;
	i3 = 1;
	
	# These will be used to rescale the log-likelihood if it becomes too small or too large for numerical tractability
	number_divided = 0;
	number_divided2 = 0;
	number_divided_LX = 0;
	number_divided2_LX = 0;
		
	# Transition probabilities for covariate and frailty when no jumps are observed
	for (j in 1:N_T_diff) {
		P_Z_nJ[,,j] = trans_prob_Z2(IZ, IZ, T_diff[j]/dpy, gamma5, gamma6, 0, 0, 0);
	}
		
	# Initialization
	fi[1] = alpha*Y0 + beta*Z0;	
		
	for (i in 2:N3) {
					
	  t2 = seq_times[i];
	  t1 = seq_times[i-1];
	  
	  diff_t = t2-t1;
	  
	  # Initialize temporary value matrix
	  temp_F = start_F;
	  
	  
	  ### Jumps
	  i2 = if (findInterval(t1, T2) == 0) 1 else findInterval(t1, T2);
	  i2 = if (T2[i2] < t1) i2 = i2+1 else i2;
	  i3 = findInterval(t2, T2)
	  
	  temp_i = if ((i3 > 0) && (T2[i3] == t2)) i3-1 else i3;
	  jumps_lambda = if (temp_i == 1) 0 else sum(gamma*exp(-kappa*(t2-T2[1:temp_i])/dpy));
	  
	  jumps[] = 1;
	  if (i3 >= i2) {				
	    # Define variable F
	    for (i5 in i2:i3) {
	      temp = if (i2 == 1) 0 else exp(-kappa*(T2[i5]-T2[1:(i5-1)])/dpy);
	      temp = gamma*sum(temp) 
	      temp2 = if ((i5 == i2) & (i2 < i3)) exp(-(lambda_nJ + temp)*(T2[i5]-t1)/dpy) else if ((i5 == i2) & (i2 == i3)) exp(-(lambda_nJ + temp)*(t2-t1)/dpy) else if (i5 == i3) exp(-(lambda_nJ + temp)*(t2-T2[i5])/dpy) else exp(-(lambda_nJ + temp)*(T2[i5+1]-T2[i5])/dpy);
	      jumps = jumps*(lambda_nJ + temp)*temp2
	    }			
	  } else {
	    temp = if (i2 == 1) 0 else exp(-kappa*(t2-T2[1:(i2-1)])/dpy);
	    temp = gamma*sum(temp) 
	    temp2 = exp(-(lambda_nJ + temp)*(t2-t1)/dpy);
	    jumps = jumps*temp2
	  }
	  temp_F = temp_F*jumps
	  
	  # Transition matrix for Z
	  what_P = findInterval(diff_t, T_diff)
	  tr_prob_Z = P_Z_nJ[,,what_P];
	  
	  
	  # Transition matrix for Y
	  
	  pre_i4 = findInterval(t1, cov[,1]);
	  post_i4 = findInterval(t2, cov[,1])
	  post_i4 = if ((t2 - cov[post_i4,1] >= precision) | (pre_i4 == post_i4)) post_i4 = post_i4 + 1 else post_i4;
	  
	  i6 = findInterval(cov[post_i4,1], T2)
	  i7 = findInterval(t1,T2);
	  i5 = findInterval(cov[pre_i4,1], T2)
	  
	  if ((post_i4 <= nrow(cov)) & (abs(cov[post_i4,1] - t2) < precision)) { 
	    # If value of covariate is observed, then deterministic transition to one of the states
	    value = if (findInterval(cov[post_i4,2], IY) == 0) 1 else findInterval(cov[post_i4,2], IY);
	    value = if ((value < J) & (IY[value] < cov[post_i4,2])) value + 1 else value;
	    tr_prob_Y = zero_mat;
	    tr_prob_Y[value,] = rep(1, times = J)
	    
	  } else if (post_i4 > nrow(cov)) {
	    
	    if (i7 > 0) {
	      jumps_Y1 = exp(-gamma4*(t1-T2[1:i7])/dpy)*D2[1:i7]
	      jumps_Y1 = sum(jumps_Y1)	
	    } else {
	      jumps_Y1 = 0	
	    }
	    
	    if (i3 > 0) {
	      jumps_Y2 = exp(-gamma4*(t2-T2[1:i3])/dpy)*D2[1:i3]
	      jumps_Y2 = sum(jumps_Y2)	
	    } else {
	      jumps_Y2 = 0	
	    }
	    
	    obs_cov_0 = if (t1 - cov[pre_i4,1] < precision) cov[pre_i4,2] else IY;
	    
	    tr_prob_Y = trans_prob_Z(IY, obs_cov_0, diff_t/dpy, gamma1, gamma2, gamma3, jumps_Y1, jumps_Y2)
	    
	  } else {
	    # Otherwise, calculate the proper transition matrix considering the partial observability of the covariates
	    
	    if (i7 > 0) {
	      jumps_Y1 = exp(-gamma4*(t1-T2[1:i7])/dpy)
	      jumps_Y1 = sum(jumps_Y1)	
	    } else {
	      jumps_Y1 = 0	
	    }
	    
	    if (i3 > 0) {
	      jumps_Y2 = exp(-gamma4*(t2-T2[1:i3])/dpy)
	      jumps_Y2 = sum(jumps_Y2)	
	    } else {
	      jumps_Y2 = 0	
	    }
	    
	    if (i6 > 0) {
	      jumps_Y3 = exp(-gamma4*(cov[post_i4,1]-T2[1:i6])/dpy)
	      jumps_Y3 = sum(jumps_Y3)	
	    } else {
	      jumps_Y3 = 0	
	    }
	    
	    obs_cov_0 = if (t1 - cov[pre_i4,1] < precision) cov[pre_i4,2] else IY;
	    
	    # Conditional normal transition distribution based on observations of the covariates
	    tr_prob_Y = trans_prob_Y2(IY, obs_cov_0, cov[post_i4,2], t1, t2, cov[post_i4,1], gamma1, gamma2, gamma3, jumps_Y1, jumps_Y2, jumps_Y3);
	  }
	  
	  # Update value vector
	  temp_P = xi %*% t(tr_prob_Z);
	  temp_P = tr_prob_Y %*% temp_P;
	  xi = temp_P * temp_F;	
	  
		num_xi = (lambda_nJ + jumps_lambda)*xi;
			
		fi[i] = sum(num_xi) / sum(xi)
			
			
		if (is.na(sum(xi))) {
			xi = zero_mat;
			print(c("NA xi", i));
			break;
		}
			
		if (sum(xi) == Inf) {
			print(c("xi == Inf", i));
			break;
		}
			
		if (sum(xi) > 1e+150) {
			xi = xi/1e+150;
			number_divided = number_divided + 1;
		}
			
		if ((sum(xi) < 1e-50) & (sum(xi) > 0)) {
			xi = xi/1e-50;
			number_divided2 = number_divided2 + 1;
		}
		
		if (sum(xi) == 0) {
			print(c("sum(xi) = 0", i));
			break;	
		} 
	
	}
		
	fi
		
}

# This function implements Algorithm 5.2 of the paper for g(V, tau) = Y for V = (X, Y, Z), and returns an approximation of the filtered path of the explanatory factor X conditional on its partial observations
filtered_Y = function(par2, tau2) {
	
	### Model parameters
  alpha = as.double(par2[1]);
  beta = as.double(par2[2]);
  gamma = as.double(par2[3]);
  kappa = as.double(par2[4]);
  gamma1 = as.double(par2[5]);
  gamma2 = as.double(par2[6]);
  gamma3 = as.double(par2[7]);
  gamma4 = as.double(par2[8]);
  gamma5 = as.double(par2[9]);
  gamma6 = as.double(par2[10]);

	Z0 = 1;
	Y0 = cov[1,2];
	
	# Value function for the point process likelihood
	xi = matrix(0, ncol = J, nrow = J);
	num_xi = xi
	
	
	# Used sequence of events
	i2 = findInterval(tau2, Tp)
	T2 = Tp[1:i2]
	D2 = Dp[1:i2]
	N = length(T2)
	
	
	### Initialization
		
	# Start the covariate and frailty location probabilities
	s1 = if (findInterval(Y0, IY) > 0) findInterval(Y0, IY) else findInterval(Y0, IY) + 1;
	s1 = if ((s1 < J) & (IY[s1] < Y0)) s1 + 1 else s1 
	s2 = if (findInterval(Z0, IZ) > 0) findInterval(Z0, IZ) else findInterval(Z0, IZ) + 1;
	s2 = if ((s2 < J) & (IZ[s2] < Z0)) s2 + 1 else s2
		
	# Intensity evaluated over the grid when no jump occurs
	lambda_nJ = matrix(ncol = J, nrow = J);
	for (i in 1:J) {
		lambda_nJ[i,] = alpha*IY[i] + beta*IZ;
	}
		
	# Value vector
	xi[s1, s2] = 1;
		
	# Number of time steps to go
	N3 = findInterval(tau2, seq_times)
	fi = rep(0, times = N3)
			
	# Initialize Hawkes event sums
	int_jumps_lambda = 0		
	
	### Iteration
	
	i2 = 1;
	i3 = 1;
	
	# These will be used to rescale the log-likelihood if it becomes too small or too large for numerical tractability
	number_divided = 0;
	number_divided2 = 0;
	number_divided_LX = 0;
	number_divided2_LX = 0;
		
	# Transition probabilities for covariate and frailty when no jumps are observed
	for (j in 1:N_T_diff) {
		P_Z_nJ[,,j] = trans_prob_Z2(IZ, IZ, T_diff[j]/dpy, gamma5, gamma6, 0, 0, 0);
	}
		
		
	# Initialization
	fi[1] = Y0;
		
	Y_mat = matrix(IY, ncol = J, nrow = J, byrow = FALSE)
				
	for (i in 2:N3) {
					
		t2 = seq_times[i];
		t1 = seq_times[i-1];
		
		diff_t = t2-t1;
			
		# Initialize temporary value matrix
		temp_F = start_F;
	
	
		### Jumps
		i2 = if (findInterval(t1, T2) == 0) 1 else findInterval(t1, T2);
		i2 = if (T2[i2] < t1) i2 = i2+1 else i2;
		i3 = findInterval(t2, T2)
		
		if (i3 >= i2) {				
			# Update temporary value matrix to include jumps
			for (i5 in i2:i3) {
				temp = if (i2 == 1) 0 else exp(-kappa*(T2[i5]-T2[1:(i5-1)])/dpy)*D2[i5-1];
				temp = gamma*sum(temp)
				temp_F = temp_F*(lambda_nJ + temp)
			}			
		} 
			
		# Transition matrix for frailty
		what_P = findInterval(diff_t, T_diff)
		tr_prob_Z = P_Z_nJ[,,what_P];
			
			
		# Transition matrix for covariate

		pre_i4 = findInterval(t1, cov[,1]);
		post_i4 = findInterval(t2, cov[,1])
		post_i4 = if ((t2 - cov[post_i4,1] >= precision) | (pre_i4 == post_i4)) post_i4 = post_i4 + 1 else post_i4;
			
		i6 = findInterval(cov[post_i4,1], T2)
		i7 = findInterval(t1,T2);
		i5 = findInterval(cov[pre_i4,1], T2)
			
		if ((post_i4 <= nrow(cov)) & (abs(cov[post_i4,1] - t2) < precision)) { 
			# If value of covariate is observed, then deterministic transition to one of the states
			value = if (findInterval(cov[post_i4,2], IY) == 0) 1 else findInterval(cov[post_i4,2], IY);
			value = if ((value < J) & (IY[value] < cov[post_i4,2])) value + 1 else value;
			tr_prob_Y = zero_mat;
			tr_prob_Y[value,] = rep(1, times = J)
			
		} else if (post_i4 > nrow(cov)) {
				
			if (i7 > 0) {
				jumps_Y1 = exp(-gamma4*(t1-T2[1:i7])/dpy)*D2[1:i7]
				jumps_Y1 = sum(jumps_Y1)	
			} else {
				jumps_Y1 = 0	
			}
				
			if (i3 > 0) {
				jumps_Y2 = exp(-gamma4*(t2-T2[1:i3])/dpy)*D2[1:i3]
				jumps_Y2 = sum(jumps_Y2)	
			} else {
				jumps_Y2 = 0	
			}
		
			obs_cov_0 = if (t1 - cov[pre_i4,1] < precision) cov[pre_i4,2] else IY;
				
			tr_prob_Y = trans_prob_Z(IY, obs_cov_0, diff_t/dpy, gamma1, gamma2, gamma3, jumps_Y1, jumps_Y2)
				
		} else {
		# Otherwise, calculate the proper transition matrix considering the partial observability of the covariates
			
			if (i7 > 0) {
				jumps_Y1 = exp(-gamma4*(t1-T2[1:i7])/dpy)*D2[1:i7]
				jumps_Y1 = sum(jumps_Y1)	
			} else {
				jumps_Y1 = 0	
			}
			
			if (i3 > 0) {
				jumps_Y2 = exp(-gamma4*(t2-T2[1:i3])/dpy)*D2[1:i3]
				jumps_Y2 = sum(jumps_Y2)	
			} else {
				jumps_Y2 = 0	
			}
				
			if (i6 > 0) {
				jumps_Y3 = exp(-gamma4*(cov[post_i4,1]-T2[1:i6])/dpy)*D2[1:i6]
				jumps_Y3 = sum(jumps_Y3)	
			} else {
				jumps_Y3 = 0	
			}
				
			obs_cov_0 = if (t1 - cov[pre_i4,1] < precision) cov[pre_i4,2] else IY;
				
			# Conditional normal transition distribution based on observations of the covariates
			tr_prob_Y = trans_prob_Y2(IY, obs_cov_0, cov[post_i4,2], t1, t2, cov[post_i4,1], gamma1, gamma2, gamma3, jumps_Y1, jumps_Y2, jumps_Y3);
		}
			
		# Update F vector
		jumps_lambda = if (i3 > 0) exp(-kappa*(t2-T2[1:i3])/dpy)*D2[1:i3] else 0
		jumps_lambda = sum(jumps_lambda)
		add_int_jumps_lambda = ((i3 - jumps_lambda)/kappa) - int_jumps_lambda;
		temp_F = temp_F*exp(-(diff_t*lambda_nJ/dpy) - gamma*add_int_jumps_lambda);
			
		# Update integral of jumps term
		int_jumps_lambda = int_jumps_lambda + add_int_jumps_lambda;
			
		# Update value vector
		temp_zeta = xi %*% t(tr_prob_Z);
		xi = tr_prob_Y %*% temp_zeta;
		xi = xi * temp_F;	
			
		num_xi = Y_mat*xi;
			
		fi[i] = sum(num_xi) / sum(xi)
			
		
		if (is.na(sum(xi))) {
		xi = zero_mat;
		print(c("NA xi", i));
			break;
		}
			
		if (sum(xi) == Inf) {
			break;
		}
			
		if (sum(xi) > 1e+150) {
			xi = xi/1e+150;
			number_divided = number_divided + 1;
		}
			
		if ((sum(xi) < 1e-50) & (sum(xi) > 0)) {
			xi = xi/1e-50;
			number_divided2 = number_divided2 + 1;
		}
			
		if (sum(xi) == 0) {
			break;	
		} 	
	}
		
	fi
		
}


# This function implements Algorithm 5.2 of the paper for g(V, tau) = Z for V = (X, Y, Z), and returns an approximation of the filtered path of the explanatory factor X conditional on its partial observations
filtered_Z = function(par2, tau2) {
	
	### Model parameters
  alpha = as.double(par2[1]);
  beta = as.double(par2[2]);
  gamma = as.double(par2[3]);
  kappa = as.double(par2[4]);
  gamma1 = as.double(par2[5]);
  gamma2 = as.double(par2[6]);
  gamma3 = as.double(par2[7]);
  gamma4 = as.double(par2[8]);
  gamma5 = as.double(par2[9]);
  gamma6 = as.double(par2[10]);

	Z0 = 1;
	Y0 = cov[1,2];
	
	# Value function for the point process likelihood
	xi = matrix(0, ncol = J, nrow = J);
	num_xi = xi
	

	# Used sequence of events
	i2 = findInterval(tau2, Tp)
	T2 = Tp[1:i2]
	D2 = Dp[1:i2]
	N = length(T2)
	
	
	### Initialization
		
	# Start the covariate and frailty location probabilities
	s1 = if (findInterval(Y0, IY) > 0) findInterval(Y0, IY) else findInterval(Y0, IY) + 1;
	s1 = if ((s1 < J) & (IY[s1] < Y0)) s1 + 1 else s1 
	s2 = if (findInterval(Z0, IZ) > 0) findInterval(Z0, IZ) else findInterval(Z0, IZ) + 1;
	s2 = if ((s2 < J) & (IZ[s2] < Z0)) s2 + 1 else s2
		
	# Intensity evaluated over the grid when no jump occurs
	lambda_nJ = matrix(ncol = J, nrow = J);
	for (i in 1:J) {
		lambda_nJ[i,] = alpha*IY[i] + beta*IZ;
	}
		
	# Value vector
	xi[s1, s2] = 1;
		
	# Number of time steps to go
	N3 = findInterval(tau2, seq_times)
	fi = rep(0, times = N3)
		
	# Initialize Hawkes event sums
	int_jumps_lambda = 0		
	
	### Iteration
	
	i2 = 1;
	i3 = 1;
	
	# These will be used to rescale the log-likelihood if it becomes too small or too large for numerical tractability
	number_divided = 0;
	number_divided2 = 0;
	number_divided_LX = 0;
	number_divided2_LX = 0;
		
	# Transition probabilities for covariate and frailty when no jumps are observed
	for (j in 1:N_T_diff) {
		P_Z_nJ[,,j] = trans_prob_Z2(IZ, IZ, T_diff[j]/dpy, gamma5, gamma6, 0, 0, 0);
	}
		
	# Initialization
	fi[1] = Z0;
		
	Y_mat = matrix(IZ, ncol = J, nrow = J, byrow = TRUE)
		
	PZ = rep(0, times = J)
	PZ[s2] = 1
		
		
	for (i in 2:N3) {
					
		t2 = seq_times[i];
		t1 = seq_times[i-1];
		
		diff_t = t2-t1;
			
		# Initialize temporary value matrix
		temp_F = start_F;
	
		### Jumps
		i2 = if (findInterval(t1, T2) == 0) 1 else findInterval(t1, T2);
		i2 = if (T2[i2] < t1) i2 = i2+1 else i2;
		i3 = findInterval(t2, T2)
			
		if (i3 >= i2) {				
			# Update temporary value matrix to include jumps
			for (i5 in i2:i3) {
				temp = if (i2 == 1) 0 else exp(-kappa*(T2[i5]-T2[1:(i5-1)])/dpy)*D2[i5-1];
				temp = gamma*sum(temp)
				temp_F = temp_F*(lambda_nJ + temp)
			}			
		} 
			
		# Transition matrix for frailty
		what_P = findInterval(diff_t, T_diff)
		tr_prob_Z = P_Z_nJ[,,what_P];
			
		PZ = tr_prob_Z %*% PZ
			
			
		# Transition matrix for covariate

		pre_i4 = findInterval(t1, cov[,1]);
		post_i4 = findInterval(t2, cov[,1])
		post_i4 = if ((t2 - cov[post_i4,1] >= precision) | (pre_i4 == post_i4)) post_i4 = post_i4 + 1 else post_i4;
			
		i6 = findInterval(cov[post_i4,1], T2)
		i7 = findInterval(t1,T2);
		i5 = findInterval(cov[pre_i4,1], T2)
			
		if ((post_i4 <= nrow(cov)) & (abs(cov[post_i4,1] - t2) < precision)) { 
			# If value of covariate is observed, then deterministic transition to one of the states
			value = if (findInterval(cov[post_i4,2], IY) == 0) 1 else findInterval(cov[post_i4,2], IY);
			value = if ((value < J) & (IY[value] < cov[post_i4,2])) value + 1 else value;
			tr_prob_Y = zero_mat;
			tr_prob_Y[value,] = rep(1, times = J)
				
		} else if (post_i4 > nrow(cov)) {
				
			if (i7 > 0) {
				jumps_Y1 = exp(-gamma4*(t1-T2[1:i7])/dpy)*D2[1:i7]
				jumps_Y1 = sum(jumps_Y1)	
			} else {
				jumps_Y1 = 0	
			}
				
			if (i3 > 0) {
				jumps_Y2 = exp(-gamma4*(t2-T2[1:i3])/dpy)*D2[1:i3]
				jumps_Y2 = sum(jumps_Y2)	
			} else {
				jumps_Y2 = 0	
			}
		
			obs_cov_0 = if (t1 - cov[pre_i4,1] < precision) cov[pre_i4,2] else IY;
				
			tr_prob_Y = trans_prob_Z(IY, obs_cov_0, diff_t/dpy, gamma1, gamma2, gamma3, jumps_Y1, jumps_Y2)
				
		} else {
		# Otherwise, calculate the proper transition matrix considering the partial observability of the covariates
			
			if (i7 > 0) {
				jumps_Y1 = exp(-gamma4*(t1-T2[1:i7])/dpy)*D2[1:i7]
				jumps_Y1 = sum(jumps_Y1)	
			} else {
				jumps_Y1 = 0	
			}
			
			if (i3 > 0) {
				jumps_Y2 = exp(-gamma4*(t2-T2[1:i3])/dpy)*D2[1:i3]
				jumps_Y2 = sum(jumps_Y2)	
			} else {
				jumps_Y2 = 0	
			}
				
			if (i6 > 0) {
				jumps_Y3 = exp(-gamma4*(cov[post_i4,1]-T2[1:i6])/dpy)*D2[1:i6]
				jumps_Y3 = sum(jumps_Y3)	
			} else {
				jumps_Y3 = 0	
			}
				
			obs_cov_0 = if (t1 - cov[pre_i4,1] < precision) cov[pre_i4,2] else IY;
				
			# Conditional normal transition distribution based on observations of the covariates
			tr_prob_Y = trans_prob_Y2(IY, obs_cov_0, cov[post_i4,2], t1, t2, cov[post_i4,1], gamma1, gamma2, gamma3, jumps_Y1, jumps_Y2, jumps_Y3);
		}
			
		# Update F vector
		jumps_lambda = if (i3 > 0) exp(-kappa*(t2-T2[1:i3])/dpy)*D2[1:i3] else 0
		jumps_lambda = sum(jumps_lambda)
		add_int_jumps_lambda = ((i3 - jumps_lambda)/kappa) - int_jumps_lambda;
		temp_F = temp_F*exp(-(diff_t*lambda_nJ/dpy) - gamma*add_int_jumps_lambda);
			
		# Update integral of jumps term
		int_jumps_lambda = int_jumps_lambda + add_int_jumps_lambda;
			
		# Update value vector
		temp_zeta = xi %*% t(tr_prob_Z);
		xi = tr_prob_Y %*% temp_zeta;
		xi = xi * temp_F;	
			
		num_xi = Y_mat*xi;
			
		fi[i] = sum(num_xi) / sum(xi)
			
			
		if (is.na(sum(xi))) {
			xi = zero_mat;
			print(c("NA xi", i));
			break;
		}
			
		if (sum(xi) == Inf) {
			print(c("xi == Inf", i));
			break;
		}
			
		if (sum(xi) > 1e+150) {
			xi = xi/1e+150;
			number_divided = number_divided + 1;
		}
			
		if ((sum(xi) < 1e-50) & (sum(xi) > 0)) {
			xi = xi/1e-50;
			number_divided2 = number_divided2 + 1;
		}
			
		if (sum(xi) == 0) {
			print(c("sum(xi) = 0", i));
			break;	
		} 
	}
		
	fi
		
}




### These codes compute the standard errors of a parameter estimator in accordance to Proposition 4.3 of the paper ###

st_errors = function(par3, tau4) {
	
	m = findInterval(tau4, cov[,1])
	Sigma_E = Sigma_F = Sigma_EF = matrix(0, nrow = length(par3), ncol = length(par3))
	
	grad_F = nabla_log_l_F_star(par3, tau4)	
	grad_E = nabla_log_l_E_star(par3, tau4)	
	
	for (i in 1:(m-1)) {
	  
	  grad_Fi = grad_F[i,]
	  grad_Ei = grad_E[i,]
	  
	  Sigma_EF = Sigma_EF + ((grad_Ei %*% t(grad_Fi))/(m-1))
	  Sigma_E = Sigma_E + ((grad_Ei %*% t(grad_Ei))/(m-1))
	  Sigma_F = Sigma_F + ((grad_Fi %*% t(grad_Fi))/(m-1))
	}
	
	temp = solve(Sigma_E + Sigma_F, tol = (.Machine$double.eps)^2)
	temp2 = temp %*% (Sigma_EF %*% temp)
		
	Sigma_0 = (temp + 2*temp2) * (freq_zeta/dpy)
	
	sqrt(diag(Sigma_0)) / sqrt(tau4/dpy)
}

# This function computes the vectors of gradients of the factor likelihood needed to evaluate the asymptotic variance-covariance matrix
nabla_log_l_F_star = function(par2, tau2) {
  gradient_here = jacobian(log_l_F_star, x = par2, method="Richardson", tau3 = tau2)
  #gradient_here = jacobian(log_l_F_star, x = par2, method="simple", tau3 = tau2)
  gradient_here
}


# This function computes the ratios in Equation (28) of Giesecke & Schwenkler (2017)
log_l_F_star = function(par3, tau3) {
  
  ### Model parameters
  alpha = as.double(par3[1]);
  beta = as.double(par3[2]);
  gamma = as.double(par3[3]);
  kappa = as.double(par3[4]);
  gamma1 = as.double(par3[5]);
  gamma2 = as.double(par3[6]);
  gamma3 = as.double(par3[7]);
  gamma4 = as.double(par3[8]);
  gamma5 = as.double(par3[9]);
  gamma6 = as.double(par3[10]);
  
  X0 = cov[1,2];
  
  if ((gamma2 < 0) | (gamma3 < 0) | (gamma4 < 0))  {
    LX = -Inf;
  } else {
    
    # Used sequence of events
    i2 = findInterval(tau3, Tp)
    T2 = Tp[1:i2]
    D2 = Dp[1:i2]
    N = length(T2)
    
    jumps_Y = 0
    
    number_divided_LX = 0;
    number_divided2_LX = 0;
    
    N3 = findInterval(tau3, cov[,1])
    
    # Likelihood function for observed covariates and frailty
    LX = rep(0, times = N3-1);
    
    for (i in 2:N3) {
      
      t2 = cov[i,1];
      t1 = cov[i-1,1];
      diff_t = t2-t1;
      
      i6 = findInterval(t2, T2)
      i5 = findInterval(t1, T2)
      
      if (i5 > 0) {
        jumps_Y0 = exp(-gamma4*(t1-T2[1:i5])/dpy)*D2[1:i5]
        jumps_Y0 = sum(jumps_Y0)	
      } else {
        jumps_Y0 = 0	
      }
      
      if (i6 > 0) {
        jumps_Y3 = exp(-gamma4*(t2-T2[1:i6])/dpy)*D2[1:i6]
        jumps_Y3 = sum(jumps_Y3)	
      } else {
        jumps_Y3 = 0	
      }
      
      
      tildeY1 = (cov[i-1,2] - gamma3*jumps_Y0)
      tildeY1 = if (tildeY1 < 0) 1e-10 else tildeY1;
      
      if (tildeY1 <= 0) {
        LY = -Inf 
        break;
      } else {
        
        tildeY2 = (cov[i,2] - gamma3*jumps_Y3)
        tildeY2 = if (tildeY2 <= 0) 1e-10 else tildeY2
        
        mu_1 = tildeY2 / tildeY1;
        mu_1 = if (is.na(mu_1)) -Inf else mu_1
        
        LX[i-1] = dlnorm(mu_1, meanlog = gamma1*diff_t/dpy, sdlog = gamma2*sqrt(diff_t/dpy), log = TRUE)			
      } 
      
      #		print(c(i, LX, tildeX1, tildeX2, mu_1))		
    }
  }
  
  ### Return end value	
  LX
  
}

# This function computes the vectors of gradients of the event likelihood needed to evaluate the asymptotic variance-covariance matrix
nabla_log_l_E_star = function(pars, tau4) {
  gradient_here = jacobian(log_l_E_star, x = pars, method="Richardson", tau2 = tau4)
  #gradient_here = jacobian(log_l_E_star, x = pars, method="simple", tau2 = tau4)
  gradient_here
}

# This function computes the rations in Equation (27) of Giesecke & Schwenker (2017)
log_l_E_star = function(par2, tau2) {
  
  ### Model parameters
  alpha = par2[1];
  beta = par2[2];
  gamma = par2[3];
  kappa = par2[4];
  gamma1 = par2[5];
  gamma2 = par2[6];
  gamma3 = par2[7];
  gamma4 = par2[8];
  gamma5 = par2[9];
  gamma6 = par2[10];
  
  Z0 = 1;
  Y0 = cov[1,2];
  
  # Value function for the point process likelihood
  xi[];
  
  # Likelihood function for event timing likelihood
  m = findInterval(tau2, cov[,1])
  L = rep(0, times = m-1);
  
  ### Initialization
  
  # Start the covariate and frailty location probabilities
  s1 = if (findInterval(Y0, IY) > 0) findInterval(Y0, IY) else findInterval(Y0, IY) + 1;
  s1 = if ((s1 < J) & (IY[s1] < Y0)) s1 + 1 else s1 
  s2 = if (findInterval(Z0, IZ) > 0) findInterval(Z0, IZ) else findInterval(Z0, IZ) + 1;
  s2 = if ((s2 < J) & (IZ[s2] < Z0)) s2 + 1 else s2
  
  # Intensity evaluated over the grid when no jump occurs
  for (i in 1:J) {
    lambda_nJ[i,] = alpha*IY[i] + beta*IZ;
  }
  
  # Value vector
  xi[s1, s2] = 1;
  
  # Initialize Hawkes event sums
  int_jumps_lambda = 0		
  
  ### Iteration
  
  i2 = 1;
  i3 = 1;
  
  # These will be used to rescale the log-likelihood if it becomes too small or too large for numerical tractability
  number_divided = 0;
  number_divided2 = 0;
  number_divided_LX = 0;
  number_divided2_LX = 0;
  
  # Transition probabilities for covariate and frailty when no jumps are observed
  for (j in 1:N_T_diff) {
    P_Z_nJ[,,j] = trans_prob_Z(IZ, IZ, T_diff[j]/dpy, gamma5, gamma6, 0, 0, 0);
  }
  
  # Number of time steps to go
  N3 = findInterval(tau2, seq_times)
  
  for (m_here in 2:m) {
    
    t_1 = cov[m_here-1,1]
    t_2 = cov[m_here,1]
    
    # Used sequence of events
    i2 = findInterval(t_2, Tp)
    T2 = Tp[1:i2]
    D2 = Dp[1:i2]
    N2 = length(T2)
    
    i1 = findInterval(t_1, Tp)
    T1 = Tp[1:i1]
    D1 = Dp[1:i1]
    N1 = length(T1)
    
    i_start = min(which(seq_times > t_1));
    i_end = max(which(seq_times <= t_2));
    
    if (i_end >= i_start) {
      for (i in i_start:i_end) {
        
        t2 = seq_times[i];
        t1 = ifelse(i == i_start, t_1, seq_times[i-1]);
        
        diff_t = t2-t1;
        
        # Initialize temporary value matrix
        temp_F = start_F;
        
        
        ### Jumps
        i2 = if (findInterval(t1, T2) == 0) 1 else findInterval(t1, T2);
        i2 = if (T2[i2] < t1) i2 = i2+1 else i2;
        i3 = findInterval(t2, T2)
        
        jumps[] = 1;
        if (i3 >= i2) {				
          # Define variable F
          for (i5 in i2:i3) {
            temp = if (i2 == 1) 0 else exp(-kappa*(T2[i5]-T2[1:(i5-1)])/dpy);
            temp = gamma*sum(temp) 
            temp2 = if ((i5 == i2) & (i2 < i3)) exp(-(lambda_nJ + temp)*(T2[i5]-t1)/dpy) else if ((i5 == i2) & (i2 == i3)) exp(-(lambda_nJ + temp)*(t2-t1)/dpy) else if (i5 == i3) exp(-(lambda_nJ + temp)*(t2-T2[i5])/dpy) else exp(-(lambda_nJ + temp)*(T2[i5+1]-T2[i5])/dpy);
            jumps = jumps*(lambda_nJ + temp)*temp2
          }			
        } else {
          temp = if (i2 == 1) 0 else exp(-kappa*(t2-T2[1:(i2-1)])/dpy);
          temp = gamma*sum(temp) 
          temp2 = exp(-(lambda_nJ + temp)*(t2-t1)/dpy);
          jumps = jumps*temp2
        }
        temp_F = temp_F*jumps
        
        # Transition matrix for Z
        what_P = findInterval(diff_t, T_diff)
        tr_prob_Z = P_Z_nJ[,,what_P];
        
        # Transition matrix for Y
        pre_i4 = findInterval(t1, cov[,1]);
        post_i4 = findInterval(t2, cov[,1])
        post_i4 = if ((t2 - cov[post_i4,1] >= precision) | (pre_i4 == post_i4)) post_i4 = post_i4 + 1 else post_i4;
        
        i6 = findInterval(cov[post_i4,1], T2)
        i7 = findInterval(t1,T2);
        i5 = findInterval(cov[pre_i4,1], T2)
        
        if ((post_i4 <= nrow(cov)) & (abs(cov[post_i4,1] - t2) < precision)) { 
          # If value of covariate is observed, then deterministic transition to one of the states
          value = if (findInterval(cov[post_i4,2], IY) == 0) 1 else findInterval(cov[post_i4,2], IY);
          value = if ((value < J) & (IY[value] < cov[post_i4,2])) value + 1 else value;
          tr_prob_Y = zero_mat;
          tr_prob_Y[value,] = rep(1, times = J)
          
        } else if (post_i4 > nrow(cov)) {
          
          if (i7 > 0) {
            jumps_Y1 = exp(-gamma4*(t1-T2[1:i7])/dpy)*D2[1:i7]
            jumps_Y1 = sum(jumps_Y1)	
          } else {
            jumps_Y1 = 0	
          }
          
          if (i3 > 0) {
            jumps_Y2 = exp(-gamma4*(t2-T2[1:i3])/dpy)*D2[1:i3]
            jumps_Y2 = sum(jumps_Y2)	
          } else {
            jumps_Y2 = 0	
          }
          
          obs_cov_0 = if (t1 - cov[pre_i4,1] < precision) cov[pre_i4,2] else IY;
          
          tr_prob_Y = trans_prob_Z(IY, obs_cov_0, diff_t/dpy, gamma1, gamma2, gamma3, jumps_Y1, jumps_Y2)
          
        } else {
          # Otherwise, calculate the proper transition matrix considering the partial observability of the covariates
          
          if (i7 > 0) {
            jumps_Y1 = exp(-gamma4*(t1-T2[1:i7])/dpy)
            jumps_Y1 = sum(jumps_Y1)	
          } else {
            jumps_Y1 = 0	
          }
          
          if (i3 > 0) {
            jumps_Y2 = exp(-gamma4*(t2-T2[1:i3])/dpy)
            jumps_Y2 = sum(jumps_Y2)	
          } else {
            jumps_Y2 = 0	
          }
          
          if (i6 > 0) {
            jumps_Y3 = exp(-gamma4*(cov[post_i4,1]-T2[1:i6])/dpy)
            jumps_Y3 = sum(jumps_Y3)	
          } else {
            jumps_Y3 = 0	
          }
          
          obs_cov_0 = if (t1 - cov[pre_i4,1] < precision) cov[pre_i4,2] else IY;
          
          # Conditional normal transition distribution based on observations of the covariates
          tr_prob_Y = trans_prob_Y2(IY, obs_cov_0, cov[post_i4,2], t1, t2, cov[post_i4,1], gamma1, gamma2, gamma3, jumps_Y1, jumps_Y2, jumps_Y3);
        }
        
        # Update value vector
        temp_P = xi %*% t(tr_prob_Z);
        temp_P = tr_prob_Y %*% temp_P;
        xi = temp_P * temp_F;	
        
        if (sum(xi) > 1e+150) {
          xi = xi/1e+150;
          number_divided = number_divided + 1;
        }
        
        if ((sum(xi) < 1e-50) & (sum(xi) > 0)) {
          xi = xi/1e-50;
          number_divided2 = number_divided2 + 1;
        }
      }
      
      if (seq_times[i_end] < t_2) {
        
        t2 = t_2;
        t1 = seq_times[i_end];
        
        diff_t = t2-t1;
        
        # Initialize temporary value matrix
        temp_F = start_F;
        
        
        ### Jumps
        i2 = if (findInterval(t1, T2) == 0) 1 else findInterval(t1, T2);
        i2 = if (T2[i2] < t1) i2 = i2+1 else i2;
        i3 = findInterval(t2, T2)
        
        jumps[] = 1;
        if (i3 >= i2) {				
          # Define variable F
          for (i5 in i2:i3) {
            temp = if (i2 == 1) 0 else exp(-kappa*(T2[i5]-T2[1:(i5-1)])/dpy);
            temp = gamma*sum(temp) 
            temp2 = if ((i5 == i2) & (i2 < i3)) exp(-(lambda_nJ + temp)*(T2[i5]-t1)/dpy) else if ((i5 == i2) & (i2 == i3)) exp(-(lambda_nJ + temp)*(t2-t1)/dpy) else if (i5 == i3) exp(-(lambda_nJ + temp)*(t2-T2[i5])/dpy) else exp(-(lambda_nJ + temp)*(T2[i5+1]-T2[i5])/dpy);
            jumps = jumps*(lambda_nJ + temp)*temp2
          }			
        } else {
          temp = if (i2 == 1) 0 else exp(-kappa*(t2-T2[1:(i2-1)])/dpy);
          temp = gamma*sum(temp) 
          temp2 = exp(-(lambda_nJ + temp)*(t2-t1)/dpy);
          jumps = jumps*temp2
        }
        temp_F = temp_F*jumps
        
        # Transition matrix for Z
        what_P = findInterval(diff_t, T_diff)
        tr_prob_Z = P_Z_nJ[,,what_P];
        
        # Transition matrix for Y
        pre_i4 = findInterval(t1, cov[,1]);
        post_i4 = findInterval(t2, cov[,1])
        post_i4 = if ((t2 - cov[post_i4,1] >= precision) | (pre_i4 == post_i4)) post_i4 = post_i4 + 1 else post_i4;
        
        i6 = findInterval(cov[post_i4,1], T2)
        i7 = findInterval(t1,T2);
        i5 = findInterval(cov[pre_i4,1], T2)
        
        if ((post_i4 <= nrow(cov)) & (abs(cov[post_i4,1] - t2) < precision)) { 
          # If value of covariate is observed, then deterministic transition to one of the states
          value = if (findInterval(cov[post_i4,2], IY) == 0) 1 else findInterval(cov[post_i4,2], IY);
          value = if ((value < J) & (IY[value] < cov[post_i4,2])) value + 1 else value;
          tr_prob_Y = zero_mat;
          tr_prob_Y[value,] = rep(1, times = J)
          
        } else if (post_i4 > nrow(cov)) {
          
          if (i7 > 0) {
            jumps_Y1 = exp(-gamma4*(t1-T2[1:i7])/dpy)*D2[1:i7]
            jumps_Y1 = sum(jumps_Y1)	
          } else {
            jumps_Y1 = 0	
          }
          
          if (i3 > 0) {
            jumps_Y2 = exp(-gamma4*(t2-T2[1:i3])/dpy)*D2[1:i3]
            jumps_Y2 = sum(jumps_Y2)	
          } else {
            jumps_Y2 = 0	
          }
          
          obs_cov_0 = if (t1 - cov[pre_i4,1] < precision) cov[pre_i4,2] else IY;
          
          tr_prob_Y = trans_prob_Z(IY, obs_cov_0, diff_t/dpy, gamma1, gamma2, gamma3, jumps_Y1, jumps_Y2)
          
        } else {
          # Otherwise, calculate the proper transition matrix considering the partial observability of the covariates
          
          if (i7 > 0) {
            jumps_Y1 = exp(-gamma4*(t1-T2[1:i7])/dpy)
            jumps_Y1 = sum(jumps_Y1)	
          } else {
            jumps_Y1 = 0	
          }
          
          if (i3 > 0) {
            jumps_Y2 = exp(-gamma4*(t2-T2[1:i3])/dpy)
            jumps_Y2 = sum(jumps_Y2)	
          } else {
            jumps_Y2 = 0	
          }
          
          if (i6 > 0) {
            jumps_Y3 = exp(-gamma4*(cov[post_i4,1]-T2[1:i6])/dpy)
            jumps_Y3 = sum(jumps_Y3)	
          } else {
            jumps_Y3 = 0	
          }
          
          obs_cov_0 = if (t1 - cov[pre_i4,1] < precision) cov[pre_i4,2] else IY;
          
          # Conditional normal transition distribution based on observations of the covariates
          tr_prob_Y = trans_prob_Y2(IY, obs_cov_0, cov[post_i4,2], t1, t2, cov[post_i4,1], gamma1, gamma2, gamma3, jumps_Y1, jumps_Y2, jumps_Y3);
        }
        
        # Update value vector
        temp_P = xi %*% t(tr_prob_Z);
        temp_P = tr_prob_Y %*% temp_P;
        xi = temp_P * temp_F;	
        
        if (sum(xi) > 1e+150) {
          xi = xi/1e+150;
          number_divided = number_divided + 1;
        }
        
        if ((sum(xi) < 1e-50) & (sum(xi) > 0)) {
          xi = xi/1e-50;
          number_divided2 = number_divided2 + 1;
        }
      }
      
      
      
      ### Termination
      jumps = sum(exp(-kappa*(t_2 - T2)/dpy))
      int_jumps_lambda = gamma*(length(T2) - jumps)/kappa
      
      ### For log-likelihood
      denominator = ifelse(m_here > 2, numerator, 0)
      numerator = if (sum(xi) > 0 ) log(sum(xi)) + log(10)*(number_divided*150 - number_divided2*50) - int_jumps_lambda else -Inf;
      
    } else {
      
      t2 = t_2
      t1 = t_1;
      
      diff_t = t2-t1;
      
      # Initialize temporary value matrix
      temp_F = start_F;
      
      ### Jumps
      i2 = if (findInterval(t1, T2) == 0) 1 else findInterval(t1, T2);
      i2 = if (T2[i2] < t1) i2 = i2+1 else i2;
      i3 = findInterval(t2, T2)
      
      jumps[] = 1;
      if (i3 >= i2) {				
        # Define variable F
        for (i5 in i2:i3) {
          temp = if (i2 == 1) 0 else exp(-kappa*(T2[i5]-T2[1:(i5-1)])/dpy);
          temp = gamma*sum(temp) 
          temp2 = if ((i5 == i2) & (i2 < i3)) exp(-(lambda_nJ + temp)*(T2[i5]-t1)/dpy) else if ((i5 == i2) & (i2 == i3)) exp(-(lambda_nJ + temp)*(t2-t1)/dpy) else if (i5 == i3) exp(-(lambda_nJ + temp)*(t2-T2[i5])/dpy) else exp(-(lambda_nJ + temp)*(T2[i5+1]-T2[i5])/dpy);
          jumps = jumps*(lambda_nJ + temp)*temp2
        }			
      } else {
        temp = if (i2 == 1) 0 else exp(-kappa*(t2-T2[1:(i2-1)])/dpy);
        temp = gamma*sum(temp) 
        temp2 = exp(-(lambda_nJ + temp)*(t2-t1)/dpy);
        jumps = jumps*temp2
      }
      temp_F = temp_F*jumps
      
      # Transition matrix for Z
      what_P = findInterval(diff_t, T_diff)
      tr_prob_Z = P_Z_nJ[,,what_P];
      
      # Transition matrix for Y
      pre_i4 = findInterval(t1, cov[,1]);
      post_i4 = findInterval(t2, cov[,1])
      post_i4 = if ((t2 - cov[post_i4,1] >= precision) | (pre_i4 == post_i4)) post_i4 = post_i4 + 1 else post_i4;
      
      i6 = findInterval(cov[post_i4,1], T2)
      i7 = findInterval(t1,T2);
      i5 = findInterval(cov[pre_i4,1], T2)
      
      if ((post_i4 <= nrow(cov)) & (abs(cov[post_i4,1] - t2) < precision)) { 
        # If value of covariate is observed, then deterministic transition to one of the states
        value = if (findInterval(cov[post_i4,2], IY) == 0) 1 else findInterval(cov[post_i4,2], IY);
        value = if ((value < J) & (IY[value] < cov[post_i4,2])) value + 1 else value;
        tr_prob_Y = zero_mat;
        tr_prob_Y[value,] = rep(1, times = J)
        
      } else if (post_i4 > nrow(cov)) {
        
        if (i7 > 0) {
          jumps_Y1 = exp(-gamma4*(t1-T2[1:i7])/dpy)*D2[1:i7]
          jumps_Y1 = sum(jumps_Y1)	
        } else {
          jumps_Y1 = 0	
        }
        
        if (i3 > 0) {
          jumps_Y2 = exp(-gamma4*(t2-T2[1:i3])/dpy)*D2[1:i3]
          jumps_Y2 = sum(jumps_Y2)	
        } else {
          jumps_Y2 = 0	
        }
        
        obs_cov_0 = if (t1 - cov[pre_i4,1] < precision) cov[pre_i4,2] else IY;
        
        tr_prob_Y = trans_prob_Z(IY, obs_cov_0, diff_t/dpy, gamma1, gamma2, gamma3, jumps_Y1, jumps_Y2)
        
      } else {
        # Otherwise, calculate the proper transition matrix considering the partial observability of the covariates
        
        if (i7 > 0) {
          jumps_Y1 = exp(-gamma4*(t1-T2[1:i7])/dpy)
          jumps_Y1 = sum(jumps_Y1)	
        } else {
          jumps_Y1 = 0	
        }
        
        if (i3 > 0) {
          jumps_Y2 = exp(-gamma4*(t2-T2[1:i3])/dpy)
          jumps_Y2 = sum(jumps_Y2)	
        } else {
          jumps_Y2 = 0	
        }
        
        if (i6 > 0) {
          jumps_Y3 = exp(-gamma4*(cov[post_i4,1]-T2[1:i6])/dpy)
          jumps_Y3 = sum(jumps_Y3)	
        } else {
          jumps_Y3 = 0	
        }
        
        obs_cov_0 = if (t1 - cov[pre_i4,1] < precision) cov[pre_i4,2] else IY;
        
        # Conditional normal transition distribution based on observations of the covariates
        tr_prob_Y = trans_prob_Y2(IY, obs_cov_0, cov[post_i4,2], t1, t2, cov[post_i4,1], gamma1, gamma2, gamma3, jumps_Y1, jumps_Y2, jumps_Y3);
      }
      
      # Update value vector
      temp_P = xi %*% t(tr_prob_Z);
      temp_P = tr_prob_Y %*% temp_P;
      xi = temp_P * temp_F;	
      
      if (sum(xi) > 1e+150) {
        xi = xi/1e+150;
        number_divided = number_divided + 1;
      }
      
      if ((sum(xi) < 1e-50) & (sum(xi) > 0)) {
        xi = xi/1e-50;
        number_divided2 = number_divided2 + 1;
      }
      
      ### Termination
      jumps = sum(exp(-kappa*(t_2 - T2)/dpy))
      int_jumps_lambda = gamma*(length(T2) - jumps)/kappa
      
      ### For log-likelihood
      denominator = ifelse(m_here > 2, numerator, 0)
      numerator = if (sum(xi) > 0 ) log(sum(xi)) + log(10)*(number_divided*150 - number_divided2*50) - int_jumps_lambda else -Inf;
      
    }
    
#    print(c(m_here, numerator - denominator, numerator, denominator)) 
    L[m_here-1] = numerator - denominator;
  }
  
  L
  
}
