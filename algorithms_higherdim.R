##############################################################################################################
### Implementation of all algorithms of the paper "Filtered Likelihood for Point Processes" by Kay Giesecke and Gustavo Schwenkler for the higher dimensional case. The copyright to these codes remains with the authors. ###

### PLEASE DO NOT CHANGE THESE CODES ###

##############################################################################################################

library("stats")
library("MASS")
library("nlme")

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


update_temp_F = function(t1, t2, T2, lambda_nJ) {
  
  diff_t = t2-t1;
  
  # Initialize temporary value matrix
  temp_F = start_F[];

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
      temp2 = if ((i5 == i2) & (i2 < i3)) exp(-(lambda_nJ[] + temp)*(T2[i5]-t1)/dpy) else if ((i5 == i2) & (i2 == i3)) exp(-(lambda_nJ[] + temp)*(t2-t1)/dpy) else if (i5 == i3) exp(-(lambda_nJ[] + temp)*(t2-T2[i5])/dpy) else exp(-(lambda_nJ[] + temp)*(T2[i5+1]-T2[i5])/dpy);
      jumps[] = jumps[]*(lambda_nJ[] + temp)*temp2
    }			
  } else {
    temp = if (i2 == 1) 0 else exp(-kappa*(t2-T2[1:(i2-1)])/dpy);
    temp = gamma*sum(temp) 
    temp2 = exp(-(lambda_nJ[] + temp)*(t2-t1)/dpy);
    jumps[] = jumps[]*temp2
  }
  temp_F = temp_F*jumps[]
  
  temp_F
}


construct_tr_prob_Y = function(t1, t2, T2, D2) {
  
  diff_t = t2-t1;
  
  ### Jumps
  i2 = if (findInterval(t1, T2) == 0) 1 else findInterval(t1, T2);
  i2 = if (T2[i2] < t1) i2 = i2+1 else i2;
  i3 = findInterval(t2, T2)
  
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
    value = if ((value < JY) & (IY[value] < cov[post_i4,2])) value + 1 else value;
    temp_Y = zero_mat_Y;
    temp_Y[value,] = rep(1, times = JY)
    
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
    
    temp_Y = trans_prob_Z(IY, obs_cov_0, diff_t/dpy, gamma1, gamma2, gamma3, jumps_Y1, jumps_Y2)
    
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
    temp_Y = trans_prob_Y2(IY, obs_cov_0, cov[post_i4,2], t1, t2, cov[post_i4,1], gamma1, gamma2, gamma3, jumps_Y1, jumps_Y2, jumps_Y3);
  }
  
  temp_Y
}

update_xi = function(l, xi, tr_prob_Y, tr_prob_Z, temp_F) {
  probs_here[] = 1;
  for (i8 in 1:dY) {
    probs_here = probs_here*tr_prob_Y[indices[,i8],indices[l,i8]]
  }
  for (i8 in (dY+1):(dY+dZ)) {
    probs_here = probs_here*tr_prob_Z[indices[,i8],indices[l,i8]]
  }
  res = temp_F[l]*as.double(probs_here %*% xi)
  res
}


update_xi_2 = function(xi, tr_prob_Y, tr_prob_Z, temp_F) {
  hat_Q[] = 1
  for (i8 in 1:dY) {
    hat_Q[] = hat_Q[]*tr_prob_Y[indices[,i8],indices[,i8]]
  }
  for (i8 in (dY+1):(dY+dZ)) {
    hat_Q[] = hat_Q[]*tr_prob_Z[indices[,i8],indices[,i8]]
  }
  hat_Q[] = temp_F*hat_Q[]
  as.double(hat_Q[] %*% xi)
}

IY_mat_to_vec = function(index) {
  mean(IY[indices_Y[index,]])
}

IZ_mat_to_vec = function(index) {
  mean(IZ[indices_Z[index,]])
}

PZ_mat_to_vec = function(k, temp_Z) {
  tr_prob_Z[,k] = 1
  for (l in 1:dZ) {
    tr_prob_Z[,k] = tr_prob_Z[,k]*temp_Z[indices_Z[,l], indices_Z[k,l]]
  }
}

filtered_likelihood_events_higherdim_2 = function(par2, tau2) {
  
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
    sY = if (findInterval(Y0, IY) > 0) findInterval(Y0, IY) else findInterval(Y0, IY) + 1;
    sY = if ((sY < JY) & (IY[sY] < Y0)) sY + 1 else sY 
    sZ = if (findInterval(Z0, IZ) > 0) findInterval(Z0, IZ) else findInterval(Z0, IZ) + 1;
    sZ = if ((sZ < JZ) & (IZ[sZ] < Z0)) sZ + 1 else sZ
    
    # Intensity evaluated over the grid when no jump occurs
    lambda_nJ_Y[] = sapply(1:(JY^dY), IY_mat_to_vec)
    lambda_nJ_Z[] = sapply(1:(JZ^dZ), IZ_mat_to_vec)
    lambda_nJ[] = alpha*lambda_nJ_Y[] + beta*t(lambda_nJ_Z[])

    # Value vector
    pos_Y = which(apply(indices_Y, 1, function(x) all.equal(x[1:dY], rep(sY, times = dY))) == "TRUE")
    pos_Z = which(apply(indices_Z, 1, function(x) all.equal(x[1:dZ], rep(sZ, times = dZ))) == "TRUE")
    xi[pos_Y, pos_Z] = 1;
    
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
      temp_F = update_temp_F(t1, t2, T2, lambda_nJ)
      
      # Compute transition probabilities
      what_P = findInterval(diff_t, T_diff)
      temp_tr_prob_Z = P_Z_nJ[,,what_P];
      temp_tr_prob_Y = construct_tr_prob_Y(t1, t2, T2, D2);
      
      for (k in 1:(JY^dY)) {
        tr_prob_Y[,k] = 1
        for (l in 1:dY) {
          tr_prob_Y[,k] = tr_prob_Y[,k]*temp_tr_prob_Y[indices_Y[,l], indices_Y[k,l]]
        }
      }
      
      for (k in 1:(JZ^dZ)) {
        tr_prob_Z[,k] = 1
        for (l in 1:dZ) {
          tr_prob_Z[,k] = tr_prob_Z[,k]*temp_tr_prob_Z[indices_Z[,l], indices_Z[k,l]]
        }
      }
      
      #temp = sapply(1:(JY^dY), PY_mat_to_vec, temp_Y = temp_tr_prob_Y)
      #temp = sapply(1:(JZ^dZ), PZ_mat_to_vec, temp_Z = temp_tr_prob_Z)
    
      #tr_prob_Y[] = matrix(sapply(1:nrow(rows_Y), P_mat_to_vec, temp_p = temp_tr_prob_Y, indices = indices_Y, rows = rows_Y, d = dY), nrow = JY^dY, ncol = JY^dY, byrow = FALSE)
      #tr_prob_Z[] = matrix(sapply(1:nrow(rows_Z), P_mat_to_vec, temp_p = temp_tr_prob_Z, indices = indices_Z, rows = rows_Z, d = dZ), nrow = JZ^dZ, ncol = JZ^dZ, byrow = FALSE)
      
      # Update value vector
      temp_P = xi %*% t(tr_prob_Z[]);
      temp_P = tr_prob_Y[] %*% temp_P;
      xi = temp_P * temp_F;	
      
      
      if (is.na(sum(xi))) {
        xi = 0;
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
      
#      print(c(i, sum(xi)))
      
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
