##############################################################################################################
### Computation of the probabilities hat p of Section 4 of the paper "Filtered Likelihood for Point Processes" by Kay Giesecke and Gustavo Schwenkler. The copyright to these codes remains with the authors. ###
### If you wish to work with a model different than the model in Section 6 of the paper, then you need to input here the transition density of the frailty Z and the conditional transition density of the partially observed factor X ###

##############################################################################################################


# This function computes the transition density of the frailty Z. In the paper, the transition density of Z is log-normal
trans_prob_Z = function(to_vec, from_obs, difft, g1, g2, g3, jumps1, jumps2) {
	
	
	len = length(to_vec);
	result = matrix(0, nrow = len, ncol = len);
	
	Del = to_vec[2] - to_vec[1]
	to_vec2 = c(to_vec[1]-Del, to_vec)


	if (length(from_obs) > 1) {
		
		for (i in 1:len) {
			mu = max(from_obs[i]-g3*jumps1, 1e-10)
		
			new_I_X = (to_vec2-g3*jumps2)/mu;
			temp_X = plnorm(new_I_X, meanlog = (g1)*difft, sdlog = g2*sqrt(difft));
			temp_X = temp_X[2:(len+1)] - temp_X[1:len]
		
			result[,i] = temp_X
			
			
		}
		
	} else {
		
		where = if (findInterval(from_obs, to_vec) == 0) 1 else findInterval(from_obs, to_vec);
		where = if ((where < J) & (to_vec[where] < from_obs)) where + 1 else where;
		
		mu = max(to_vec[where]-g3*jumps1, 1e-10)
			
		new_I_X = (to_vec2-g3*jumps2)/mu;
		temp_X = plnorm(new_I_X, meanlog = (g1)*difft, sdlog = g2*sqrt(difft));
		temp_X = temp_X[2:(len+1)] - temp_X[1:len]
	
		result[,where] = temp_X
		
	}
	
	result	
}


# This function computes the transition density of the frailty Z. In the paper, the transition density of Z is log-normal
trans_prob_Z2 = function(to_vec, from_obs, difft, g1, g2, g3, jumps1, jumps2) {
  
  
  len = length(to_vec);
  result = matrix(0, nrow = len, ncol = len);
  
  Del = to_vec[2] - to_vec[1]
  to_vec2 = c(to_vec[1]-Del/2, to_vec+Del/2)
  
  
  if (length(from_obs) > 1) {
    
    for (i in 1:len) {
      mu = max(from_obs[i]-g3*jumps1, 1e-10)
      
      new_I_X = (to_vec2-g3*jumps2)/mu;
      temp_X = plnorm(new_I_X, meanlog = (g1)*difft, sdlog = g2*sqrt(difft));
      temp_X = temp_X[2:(len+1)] - temp_X[1:len]
      
      result[,i] = temp_X
      
      
    }
    
  } else {
    
    where = if (findInterval(from_obs, to_vec) == 0) 1 else findInterval(from_obs, to_vec);
    where = if ((where < J) & (to_vec[where] < from_obs)) where + 1 else where;
    
    mu = max(to_vec[where]-g3*jumps1, 1e-10)
    
    new_I_X = (to_vec2-g3*jumps2)/mu;
    temp_X = plnorm(new_I_X, meanlog = (g1)*difft, sdlog = g2*sqrt(difft));
    temp_X = temp_X[2:(len+1)] - temp_X[1:len]
    
    result[,where] = temp_X
    
  }
  
  result	
}


# This function computes the conditional transition density of the partially observed factor X given the partial data on X. In the paper, conditional on the partial data on X and the jumps of X, the transition density of X is log-normal
trans_prob_Y2 = function(to_vec, from_obs, next_obs, t_1, t_2, t_3, g1, g2, g3, jumps1, jumps2, jumps3) {
	
	len = length(to_vec)
	result = matrix(0, nrow = len, ncol = len);
	
	
	rho = sqrt( (t_2 - t_1) / (t_3 - t_1) )
	stdev2 = sqrt((t_2 - t_1)/dpy)
	stdev3 = sqrt((t_3 - t_1)/dpy)
	
	Del = to_vec[2] - to_vec[1]
	to_vec2 = c(to_vec[1]-Del, to_vec)
	
	if ((length(from_obs) > 1) & (length(to_vec) > 1)) {
		
		for (j in 1:J) {
			
			tildeY0 = (from_obs[j] - g3*jumps1)
			if (tildeY0 < 0) {
				
				temp_X = rep(0, times = len)
				temp_X[1] = 1
				
			} else {
				
				mu1 = (log((to_vec2 - g3*jumps2) / tildeY0) - g1*(t_2 - t_1)/dpy) / g2;	
				mu1[which(is.na(mu1))] = -Inf	
				
				land_in = findInterval(next_obs, to_vec)
				land_in = if ((land_in < J) & (to_vec[land_in] < next_obs)) land_in + 1 else land_in;
				
				mu2 = (log((to_vec[land_in] - g3*jumps3) / tildeY0) - g1*(t_3 - t_1)/dpy) / g2;
				mu2 = if (is.na(mu2)) -Inf else mu2

				stdev = stdev2*sqrt(1 - rho^2)
		
				temp_X = pnorm(mu1, mean = stdev2*rho*mu2/stdev3, sd = stdev);
				temp_X = temp_X[2:(len+1)] - temp_X[1:len]
				if (mu2 < 0) {
				  temp_X[] = 0
				  temp_X[1] = 1
				} else if (mu2 == Inf) {
				  temp_X[] = 0
				  temp_X[length(temp_X)] = 1
				} else if ((mu2 > 0) & (max(mu1) == -Inf)) {
				  temp_X[] = 1/length(temp_X)
				}
			
			} 
			
			result[,j] = temp_X
		}
		
	} else {
		
		where = if (findInterval(from_obs, to_vec) == 0) 1 else findInterval(from_obs, to_vec);
		where = if ((where < J) & (to_vec[where] < from_obs)) where + 1 else where;
		
		tildeY0 = (to_vec[where] - g3*jumps1)
		
		if (tildeY0 < 0) {
				
			temp_X = rep(0, times = len)
			temp_X[1] = 1;
				
		} else {
				
			mu1 = (log((to_vec2 - g3*jumps2) / tildeY0) - g1*(t_2 - t_1)/dpy) / g2;	
			mu1[which(is.na(mu1))] = -Inf
			
			land_in = findInterval(next_obs, to_vec)
			land_in = if ((land_in < J) & (to_vec[land_in] < next_obs)) land_in + 1 else land_in;
				
			mu2 = (log((to_vec[land_in] - g3*jumps3) / tildeY0) - g1*(t_3 - t_1)/dpy) / g2;
			mu2 = if (is.na(mu2)) -Inf else mu2

			stdev = stdev2*sqrt(1 - rho^2)
		
			temp_X = pnorm(mu1, mean = stdev2*rho*mu2/stdev3, sd = stdev);
			temp_X = temp_X[2:(len+1)] - temp_X[1:len]
			if (mu2 < 0) {
			  temp_X[] = 0
			  temp_X[1] = 1;
			} else if (mu2 == Inf) {
			  temp_X[] = 0
			  temp_X[length(temp_X)] = 1
			} else if ((mu2 > 0) & (max(mu1) == -Inf)) {
			  temp_X[] = 1/length(temp_X)
			}
			
		} 

		result[,where] = temp_X
		
	}
		
	result	
}


