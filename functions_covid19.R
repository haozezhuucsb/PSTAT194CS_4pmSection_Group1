########################
# functions for SIRDC model for U.S. counties
#
# Author: Hanmo Li & Mengyang Gu 
#
# Reference: Robust estimation of SARS-CoV-2 epidemic at US counties (https://arxiv.org/abs/2010.11514)
# COVID-19 dashboard: https://covid19-study.pstat.ucsb.edu/
########################

eliminate_abnormal_death = function(data){
  abnormal_index = which(diff(data)<0)
  if(length(abnormal_index)==0){
    return(data)
  } else{
    for (i in 1:length(abnormal_index)){
      data[abnormal_index[i]] = data[abnormal_index[i]+1]
    }
  }
  return(data)
}

data_seven_day_smoothing = function(data, default_method = T){

  if(!default_method){
    data_length = length(data)
    
    data_seven_day_average = rep(NA, data_length)
    
    data_seven_day_average[(1+3):(data_length-3)] = rollapply(data, width = 7, by = 1, FUN = mean, align = "right")
    
    data_seven_day_average[1+2] = mean(data[1:5], na.rm=T)
    data_seven_day_average[data_length-2] = mean(data[(data_length-4):(data_length)], na.rm=T)
    
    data_seven_day_average[1+1] = mean(data[1:3], na.rm=T)
    data_seven_day_average[data_length-1] = mean(data[(data_length-2):(data_length)], na.rm=T)
    
    # data_seven_day_average[1] = mean(data[1:4], na.rm = T)
    # data_seven_day_average[data_length] = mean(data[(data_length-3):data_length], na.rm=T)
    
    data_seven_day_average[1] = mean(data[1], na.rm=T)
    data_seven_day_average[data_length] = mean(data[data_length], na.rm=T)
  } else{
    data_length = length(data)
    
    data_seven_day_average = rep(NA, data_length)
    
    data_seven_day_average[(1+3):(data_length-3)] = rollapply(data, width = 7, by = 1, FUN = mean, align = "right")
    
    data_seven_day_average[1+2] = mean(data[1:6], na.rm=T)
    data_seven_day_average[data_length-2] = mean(data[(data_length-5):data_length], na.rm=T)
    
    data_seven_day_average[1+1] = mean(data[1:5], na.rm=T)
    data_seven_day_average[data_length-1] = mean(data[(data_length-4):data_length], na.rm=T)
    
    # data_seven_day_average[1] = mean(data[1:4], na.rm = T)
    # data_seven_day_average[data_length] = mean(data[(data_length-3):data_length], na.rm=T)
    
    data_seven_day_average[1] = mean(data[1:4], na.rm=T)
    data_seven_day_average[data_length] = mean(data[(data_length-3):data_length], na.rm=T)
  }
  
  return(data_seven_day_average)
}






clean_us_death = function(data){
  
  
  # data cleaning
  names(data)[names(data) == "Long_"] = "Long"
  names(data)[names(data) == "FIPS"] = "county_fips"
  data$county_fips = as.character(data$county_fips)
  # there are some county fips missing, just redistribute or delete them
  missing_fips = which(is.na(data$county_fips))
  # data$county_fips[3148] = "25007"
  data = data[-missing_fips,]
  
  k = dim(data)[1]
  n = dim(data)[2]
  

  return(data)
}

clean_JHU_data_for_map = function(confirm_data, death_data){
  
  
  # clean the county fips
  names(confirm_data)[names(confirm_data) == "Long_"] = "Long"
  names(confirm_data)[names(confirm_data) == "FIPS"] = "county_fips"
  confirm_data$county_fips = as.character(confirm_data$county_fips)
  
  names(death_data)[names(death_data) == "Long_"] = "Long"
  names(death_data)[names(death_data) == "FIPS"] = "county_fips"
  death_data$county_fips = as.character(death_data$county_fips)

  # there are some county fips missing, just delete them
  missing_fips = which(is.na(confirm_data$county_fips))
  if(length(missing_fips)>0){
    confirm_data = confirm_data[-missing_fips,]
  }
  
  missing_fips = which(is.na(death_data$county_fips))
  if(length(missing_fips)>0){
    death_data = death_data[-missing_fips,]
  }

  # delete the records with population = 0
  # impute "0" in the front of "county_fips" to make sure they are in 6-digits

  no_population_idx = c()
  for (i in 1:dim(confirm_data)[1]){
    county_fips_iter = confirm_data$county_fips[i]
    county_fips_length = str_length(county_fips_iter)
    diff = 5-county_fips_length
    if (diff<5){
      confirm_data$county_fips[i] = paste0(do.call(paste0,as.list((rep(0,diff)))), county_fips_iter)
    }
    if(death_data$Population[i]==0){
      no_population_idx = c(no_population_idx, i)
    }
  }
  if(length(no_population_idx)>0){
    confirm_data = confirm_data[-no_population_idx,]
    death_data = death_data[-no_population_idx,]
  }
  
  results = list()
  results[[1]] = confirm_data
  results[[2]] = death_data
  return(results)
}


get_output_same_time_zone = function(data_type = "death_rate", state_name, state_name_short,start_date,training_length, duration = 90, criterion_death = 10, smoothness = TRUE){
  
  state_death = us_death %>%
    filter(Province_State == state_name, Admin2 != "Unassigned", Admin2 != paste0("Out of ",state_name), Admin2 != paste0("Out of ",state_name_short))
  
  state_death_rate = us_death_rate %>%
    filter(Province_State == state_name, Admin2 != "Unassigned", Admin2 != paste0("Out of ",state_name), Admin2 != paste0("Out of ",state_name_short))
  
  n_cols = dim(state_death_rate)[2]
  
  if(duration > (n_cols-12)){
    return("Error: the argument Duration is too large")
  }
  
  start_date_index = which(all_dates==start_date)
  
  end_training_date_index = start_date_index + training_length - 1
  
  end_date_index = start_date_index + duration - 1
  
  select_criterion = which(state_death[, (12 + end_training_date_index)] >=criterion_death)
  
  # select_criterion = which(state_death[, dim(state_death)[2]] >=criterion_death)
  
  state_death_rate_selected = state_death_rate[select_criterion,]
  county_name_all = state_death_rate_selected$Admin2
  
  if (data_type == "death_rate"){
    
    output = matrix(0,length(select_criterion), duration)
    
    for (i in 1:length(select_criterion)){
      output[i,] = as.numeric(state_death_rate_selected[ i , (12 + start_date_index):(12 + end_date_index)])
    }
  } else if(data_type == "death"){
    
    state_death_selected = state_death[select_criterion,]
    county_name_all = state_death_selected$Admin2
    
    output = matrix(0,length(select_criterion), duration)
    
    for (i in 1:length(select_criterion)){
      output[i,] = as.numeric(state_death_selected[ i ,  (12 + start_date_index):(12 + end_date_index)])
    }
  }
  
  output[is.na(output)] = 0
  output[is.infinite(output)] = 0
  
  k = dim(output)[1]
  for(i in 1:k){
    for(j in 1:(dim(output)[2]-1)){
      if(output[i, j+1]<output[i, j]){
        output[i, j+1] = output[i, j]
      }
    }
  }
  
  if (smoothness){
    output_smooth_all = matrix(0,k,duration)
    for(i in 1:k){
      output_each = output[i,]
      output_window = zoo(output_each)
      output_smoothed = c(output_each[1:2], rollapply(output_window, width = 5, by = 1, FUN = mean, align = "left"),output_each[(length(output_each)-1):length(output_each)])
      output_smoothed[length(output_smoothed)-1] = mean(c(output_smoothed[length(output_smoothed)-2], output_smoothed[length(output_smoothed)]))
      output_smoothed[2] = mean(c(output_smoothed[1],output_smoothed[3]))
      output_smooth_all[i,] = output_smoothed
    }
    
    
    results = list(1:2)
    results[[1]] = output_smooth_all
    results[[2]] = county_name_all
    
    return(results)
  } else {
    results = list(1:2)
    results[[1]] = output
    results[[2]] = county_name_all
    
    return(results)
  }
}





SIRDC <- function(time, state, parms) {
  gamma = parms[[1]]
  theta = parms[[2]]
  delta = parms[[3]]
  N = parms[[4]]
  betafun = parms[[5]]
  
  par <- as.list(c(state))
  with(par, {
    beta <- betafun(time)
    dS <- -beta/N * I * S
    dI <- beta/N * I * S - gamma * I
    dR <- gamma * I - theta * R
    dD <- delta * theta * R
    dC <- (1-delta) * theta * R
    list(c(dS, dI, dR, dD, dC))
  })
}




clean_data = function(data){
  data_clean = data
  for(i in 1:(length(data)-1)){
    if (data_clean[i+1]<data_clean[i]){
      data_clean[i+1] = data_clean[i]
    }
  }
  return(data_clean)
}

##########################
# functions for estimating beta_t in approximation method
##########################



find_root_beta = function(beta_t_1, param,N, gamma){
  S_t_1 = param[1]
  S_t_2 = param[2]
  I_t_1 = param[3]
  
  I_t_2 = I_t_1 * exp(beta_t_1 * (S_t_1 + S_t_2) / (2*N) - gamma)
  
  results = S_t_2/S_t_1 - exp((beta_t_1/N* -(I_t_1+I_t_2)/2))
}


loss_approx_beta = function(param, death_cases, confirmed_cases,unadjusted_confirm_selected_smoothed, N_population, trianing_length, fixed_global_params, penalty = F,
                            fitted_days=length(death_cases),weight_loss=F){
  # ratio = exp(param[1])
  I_0 = param[1]
  R_0 = param[2]
  
  gamma = fixed_global_params[1]
  theta = fixed_global_params[2]
  delta = fixed_global_params[3]
  
  ratio = confirmed_cases[1]/(I_0+ R_0+ death_cases[1])
  
  estimated_confirm = confirmed_cases/ratio
  
  for(i in 1:length(estimated_confirm)){
    if (estimated_confirm[i] < unadjusted_confirm_selected_smoothed[i]){
      estimated_confirm[i] = unadjusted_confirm_selected_smoothed[i]
    }
  }
  
  S_t_seq = N_population - estimated_confirm
  
  # S_t_seq = N_population - confirmed_cases/ratio
  
  init_for_beta = c(S_t_seq[1], I_0, R_0, death_cases[1], 0)
  
  param_record_approx_for_beta = matrix(0, 5, trianing_length)
  param_record_approx_for_beta[,1] = init_for_beta
  param_record_approx_for_beta[1,] = S_t_seq
  
  approx_beta_seq = rep(0, trianing_length-1)
  
  for (i in 1:(trianing_length-1)){
    S_t_1 = param_record_approx_for_beta[1,i]
    S_t_2 = param_record_approx_for_beta[1,i+1]
    I_t_1 = param_record_approx_for_beta[2,i]
    R_t_1 = param_record_approx_for_beta[3,i]
    D_t_1 = param_record_approx_for_beta[4,i]
    C_t_1 = param_record_approx_for_beta[5,i]
    
    if(I_t_1<1){
      I_t_1 = 1
    }
    if(R_t_1<1){
      R_t_1 = 1
    }
    
    beta_t_1_2 = uniroot(find_root_beta, c(0, 10^6), tol = 0.0001, param = c(S_t_1, S_t_2, I_t_1), N = N_population, gamma=gamma)
    # I_t_2 = uniroot(find_root_I_t_2, c(0, N_population), tol = 0.0001, param = c(S_t_1, beta_t_1_2$root, I_t_1), N = N_population, gamma=gamma)
    
    I_t_2 = I_t_1 * exp(beta_t_1_2$root*(S_t_1 + S_t_2)/(2*N_population) - gamma)
    R_t_2 = (2-theta)/(2+theta)*R_t_1 + gamma/(2+theta)*(I_t_1+I_t_2)
    D_t_2 = D_t_1 + delta*theta*(R_t_1+R_t_2)/2
    C_t_2 = C_t_1 + (1-delta)*theta*(R_t_1+R_t_2)/2
    
    param_record_approx_for_beta[2:5, i+1] = c(I_t_2, R_t_2, D_t_2, C_t_2)
    approx_beta_seq[i] = beta_t_1_2$root
  }
  
  if(penalty){
    negative_index = which(param_record_approx_for_beta[4,] <= death_cases)
    positive_index=  which(param_record_approx_for_beta[4,] > death_cases)
    
    loss = sqrt( 10 * sum(((param_record_approx_for_beta[4,] - death_cases)[negative_index])^2) + 1 * sum(((param_record_approx_for_beta[4,] - death_cases)[positive_index])^2))
    return(loss)
    
  } else{
    #return(sqrt(sum((param_record_approx_for_beta[4,] - death_cases)^2)))
    index_selected_fitted=(length(death_cases)-fitted_days+1):length(death_cases)
    if(fitted_days>length(death_cases)){
      index_selected_fitted=1:length(death_cases)
    }
    
    if(!weight_loss){
      return(sqrt(sum((param_record_approx_for_beta[4, index_selected_fitted] - death_cases[index_selected_fitted])^2)))
    }else{
      weights=seq(length(index_selected_fitted),1,-1)
      return(sqrt(sum( ((param_record_approx_for_beta[4, index_selected_fitted] - death_cases[index_selected_fitted])^2)/weights^2)))
    }
  }
  
}
find_root_S_t_2 = function(S_t_2, param, N, gamma){
  S_t_1 = param[1]
  beta_t_1 = param[2]
  I_t_1 = param[3]
  
  I_t_2 = I_t_1 * exp(beta_t_1 * (S_t_1 + S_t_2) / (2*N) - gamma)
  
  results = S_t_2/S_t_1 - exp((beta_t_1/N* -(I_t_1+I_t_2)/2))
}

