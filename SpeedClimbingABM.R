
# FUNCTIONS ---------------------------------------------------------------

#bounded exponential function for athletic improvement
bounded_exp <- function(x, rate, min){return((1-min)*(rate/rate^x)+min)}

#logit functions for interactions
logit <- function(p){return(log(p/(1-p)))}
inv_logit <- function(l){return(exp(l)/(1+exp(l)))}

# #parameter definition for manual debugging
# n_holds <- 20
# beta_true_prob <- 1
# innov_prob <- 0.2
# innov_x_times <- 0
# innov_x_pop <- 0
# learn_prob <- 0.2
# learn_x_times <- 0
# learn_x_pop <- 0
# n_top <- 10
# max_dist <- 2
# constraint_a <- 1
# constraint_b <- 1
# improve_rate_m <- 2
# improve_rate_sd <- 0.5
# improve_min <- 0.8
# sd_multiplier <- 0.5
# sum_stats <- TRUE
# plot <- TRUE
# bw <- 1
# ylim <- 0.3
# quant_by <- 0.2

#function for the model
SpeedClimbingABM <- function(n, years, pop_data, n_holds, beta_true_prob, 
                             learn_prob, n_top, learn_x_times = 0, learn_x_pop = 0, learn_x_year = 0,
                             innov_prob, max_dist, constraint_a = 0, constraint_b = 0, grid,
                             innov_x_times = 0, innov_x_pop = 0, innov_x_year = 0, 
                             improve_rate_m, improve_rate_sd, improve_min, sd_multiplier = 0.5, 
                             sum_stats = TRUE, plot = TRUE, raw = FALSE, bw = 1, ylim = 0.3, quant_by = 0.1){
  #handle output booleans
  if(raw){
    sum_stats <- FALSE
  }
  
  #round integer params
  n_top <- round(n_top)

  #colors for plotting after each timepoint
  colors <- rainbow((length(n)-1)*1.25) #times 1.2 so it doesn't loop back around
  
  #get distances between each hold from grids object
  dists <- sapply(1:n_holds, function(x){sqrt((grid$x[x+1] - grid$x[x])^2 + (grid$y[x+1] - grid$y[x])^2)})

  #get booleans for whether learning and innovation needs to be individually calculated, and scale and center population size
  learn_bool <- learn_x_times != 0 | learn_x_pop != 0 | learn_x_year != 0
  innov_bool <- innov_x_times != 0 | innov_x_pop != 0 | innov_x_year != 0
  n_scale <- scale(n)
  y_scale <- scale(1:length(n))
  
  #initialize starting beta for all agents
  beta <- sample(c(TRUE, FALSE), n_holds, prob = c(beta_true_prob, 1-beta_true_prob), replace = TRUE)
  
  #initialize data table of agents
  climbers <- data.table::data.table(ID = pop_data$ID[which(pop_data$start == years[1])],
                                     ref_times = pop_data$time[which(pop_data$start == years[1])]/sum(beta),
                                     beta = lapply(1:n[1], function(x){beta}),
                                     seq_ratios = lapply(1:n[1], function(x){sort(truncnorm::rtruncnorm(n_holds, a = 0, mean = 1, sd = sd_multiplier))[rank(dists, ties.method = "first")]}),
                                     ath_imp = lapply(1:n[1], function(x){bounded_exp(1:length(n), truncnorm::rtruncnorm(1, a = 1, mean = improve_rate_m, sd = improve_rate_sd), improve_min)}),
                                     current_record = pop_data$time[which(pop_data$start == years[1])])
  
  #create output list
  output <- list()
  
  #store the output from initialization
  if(sum_stats){output[[1]] <- quantile(climbers$current_record, probs = seq(0, 1, quant_by))}
  if(!sum_stats & !raw){output[[1]] <- sort(climbers$current_record)}
  
  #if raw output is stored, also create counter for actualized learning and innovation rate
  if(raw){
    act_learn <- 0
    act_innov <- 0
  }
  
  #initialize plot
  if(plot){
    par(mar = c(4, 4, 1, 1))
    plot(density(climbers$current_record, bw = bw), xlab = "Time (s)", main = "", xlim = c(0, mean(climbers$current_record)), ylim = c(0, ylim))
  }
  
  #iterate over time
  for(i in 2:length(n)){
    #get top n climbers for each climber to compare themselves with
    top_climbers <- order(climbers$current_record)[1:n_top]
    
    #if needed, calculate the unique learn_prob and innov_prob for each old and new climber in this timestep
    t_scale <- scale(climbers$current_record)
    if(learn_bool){learn_prob_ind <- sapply(1:nrow(climbers), function(x){inv_logit(logit(learn_prob) + learn_x_times*t_scale[x] + learn_x_pop*n_scale[i] + learn_x_year*y_scale[i])})}
    if(innov_bool){innov_prob_ind <- sapply(1:nrow(climbers), function(x){inv_logit(logit(innov_prob) + innov_x_times*t_scale[x] + innov_x_pop*n_scale[i] + innov_x_year*y_scale[i])})}
    
    #get who will learn
    if(learn_bool){to_learn <- sapply(1:nrow(climbers), function(x){sample(c(TRUE, FALSE), 1, prob = c(learn_prob_ind[x], 1-learn_prob_ind[x]))})}
    if(!learn_bool){to_learn <- sample(c(TRUE, FALSE), nrow(climbers), prob = c(learn_prob, 1-learn_prob), replace = TRUE)}
    
    #for each climber who is sampled
    for(k in c(1:nrow(climbers))[which(to_learn)]){
      #find out which top climbers have different beta than climber k
      diff_top_climbers <- top_climbers[which(sapply(top_climbers, function(x){!identical(climbers$beta[[k]], climbers$beta[[top_climbers[x]]])}))]
      
      #if there are top climbers with different beta
      if(length(diff_top_climbers) > 0){
        #get possible new times for climber k assuming the beta and seq_ratios of the top climbers with different beta
        poss_new_times <- sapply(diff_top_climbers, function(x){sum((climbers$ref_times[x]*climbers$seq_ratios[[x]]*climbers$ath_imp[[x]][i])[climbers$beta[[x]]])})
        
        #if the lowest possible new time is better than the current record, then replace the beta and seq_ratio of climber k with beta from best beta from top climbers
        if(min(poss_new_times) < climbers$current_record[k]){
          climbers$beta[[k]] <- climbers$beta[[diff_top_climbers[which.min(poss_new_times)]]]
          climbers$seq_ratios[[k]] <- climbers$seq_ratios[[diff_top_climbers[which.min(poss_new_times)]]]
          if(raw){act_learn <- act_learn + 1}
        }
        
        rm(list = c("poss_new_times"))
      }
      
      rm(list = c("diff_top_climbers"))
    }
    
    #get who will innovate
    if(innov_bool){to_flip <- which(sapply(1:nrow(climbers), function(x){sample(c(TRUE, FALSE), 1, prob = c(innov_prob_ind[x], 1-innov_prob_ind[x]))}))}
    if(!innov_bool){to_flip <- which(sample(c(TRUE, FALSE), nrow(climbers), prob = c(innov_prob, 1-innov_prob), replace = TRUE))}
    
    #if individuals are going to change their beta
    if(length(to_flip) > 0){
      #go through individuals
      for(j in to_flip){
        #store original beta
        beta_a <- climbers$beta[[j]]
        seq_ratios_a <- climbers$seq_ratios[[j]]
        
        #store holds used (skipping first hold since it cannot be flipped)
        ok_holds <- which(beta_a)[-1]
        
        #if any holds are available
        if(length(ok_holds) >= 1){
          #get euclidean distances between adjacent holds for each option (column 1), and ratios between original and new path distances (column 2)
          euc_dists <- t(sapply(1:length(ok_holds), function(x){
            #store adjacent distances (with added TRUEs for first hold and final buzzer)
            adj_dists <- c(1, ok_holds, n_holds+1)-ok_holds[x]
            
            #get used hold below and above
            lower_adj <- c(1, ok_holds, n_holds+1)[which.min(abs(adj_dists[which(adj_dists < 0)]))]
            upper_adj <- c(1, ok_holds, n_holds+1)[which(adj_dists > 0)[which.min(adj_dists[which(adj_dists > 0)])]]
            
            #return euclidean distances of original and new paths
            orig_dist <- sqrt((grid$x[ok_holds[x]]-grid$x[lower_adj])^2+(grid$y[ok_holds[x]]-grid$y[lower_adj])^2) + sqrt((grid$x[upper_adj]-grid$x[ok_holds[x]])^2+(grid$y[upper_adj]-grid$y[ok_holds[x]])^2)
            new_dist <- sqrt((grid$x[upper_adj]-grid$x[lower_adj])^2+(grid$y[upper_adj]-grid$y[lower_adj])^2)
            
            #return euclidean distance and ratio
            return(c(new_dist, orig_dist/new_dist))
          }))
          
          #only include holds less than maximum distance
          ok_holds <- ok_holds[which(euc_dists[, 1] <= max_dist)]
          
          #store length for boolean statements
          ok_length <- length(ok_holds)
          
          #if there is anything left
          if(ok_length > 0){
            if(raw){act_innov <- act_innov + 1}
            
            #if there's only one
            if(ok_length == 1){
              #go ahead and flip it
              beta_to_flip <- ok_holds
            }
            
            if(ok_length > 1){
              #subset distances for weighted sampling
              euc_dists <- euc_dists[which(euc_dists[, 1] <= max_dist), ]
              
              #choose position to flip
              beta_to_flip <- sample(ok_holds, 1, prob = ((1/rank(euc_dists[, 1]))^constraint_a)*((1/rank(euc_dists[, 2]))^constraint_b))
            }
            
            #get used hold below and above flipped beta for resampling times
            adj_dists <- which(c(beta_a, TRUE))-beta_to_flip
            lower_adj <- which(c(beta_a, TRUE))[which.min(abs(adj_dists[which(adj_dists < 0)]))]
            upper_adj <- which(c(beta_a, TRUE))[which(adj_dists > 0)[which.min(adj_dists[which(adj_dists > 0)])]]
            
            #flip a position into copied beta
            beta_a[beta_to_flip] <- FALSE
            
            #resample adjacent seq_ratios
            #seq_ratios_a[lower_adj] <- sort(truncnorm::rtruncnorm(n_holds, a = 0, mean = 1, sd = sd_multiplier))[rank(dists, ties.method = "first")][lower_adj]
            #if(upper_adj < n_holds+1){seq_ratios_a[upper_adj] <- sort(truncnorm::rtruncnorm(n_holds, a = 0, mean = 1, sd = sd_multiplier))[rank(dists, ties.method = "first")][upper_adj]}
            
            #overwrite current
            climbers$beta[[j]] <- beta_a
            climbers$seq_ratios[[j]] <- seq_ratios_a
            
            rm(list = c("euc_dists", "beta_to_flip", "adj_dists", "lower_adj", "upper_adj"))
          }
        }
        
        rm(list = c("beta_a", "seq_ratios_a", "ok_holds"))
      }
    }
    
    #get indices of climbers whose beta and sequence ratios will be randomly sampled for new climbers
    add_inds <- sample(1:nrow(climbers), length(which(pop_data$start == years[i])), replace = TRUE)
    
    #generate data table of new climbers
    new_climbers <- data.table::data.table(ID = pop_data$ID[which(pop_data$start == years[i])],
                                           ref_times = pop_data$time[which(pop_data$start == years[i])]/sapply(1:length(climbers$beta[add_inds]), function(x){sum(climbers$beta[add_inds][[x]])}),
                                           beta = climbers$beta[add_inds],
                                           seq_ratios = climbers$seq_ratios[add_inds],
                                           ath_imp = lapply(1:length(which(pop_data$start == years[i])), function(x){bounded_exp(1:length(n), truncnorm::rtruncnorm(1, a = 1, mean = improve_rate_m, sd = improve_rate_sd), improve_min)}),
                                           current_record = pop_data$time[which(pop_data$start == years[i])])
    
    #scale athletic improvement of new climbers to keep them on same trajectory as old climbers
    new_climbers$ath_imp <- lapply(1:nrow(new_climbers), function(x){new_climbers$ath_imp[[x]]/new_climbers$ath_imp[[x]][i]})
    
    #simulate leaving the sport
    climbers <- climbers[-which(climbers$ID %in% pop_data$ID[which(pop_data$end == years[i-1])]), ]
    
    #store current record of old climbers
    climbers$current_record <- sapply(1:nrow(climbers), function(x){sum((climbers$ref_times[x]*climbers$seq_ratios[[x]]*climbers$ath_imp[[x]][i])[climbers$beta[[x]]])})
    
    #combine data tables of old and new climbers
    climbers <- data.table::rbindlist(list(climbers, new_climbers))
    
    #plot
    if(plot){
      lines(density(climbers$current_record, bw = bw)$x, density(climbers$current_record, bw = bw)$y, type = "l", col = colors[i-1])
    }
    
    #store the output
    if(sum_stats){output[[i]] <- quantile(climbers$current_record, probs = seq(0, 1, quant_by))}
    if(!sum_stats){output[[i]] <- sort(climbers$current_record)}
    
    rm(list = c("top_climbers", "add_inds", "new_climbers", "to_learn", "to_flip"))
    if(learn_bool){rm(learn_prob_ind)}
    if(innov_bool){rm(innov_prob_ind)}
  }
  
  #return output
  if(sum_stats){return(do.call(rbind, output))}
  if(!sum_stats & !raw){return(output)}
  if(raw){return(list(climbers, act_learn = act_learn, act_innov = act_innov))}
}
