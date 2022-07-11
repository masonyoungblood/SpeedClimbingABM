
# FUNCTIONS ---------------------------------------------------------------

#function for athletic improvement
athletic_improvement <- function(dist_stats, init_sd, brakes){
  #replace mean with point sample
  dist_stats[1] <- truncnorm::rtruncnorm(1, a = 0, b = 1, mean = dist_stats[1], sd = dist_stats[2])
  
  #replace sd with initial sd multiplied by the new mean to the power of brakes
  dist_stats[2] <- init_sd*(dist_stats[1]^brakes)
  
  #return new distribution
  return(dist_stats)
}

#function for the model
SpeedClimbingABM <- function(n, leave_prob, init_times, n_holds, beta_true_prob, learn_prob, n_top, learn_x_times = 0, learn_x_pop = 0, innov_prob, adj_poss, innov_x_times = 0, innov_x_pop = 0, init_sd, brakes, sd_multiplier = 0.5, sum_stats = TRUE, plot = TRUE, bw = 1, ylim = 0.3){
  #colors for plotting after each timepoint
  colors <- rainbow((length(n)-1)*1.25) #times 1.2 so it doesn't loop back around
  
  #get booleans for whether learning and innovation needs to be individually calculated, and normalize population size
  learn_bool <- learn_x_times > 0 | learn_x_pop > 0
  innov_bool <- innov_x_times > 0 | innov_x_pop > 0
  n_norm <- BBmisc::normalize(n, "range", c(-0.5, 0.5))
  
  #initialize starting beta for all agents
  beta <- sample(c(TRUE, FALSE), n_holds, prob = c(beta_true_prob, 1-beta_true_prob), replace = TRUE)
  
  #initialize data table of agents
  climbers <- data.table::data.table(age = rep(1, n[1]),
                                     beta = lapply(1:n[1], function(x){beta}),
                                     init_times = init_times,
                                     seq_ratios = lapply(1:n[1], function(x){truncnorm::rtruncnorm(n_holds, a = 0, mean = 1, sd = sd_multiplier)}),
                                     dists = lapply(1:n[1], function(x){c(1, init_sd)}),
                                     current_record = init_times)
  
  #create output list
  output <- list()
  
  #store the output from initialization
  if(sum_stats){output[[1]] <- quantile(climbers$current_record, probs = seq(0, 1, 0.1))}
  if(!sum_stats){output[[1]] <- sort(climbers$current_record)}
  
  #initialize plot
  if(plot){
    par(mar = c(4, 4, 1, 1))
    plot(density(climbers$current_record, bw = bw), xlab = "Time (s)", main = "", xlim = c(0, mean(init_times)), ylim = c(0, ylim))
  }
  
  #iterate over time
  for(i in 2:length(n)){
    #update distributions with athletic improvement
    climbers$dists <- lapply(1:nrow(climbers), function(x){athletic_improvement(climbers$dists[[x]], init_sd, brakes)})
    
    #get top n climbers for each climber to compare themselves with
    top_climbers <- order(climbers$current_record)[1:n_top]
    
    #if needed, calculate the unique learn_prob and innov_prob for each climber in this timestep
    t_norm <- BBmisc::normalize(1/climbers$current_record, "range", c(-0.5, 0.5)) #take inverse of current record so smaller is better
    if(learn_bool){learn_prob_ind <- sapply(1:nrow(climbers), function(x){learn_prob + learn_prob*t_norm[x]*learn_x_times + learn_prob*n_norm[i]*learn_x_pop})}
    if(innov_bool){innov_prob_ind <- sapply(1:nrow(climbers), function(x){innov_prob + innov_prob*t_norm[x]*innov_x_times + innov_prob*n_norm[i]*innov_x_pop})}
    
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
        poss_new_times <- sapply(diff_top_climbers, function(x){sum(((climbers$init_times[x]/n_holds)*climbers$seq_ratios[[x]]*truncnorm::rtruncnorm(n_holds, a = 0, b = 1, mean = climbers$dists[[x]][1], sd = climbers$dists[[x]][2]))[climbers$beta[[x]]])})
        
        #if the lowest possible new time is better than the current record, then replace the beta and seq_ratio of climber k with beta from best beta from top climbers
        if(min(poss_new_times) < climbers$current_record[k]){
          climbers$beta[[k]] <- climbers$beta[[diff_top_climbers[which.min(poss_new_times)]]]
          climbers$seq_ratios[[k]] <- climbers$seq_ratios[[diff_top_climbers[which.min(poss_new_times)]]]
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
        
        #if there are any FALSE positions
        if(length(which(beta_a == FALSE) > 0)){
          #iterate through beta and figure out which positions are okay to flip (do not create strings of FALSE that exceed adj_poss)
          ok_holds <- which(sapply(1:length(beta_a), function(x){
            #make copy of beta_a to test whether position is okay
            temp <- beta_a
            temp[x] <- !temp[x]
            
            #if the flip doesn't change everything to TRUE (rare case that throws warning)
            if(length(which(temp == FALSE) > 0)){
              #if the maximum string of FALSEs, assuming this position flips, is less than or equal to adj_poss, then it's okay to flip
              if(max(rle(temp)$lengths[which(rle(temp)$values == FALSE)]) <= adj_poss){
                return(TRUE)
              } else{
                return(FALSE)
              }
            } else{
              return(TRUE)
            }
            
            #remove temp object
            rm(temp)
          }))
        } else{
          ok_holds <- c(1:n_holds)
        }
        
        #make a copy
        beta_b <- beta_a
        seq_ratios_b <- seq_ratios_a
        
        #choose position to flip
        beta_to_flip <- sample(ok_holds, 1)
        
        #flip a position into copied beta
        beta_b[beta_to_flip] <- !beta_b[beta_to_flip]
        
        #resample adjacent seq_ratios
        seq_ratios_b[ok_holds[which(ok_holds == beta_to_flip)-1]] <- truncnorm::rtruncnorm(1, a = 0, mean = 1, sd = sd_multiplier)
        seq_ratios_b[ok_holds[which(ok_holds == beta_to_flip)+1]] <- truncnorm::rtruncnorm(1, a = 0, mean = 1, sd = sd_multiplier)
        
        #decide which is better
        sum_a <- sum(seq_ratios_a[beta_a])
        sum_b <- sum(seq_ratios_b[beta_b])
        
        #if new one is better, then overwrite current
        if(sum_b < sum_a){
          climbers$beta[[j]] <- beta_b
          climbers$seq_ratios[[j]] <- seq_ratios_b
        }
        
        #remove temporary objects
        rm(list = c("beta_a", "seq_ratios_a", "ok_holds", "beta_b", "seq_ratios_b", "beta_to_flip", "sum_a", "sum_b"))
      }
    }
    
    #store current record in climbers table
    climbers$current_record <- sapply(1:nrow(climbers), function(x){sum(((climbers$init_times[x]/n_holds)*climbers$seq_ratios[[x]]*truncnorm::rtruncnorm(n_holds, a = 0, b = 1, mean = climbers$dists[[x]][1], sd = climbers$dists[[x]][2]))[climbers$beta[[x]]])})
    
    #plot
    if(plot){
      lines(density(climbers$current_record, bw = bw)$x, density(climbers$current_record, bw = bw)$y, type = "l", col = colors[i-1])
    }
    
    #age up each climber
    climbers$age <- climbers$age+1
    
    #simulate climbers leaving the sport, where probability of leaving the sport is based function of age and slowness of current record
    climbers <- climbers[-sample(nrow(climbers), nrow(climbers)*leave_prob[i-1], prob = (climbers$current_record/min(climbers$current_record))*climbers$age), ]
    
    #get number of climbers to add
    add_n <- n[i]-nrow(climbers)

    #get indices of climbers whose beta and ref times will be randomly sampled for new climbers
    add_inds <- sample(1:nrow(climbers), add_n, replace = TRUE)

    #construct data table of new climbers, with beta and ref times from existing climbers and randomly chosen dists
    new_climbers <- data.table::data.table(age = rep(1, add_n),
                                           beta = climbers$beta[add_inds],
                                           init_times = climbers$init_times[add_inds],
                                           seq_ratios = climbers$seq_ratios[add_inds],
                                           dists = climbers$dists[sample(1:nrow(climbers), add_n, replace = TRUE)],
                                           current_record = rep(NA, add_n))
    
    #add current records for new climbers
    new_climbers$current_record <- sapply(1:nrow(new_climbers), function(x){sum(((new_climbers$init_times[x]/n_holds)*new_climbers$seq_ratios[[x]]*truncnorm::rtruncnorm(n_holds, a = 0, b = 1, mean = new_climbers$dists[[x]][1], sd = new_climbers$dists[[x]][2]))[new_climbers$beta[[x]]])})
    
    #combine data tables
    climbers <- data.table::rbindlist(list(climbers, new_climbers))
    
    #store the output
    if(sum_stats){output[[i]] <- quantile(climbers$current_record, probs = seq(0, 1, 0.1))}
    if(!sum_stats){output[[i]] <- sort(climbers$current_record)}
    
    #remove objects
    rm(list = c("top_climbers", "add_n", "add_inds", "new_climbers", "to_learn", "to_flip"))
    if(learn_bool){rm(learn_prob_ind)}
    if(innov_bool){rm(innov_prob_ind)}
  }
  
  #return output
  if(sum_stats){return(do.call(rbind, output))}
  if(!sum_stats){return(output)}
}
