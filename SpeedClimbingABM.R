
# FUNCTIONS ---------------------------------------------------------------

#function for calculating beta parameters
est_beta_params <- function(mean, var){
  #calculate alpha from mean and var
  alpha <- (((1-mean)/var)-(1/mean))*(mean^2)
  
  #calculate beta from alpha and mean
  beta <- alpha*((1/mean)-1)
  
  #return alpha and beta
  return(params = list(alpha = alpha, beta = beta))
}

#function for athletic improvement
athletic_improvement <- function(dist_range, dist_var, dist, rate){
  #sample the max bound of exponential distribution from the current distribution
  sampled_max <- sample(dist_range, 1, prob = dist)
  
  #sample the exponential distribution
  exp <- dist_range[which(dist_range <= sampled_max)]^rate
  
  #draw random point sample from exponential distribution
  point_sample <- sample(dist_range[which(dist_range <= sampled_max)], 1, prob = exp)
  
  #create new dist with point_sample as mean
    #the variance is multiplied by the peak of the current distribution to the power of 2
    #in order to "hit the brakes" on directional improvement
  beta_params <- est_beta_params(point_sample, dist_var*(dist_range[which.max(dist)]^2))
  new_dist <- dbeta(dist_range, beta_params$alpha, beta_params$beta)
  new_dist <- (length(range)*new_dist)/sum(new_dist) #normalize
  
  #return new distribution
  return(new_dist)
}

#function for the model
SpeedClimbingABM <- function(n, leave_prob, init_times, n_holds, beta_true_prob, learn_prob, n_top, innov_prob, adj_poss, rate_sd, sd_divider, sum_stats = TRUE, plot = TRUE, bw = 1, ylim = 0.3){
  #set t parameter
  t <- length(n)
  
  #colors for plotting after each timepoint
  colors <- rainbow(t*1.25) #times 1.2 so it doesn't loop back around
  
  #set range of rates of exp curves when sampling
  #which model individual athletic improvement rate (lower rate is better, larger SD allows for lower rates)
  rate_range <- c(0, 100)
  
  #initialize starting beta for all agents
  beta <- sample(c(TRUE, FALSE), n_holds, prob = c(beta_true_prob, 1-beta_true_prob), replace = TRUE)
  
  #initialize everything for the starting distributions
  dist_range <- seq(0.01, 0.99, length.out = 1000) #generate range
  dist_mean <- 0.99 #set mean of beta distribution
  dist_var <- 0.001 #set variance of beta distribution
  beta_params <- est_beta_params(dist_mean, dist_var) #get beta parameters
  dist <- dbeta(dist_range, beta_params$alpha, beta_params$beta) #set dist
  dist <- (length(dist_range)*dist)/sum(dist) #normalize dist
  
  #initialize ref times from observed times
  ref_times <- lapply(1:length(init_times), function(x){truncnorm::rtruncnorm(n_holds, a = 0, mean = (init_times[x])/n_holds, sd = ((init_times[x])/n_holds)/sd_divider)})
  
  #initialize data table of agents
  climbers <- data.table::data.table(age = rep(1, n[1]),
                                     beta = lapply(1:(n[1]), function(x){beta}),
                                     ref_times = ref_times,
                                     dists = lapply(1:(n[1]), function(x){dist}),
                                     rate = truncnorm::rtruncnorm(n[1], a = rate_range[1], b = rate_range[2], mean = rate_range[2], sd = rate_sd),
                                     current_record = rep(NA, n[1]))
  
  #store initial records before iterating over time
  climbers$current_record <- sapply(1:nrow(climbers), function(x){sum((climbers$ref_times[[x]]*sample(dist_range, n_holds, prob = climbers$dists[[x]], replace = TRUE))[climbers$beta[[x]]])})
  
  #create output list
  output <- list()
  
  #iterate over time
  for(i in 1:t){
    #if i is not the final timepoint
    if(i %in% 1:(t-1)){
      #store the output
      if(sum_stats){output[[i]] <- quantile(climbers$current_record)}
      if(!sum_stats){output[[i]] <- sort(climbers$current_record)}
      
      #get top n climbers for each climber to compare themselves with
      top_climbers <- order(climbers$current_record)[1:n_top]
      
      #for each climber who is sampled to learn and is not in the top n
      for(k in c(1:nrow(climbers))[-top_climbers][which(sample(c(TRUE, FALSE), nrow(climbers)-n_top, prob = c(learn_prob, 1-learn_prob), replace = TRUE))]){
        #find out which top climbers have different beta than climber k
        diff_top_climbers <- top_climbers[which(sapply(top_climbers, function(x){!identical(climbers$beta[[k]], climbers$beta[[top_climbers[x]]])}))]
        
        #get possible new times for climber k assuming the beta and ref_times of the top climbers with different beta
        poss_new_times <- sapply(diff_top_climbers, function(x){sum((climbers$ref_times[[x]]*sample(dist_range, n_holds, prob = climbers$dists[[k]], replace = TRUE))[climbers$beta[[x]]])})
        
        #if the lowest possible new time is better than the current record, then replace the beta and ref_time of climber k with beta from best beta from top climbers
        if(min(poss_new_times) < climbers$current_record[k]){
          climbers$beta[[k]] <- climbers$beta[[diff_top_climbers[which.min(poss_new_times)]]]
          climbers$ref_times[[k]] <- climbers$ref_times[[diff_top_climbers[which.min(poss_new_times)]]]
        }
        
        rm(list = c("diff_top_climbers", "poss_new_times"))
      }
      
      order(climbers$current_record)[1:n_top]
      
      #update distributions with athletic improvement
      climbers$dists <- lapply(1:nrow(climbers), function(x){athletic_improvement(dist_range, dist_var, climbers$dists[[x]], climbers$rate[[x]])})
      
      #who is going to innovate
      to_flip <- which(sample(c(TRUE, FALSE), nrow(climbers), prob = c(innov_prob, 1-innov_prob), replace = TRUE))
      
      #if individuals are going to change their beta
      if(length(to_flip) > 0){
        #go through individuals
        for(j in to_flip){
          #store original beta
          beta_a <- climbers$beta[[j]]
          ref_times_a <- climbers$ref_times[[j]]
          
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
          ref_times_b <- ref_times_a
          
          #choose position to flip
          beta_to_flip <- sample(ok_holds, 1)
          
          #flip a position into copied beta
          beta_b[beta_to_flip] <- !beta_b[beta_to_flip]
          
          #resample adjacent ref_times
          ref_times_b[ok_holds[which(ok_holds == beta_to_flip)-1]] <- truncnorm::rtruncnorm(1, a = 0, mean = mean(init_times)/n_holds, sd = (mean(init_times)/n_holds)/sd_divider)
          ref_times_b[ok_holds[which(ok_holds == beta_to_flip)+1]] <- truncnorm::rtruncnorm(1, a = 0, mean = mean(init_times)/n_holds, sd = (mean(init_times)/n_holds)/sd_divider)
          
          #decide which is better
          sum_a <- sum(ref_times_a[beta_a])
          sum_b <- sum(ref_times_b[beta_b])
          
          #if new one is better, then overwrite current
          if(sum_b < sum_a){
            climbers$beta[[j]] <- beta_b
            climbers$ref_times[[j]] <- ref_times_b
          }
          
          #remove temporary objects
          rm(list = c("beta_a", "ok_holds", "beta_b", "beta_to_flip", "sum_a", "sum_b"))
        }
      }
      
      #store current record in climbers table
      climbers$current_record <- sapply(1:nrow(climbers), function(x){sum((climbers$ref_times[[x]]*sample(dist_range, n_holds, prob = climbers$dists[[x]], replace = TRUE))[climbers$beta[[x]]])})
      
      #plot
      if(plot){
        if(i == 1){
          par(mar = c(4, 4, 1, 1))
          plot(density(climbers$current_record, bw = bw), xlab = "Time (s)", main = "", xlim = c(0, mean(init_times)), ylim = c(0, ylim))
        }
        if(i > 1){lines(density(climbers$current_record, bw = bw)$x, density(climbers$current_record, bw = bw)$y, type = "l", col = colors[i-1])}
      }
      
      #age up each climber
      climbers$age <- climbers$age+1
      
      #simulate climbers leaving the sport, where probability of leaving the sport is based function of age and slowness of current record
      climbers <- climbers[-sample(1:nrow(climbers), nrow(climbers)*(leave_prob[i]), prob = (climbers$current_record/min(climbers$current_record))*climbers$age), ]
      
      #get number of climbers to add
      add_n <- n[i+1]-nrow(climbers)
      
      #get indices of climbers whose beta and ref times will be randomly sampled for new climbers
      add_inds <- sample(1:nrow(climbers), add_n, replace = TRUE)
      
      #construct data table of new climbers, with beta and ref times from existing climbers and randomly chosen dists
      new_climbers <- data.table::data.table(age = rep(1, add_n),
                                             beta = climbers$beta[add_inds],
                                             ref_times = climbers$ref_times[add_inds],
                                             dists = climbers$dists[sample(1:nrow(climbers), add_n, replace = TRUE)],
                                             rate = truncnorm::rtruncnorm(add_n, a = rate_range[1], b = rate_range[2], mean = rate_range[2], sd = rate_sd),
                                             current_record = rep(NA, add_n))
      
      #add current records for new climbers
      new_climbers$current_record <- sapply(1:nrow(new_climbers), function(x){sum((new_climbers$ref_times[[x]]*sample(dist_range, n_holds, prob = new_climbers$dists[[x]], replace = TRUE))[new_climbers$beta[[x]]])})
      
      #combine data tables
      climbers <- data.table::rbindlist(list(climbers, new_climbers))
      
      #remove objects
      rm(list = c("top_climbers", "add_n", "add_inds", "new_climbers"))
    }
    
    #if i is the final timepoint
    if(i == t){
      #store the output
      if(sum_stats){output[[i]] <- quantile(climbers$current_record)}
      if(!sum_stats){output[[i]] <- sort(climbers$current_record)}
    }
  }
  
  #return output
  if(sum_stats){return(do.call(rbind, output))}
  if(!sum_stats){return(output)}
}
