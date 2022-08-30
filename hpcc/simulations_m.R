#set working directory, load data, source code
setwd(system("pwd", intern = T))
load("data.RData")
grid <- read.csv("grid.csv")/1000
source("SpeedClimbingABM.R")

#euclidean distance function where a and b are the unlisted vectors of times
euclidean <- function(a, b){
  return(sqrt(sum((a-b)^2)))
}

#random seed
set.seed(12345)

#subset data by gender
data <- data[which(data$gender == "M"), ]

#get all unique climbers
uniq_climbers <- unique(data$athlete)

#separate continuous sequences of years
seqs <- lapply(uniq_climbers, function(x){split(data$year[which(data$athlete == x)], cumsum(seq_along(data$year[which(data$athlete == x)]) %in% (which(diff(data$year[which(data$athlete == x)]) > 1) + 1)))})

#for each unique climber, iterate through their sequences, and and extract their ID, start year, end year, and time in start year (separate row per sequence)
pop_data <- data.table::data.table(do.call(rbind, lapply(1:length(uniq_climbers), function(i){
  t(sapply(1:length(seqs[[i]]), function(j){
    c(uniq_climbers[i], min(unlist(seqs[[i]][j])), max(unlist(seqs[[i]][j])), data$time[which(data$athlete == uniq_climbers[i] & data$year == min(unlist(seqs[[i]][j])))])
  }))
})))
colnames(pop_data) <- c("ID", "start", "end", "time")
pop_data

#get years
years <- sort(unique(data$year))

#get population sizes
n <- unlist(lapply(1:length(years), function(x){nrow(data[which(data$year == years[x]), ])}))

#get observed summary statistics
obs_stats <- lapply(years, function(x){sort(data$time[which(data$year == x)])})

#wrap SpeedClimbingABM in a simpler function for slurm, that outputs the sum of the euclidean distances between the distributions in each timepoint
SpeedClimbingABM_slurm <- function(innov_prob, innov_x_times, innov_x_pop, learn_prob, learn_x_times, learn_x_pop, n_top, adj_poss, improve_rate_m, improve_rate_sd, improve_min){
  temp <- SpeedClimbingABM(n = n, years = years, pop_data = pop_data, grid = grid, n_holds = 20,
                           beta_true_prob = 1, innov_prob = innov_prob, learn_prob = learn_prob, n_top = n_top, adj_poss = adj_poss, 
                           improve_rate_m = improve_rate_m, improve_rate_sd = improve_rate_sd, improve_min = improve_min,
                           sum_stats = FALSE, plot = FALSE)
  euclidean(unlist(obs_stats[-1]), unlist(temp[-1]))
}

#store required packages
pkgs <- unique(getParseData(parse("SpeedClimbingABM.R"))$text[getParseData(parse("SpeedClimbingABM.R"))$token == "SYMBOL_PACKAGE"])

#number of simulations per round
n_sim <- 10000

#tolerance level per round
tol <- 0.1

#number of rounds
rounds <- 500

for(i in 1:rounds){
  if(i == 1){
    #set priors
    priors <- data.frame(innov_prob = runif(n_sim, 0, 0.2),
                         innov_x_times = rnorm(n_sim, 0, 0.5), 
                         innov_x_pop = rnorm(n_sim, 0, 0.5), 
                         learn_prob = runif(n_sim, 0, 1),
                         learn_x_times = rnorm(n_sim, 0, 0.5),
                         learn_x_pop = rnorm(n_sim, 0, 0.5),
                         n_top = runif(n_sim, 1, 34),
                         constraint = runif(n_sim, 0, 6),
                         improve_rate_m = runif(n_sim, 1, 3),
                         improve_rate_sd = runif(n_sim, 0, 0.35),
                         improve_min = runif(n_sim, 0.15, 0.35))
  } else{
    #load parameters from previous round
    params <- readRDS(paste0("_rslurm_", i-1, "/params.RDS"))
    
    #get posteriors for each parameter for prior sampling
    innov_prob_post <- density(params$innov_prob[order(results)[1:(n_sim*tol)]], from = min(params$innov_prob), to = max(params$innov_prob))
    innov_x_times_post <- density(params$innov_x_times[order(results)[1:(n_sim*tol)]], from = min(params$innov_x_times), to = max(params$innov_x_times))
    innov_x_pop_post <- density(params$innov_x_pop[order(results)[1:(n_sim*tol)]], from = min(params$innov_x_pop), to = max(params$innov_x_pop))
    learn_prob_post <- density(params$learn_prob[order(results)[1:(n_sim*tol)]], from = min(params$learn_prob), to = max(params$learn_prob))
    learn_x_times_post <- density(params$learn_x_times[order(results)[1:(n_sim*tol)]], from = min(params$learn_x_times), to = max(params$learn_x_times))
    learn_x_pop_post <- density(params$learn_x_pop[order(results)[1:(n_sim*tol)]], from = min(params$learn_x_pop), to = max(params$learn_x_pop))
    n_top_post <- density(params$n_top[order(results)[1:(n_sim*tol)]], from = min(params$n_top), to = max(params$n_top))
    constraint_post <- density(params$constraint[order(results)[1:(n_sim*tol)]], from = min(params$constraint), to = max(params$constraint))
    improve_rate_m_post <- density(params$improve_rate_m[order(results)[1:(n_sim*tol)]], from = min(params$improve_rate_m), to = max(params$improve_rate_m))
    improve_rate_sd_post <- density(params$improve_rate_sd[order(results)[1:(n_sim*tol)]], from = min(params$improve_rate_sd), to = max(params$improve_rate_sd))
    improve_min_post <- density(params$improve_min[order(results)[1:(n_sim*tol)]], from = min(params$improve_min), to = max(params$improve_min))
    
    #remove objects
    rm(list = c("params", "results"))
    
    #set new priors by sampling from posteriors
    priors <- data.frame(innov_prob = sample(innov_prob_post$x, n_sim, replace = TRUE, prob = innov_prob_post$y),
                         innov_x_times = sample(innov_x_times_post$x, n_sim, replace = TRUE, prob = innov_x_times_post$y),
                         innov_x_pop = sample(innov_x_pop_post$x, n_sim, replace = TRUE, prob = innov_x_pop_post$y),
                         learn_prob = sample(learn_prob_post$x, n_sim, replace = TRUE, prob = learn_prob_post$y),
                         learn_x_times = sample(learn_x_times_post$x, n_sim, replace = TRUE, prob = learn_x_times_post$y),
                         learn_x_pop = sample(learn_x_pop_post$x, n_sim, replace = TRUE, prob = learn_x_pop_post$y),
                         n_top = sample(n_top_post$x, n_sim, replace = TRUE, prob = n_top_post$y),
                         constraint = sample(constraint_post$x, n_sim, replace = TRUE, prob = constraint_post$y),
                         improve_rate_m = sample(improve_rate_m_post$x, n_sim, replace = TRUE, prob = improve_rate_m_post$y),
                         improve_rate_sd = sample(improve_rate_sd_post$x, n_sim, replace = TRUE, prob = improve_rate_sd_post$y),
                         improve_min = sample(improve_min_post$x, n_sim, replace = TRUE, prob = improve_min_post$y))
    
    #remove objects
    rm(list = c("innov_prob_post", "innov_x_times_post", "innov_x_pop_post",
                "learn_prob_post", "learn_x_times_post", "learn_x_pop_post",
                "n_top_post", "constraint_post", "improve_rate_m_post", "improve_rate_sd_post", "improve_min_post"))
  }
  
  #run simulations
  slurm <- rslurm::slurm_apply(SpeedClimbingABM_slurm, priors, jobname = as.character(i),
                               nodes = 5, cpus_per_node = 48, pkgs = pkgs,
                               global_objects = objects(), slurm_options = list(mem = 0))
  
  #get simulation output
  results <- rslurm::get_slurm_out(slurm)
  results <- unlist(results)
  
  #remove objects
  rm(list = c("priors", "slurm"))
}
