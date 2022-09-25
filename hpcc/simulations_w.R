#set working directory, load data, source code
setwd(system("pwd", intern = T))
load("data.RData")
grid <- read.csv("grid.csv")/1000
source("SpeedClimbingABM.R")

#euclidean distance function
euclidean <- function(obs_stats, sum_stats){
  return(sum(sapply(2:length(obs_stats), function(x){sqrt(sum((obs_stats[[x]]-sum_stats[[x]])^2))*(x/length(obs_stats))})))
}

#subset data by gender
data <- data[which(data$gender == "W"), ]

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
SpeedClimbingABM_slurm <- function(innov_prob, innov_x_times, innov_x_pop, innov_x_year, learn_prob, learn_x_times, learn_x_pop, learn_x_year, n_top, constraint_a, constraint_b, improve_rate_m, improve_rate_sd){
  temp <- SpeedClimbingABM(n = n, years = years, pop_data = pop_data, grid = grid, n_holds = 20,
                           beta_true_prob = 1, innov_prob = innov_prob, learn_prob = learn_prob,
                           innov_x_times = innov_x_times, innov_x_pop = innov_x_pop, innov_x_year = innov_x_year,
                           learn_x_times = learn_x_times, learn_x_pop = learn_x_pop, learn_x_year = learn_x_year,
                           n_top = n_top, max_dist = 1.645, constraint_a = constraint_a, constraint_b = constraint_b,
                           improve_rate_m = improve_rate_m, improve_rate_sd = improve_rate_sd, 
                           improve_min = 0.4114153, sum_stats = FALSE, plot = FALSE)
  euclidean(obs_stats, temp)
}

#store required packages
pkgs <- unique(getParseData(parse("SpeedClimbingABM.R"))$text[getParseData(parse("SpeedClimbingABM.R"))$token == "SYMBOL_PACKAGE"])

#number of simulations per round
n_sim <- 2000

#tolerance level per round
tol <- 0.5

#number of rounds
rounds <- 500

for(i in 1:rounds){
  if(i == 1){
    #set priors
    priors <- data.frame(innov_prob = rbeta(n_sim, 1, 2),
                         innov_x_times = rnorm(n_sim, 0, 0.5),
                         innov_x_pop = rnorm(n_sim, 0, 0.5),
                         innov_x_year = rnorm(n_sim, 0, 0.5),
                         learn_prob = rbeta(n_sim, 1, 2),
                         learn_x_times = rnorm(n_sim, 0, 0.5),
                         learn_x_pop = rnorm(n_sim, 0, 0.5),
                         learn_x_year = rnorm(n_sim, 0, 0.5),
                         n_top = runif(n_sim, 1, 28),
                         constraint_a = rnorm(n_sim, 0, 1),
                         constraint_b = rnorm(n_sim, 0, 1),
                         improve_rate_m = runif(n_sim, 1, 3),
                         improve_rate_sd = rexp(n_sim, 10))
  } else{
    #load parameters from previous round
    params <- readRDS(paste0("_rslurm_", i-1, "/params.RDS"))
    
    #get posteriors for each parameter for prior sampling
    innov_prob_post <- density(params$innov_prob[order(results)[1:(n_sim*tol)]], from = min(params$innov_prob[order(results)[1:(n_sim*tol)]]), to = max(params$innov_prob[order(results)[1:(n_sim*tol)]]), n = 2^12, bw = "SJ")
    innov_x_times_post <- density(params$innov_x_times[order(results)[1:(n_sim*tol)]], from = min(params$innov_x_times[order(results)[1:(n_sim*tol)]]), to = max(params$innov_x_times[order(results)[1:(n_sim*tol)]]), n = 2^12, bw = "SJ")
    innov_x_pop_post <- density(params$innov_x_pop[order(results)[1:(n_sim*tol)]], from = min(params$innov_x_pop[order(results)[1:(n_sim*tol)]]), to = max(params$innov_x_pop[order(results)[1:(n_sim*tol)]]), n = 2^12, bw = "SJ")
    innov_x_year_post <- density(params$innov_x_year[order(results)[1:(n_sim*tol)]], from = min(params$innov_x_year[order(results)[1:(n_sim*tol)]]), to = max(params$innov_x_year[order(results)[1:(n_sim*tol)]]), n = 2^12, bw = "SJ")
    learn_prob_post <- density(params$learn_prob[order(results)[1:(n_sim*tol)]], from = min(params$learn_prob[order(results)[1:(n_sim*tol)]]), to = max(params$learn_prob[order(results)[1:(n_sim*tol)]]), n = 2^12, bw = "SJ")
    learn_x_times_post <- density(params$learn_x_times[order(results)[1:(n_sim*tol)]], from = min(params$learn_x_times[order(results)[1:(n_sim*tol)]]), to = max(params$learn_x_times[order(results)[1:(n_sim*tol)]]), n = 2^12, bw = "SJ")
    learn_x_pop_post <- density(params$learn_x_pop[order(results)[1:(n_sim*tol)]], from = min(params$learn_x_pop[order(results)[1:(n_sim*tol)]]), to = max(params$learn_x_pop[order(results)[1:(n_sim*tol)]]), n = 2^12, bw = "SJ")
    learn_x_year_post <- density(params$learn_x_year[order(results)[1:(n_sim*tol)]], from = min(params$learn_x_year[order(results)[1:(n_sim*tol)]]), to = max(params$learn_x_year[order(results)[1:(n_sim*tol)]]), n = 2^12, bw = "SJ")
    n_top_post <- density(params$n_top[order(results)[1:(n_sim*tol)]], from = min(params$n_top[order(results)[1:(n_sim*tol)]]), to = max(params$n_top[order(results)[1:(n_sim*tol)]]), n = 2^12, bw = "SJ")
    constraint_a_post <- density(params$constraint_a[order(results)[1:(n_sim*tol)]], from = min(params$constraint_a[order(results)[1:(n_sim*tol)]]), to = max(params$constraint_a[order(results)[1:(n_sim*tol)]]), n = 2^12, bw = "SJ")
    constraint_b_post <- density(params$constraint_b[order(results)[1:(n_sim*tol)]], from = min(params$constraint_b[order(results)[1:(n_sim*tol)]]), to = max(params$constraint_b[order(results)[1:(n_sim*tol)]]), n = 2^12, bw = "SJ")
    improve_rate_m_post <- density(params$improve_rate_m[order(results)[1:(n_sim*tol)]], from = min(params$improve_rate_m[order(results)[1:(n_sim*tol)]]), to = max(params$improve_rate_m[order(results)[1:(n_sim*tol)]]), n = 2^12, bw = "SJ")
    improve_rate_sd_post <- density(params$improve_rate_sd[order(results)[1:(n_sim*tol)]], from = min(params$improve_rate_sd[order(results)[1:(n_sim*tol)]]), to = max(params$improve_rate_sd[order(results)[1:(n_sim*tol)]]), n = 2^12, bw = "SJ")
    
    rm(list = c("params", "results"))
    
    #set new priors by sampling from posteriors
    priors <- data.frame(innov_prob = sample(innov_prob_post$x, n_sim, replace = TRUE, prob = innov_prob_post$y),
                         innov_x_times = sample(innov_x_times_post$x, n_sim, replace = TRUE, prob = innov_x_times_post$y),
                         innov_x_pop = sample(innov_x_pop_post$x, n_sim, replace = TRUE, prob = innov_x_pop_post$y),
                         innov_x_year = sample(innov_x_year_post$x, n_sim, replace = TRUE, prob = innov_x_year_post$y),
                         learn_prob = sample(learn_prob_post$x, n_sim, replace = TRUE, prob = learn_prob_post$y),
                         learn_x_times = sample(learn_x_times_post$x, n_sim, replace = TRUE, prob = learn_x_times_post$y),
                         learn_x_pop = sample(learn_x_pop_post$x, n_sim, replace = TRUE, prob = learn_x_pop_post$y),
                         learn_x_year = sample(learn_x_year_post$x, n_sim, replace = TRUE, prob = learn_x_year_post$y),
                         n_top = sample(n_top_post$x, n_sim, replace = TRUE, prob = n_top_post$y),
                         constraint_a = sample(constraint_a_post$x, n_sim, replace = TRUE, prob = constraint_a_post$y),
                         constraint_b = sample(constraint_b_post$x, n_sim, replace = TRUE, prob = constraint_b_post$y),
                         improve_rate_m = sample(improve_rate_m_post$x, n_sim, replace = TRUE, prob = improve_rate_m_post$y),
                         improve_rate_sd = sample(improve_rate_sd_post$x, n_sim, replace = TRUE, prob = improve_rate_sd_post$y))
    
    rm(list = c("innov_prob_post", "innov_x_times_post", "innov_x_pop_post", "innov_x_year_post",
                "learn_prob_post", "learn_x_times_post", "learn_x_pop_post", "learn_x_year_post",
                "n_top_post", "constraint_a", "constraint_b",
                "improve_rate_m_post", "improve_rate_sd_post"))
  }
  
  #run simulations
  slurm <- rslurm::slurm_apply(SpeedClimbingABM_slurm, priors, jobname = as.character(i),
                               nodes = 5, cpus_per_node = 48, pkgs = pkgs,
                               global_objects = objects(), slurm_options = list(mem = 0))
  
  #get simulation output
  results <- rslurm::get_slurm_out(slurm)
  results <- unlist(results)
  
  rm(list = c("priors", "slurm"))
}
