#set working directory, load data, source code
setwd(system("pwd", intern = T))
load("data.RData")
grid <- read.csv("grid.csv")/1000
source("SpeedClimbingABM.R")

#euclidean distance function
euclidean <- function(obs_stats, sum_stats){
  return(sum(sapply(2:length(obs_stats), function(x){sqrt(sum((obs_stats[[x]]-sum_stats[[x]])^2))*(x/length(obs_stats))})))
}

#sample from arbitrary distribution function using grid search: https://jsta.github.io/rv/reference/rvdens.html
arb_sample <- function(n, object){
  #get grids and probabilities from bde object
  x <- bde::getdataPointsCache(object)
  y <- bde::getdensityCache(object)
  y[which(y < 0)] <- 0
  
  #sample from the grid and add scaled noise
  s <- sample(x, size = n, prob = y, replace = TRUE)
  s <- s + runif(length(s), (x[1]-x[2])/2, (x[2]-x[1])/2)
  
  #replace numbers out of bounds
  s[which(s < min(x))] <- min(x)
  s[which(s > max(x))] <- max(x)
  
  return(s)
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
SpeedClimbingABM_slurm <- function(innov_prob, innov_x_times, innov_x_pop, learn_prob, learn_x_times, learn_x_pop, n_top, constraint_a, constraint_b, improve_rate_m, improve_rate_sd){
  temp <- SpeedClimbingABM(n = n, years = years, pop_data = pop_data, grid = grid, n_holds = 20,
                           beta_true_prob = 1, innov_prob = innov_prob, learn_prob = learn_prob, 
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

#set resolution of density estimation
dens_res <- 1000

#set bounds of priors
bounds <- matrix(ncol = 2, byrow = TRUE,
                 c(0, 1, #innov_prob
                   -0.5, 0.5, #innov_x_times
                   -0.5, 0.5, #innov_x_pop
                   0, 1, #learn_prob
                   -0.5, 0.5, #learn_x_times
                   -0.5, 0.5, #learn_x_pop
                   1, 28, #n_top
                   -5, 5, #constraint_a
                   -5, 5, #constraint_b
                   1, 3, #improve_rate_m
                   0, 0.2)) #improve_rate_sd

for(i in 1:rounds){
  if(i == 1){
    #set priors
    priors <- data.frame(innov_prob = runif(n_sim, bounds[1, 1], bounds[1, 2]),
                         innov_x_times = runif(n_sim, bounds[2, 1], bounds[2, 2]),
                         innov_x_pop = runif(n_sim, bounds[3, 1], bounds[3, 2]),
                         learn_prob = runif(n_sim, bounds[4, 1], bounds[4, 2]),
                         learn_x_times = runif(n_sim, bounds[5, 1], bounds[5, 2]),
                         learn_x_pop = runif(n_sim, bounds[6, 1], bounds[6, 2]),
                         n_top = runif(n_sim, bounds[7, 1], bounds[7, 2]),
                         constraint_a = runif(n_sim, bounds[8, 1], bounds[8, 2]),
                         constraint_b = runif(n_sim, bounds[9, 1], bounds[9, 2]),
                         improve_rate_m = runif(n_sim, bounds[10, 1], bounds[10, 2]),
                         improve_rate_sd = runif(n_sim, bounds[11, 1], bounds[11, 2]))
  } else{
    #load parameters from previous round
    params <- readRDS(paste0("_rslurm_", i-1, "/params.RDS"))
    
    #get closes parameter values for density estimation
    innov_prob_post <- params$innov_prob[order(results)[1:(n_sim*tol)]]
    innov_x_times_post <- params$innov_x_times[order(results)[1:(n_sim*tol)]]
    innov_x_pop_post <- params$innov_x_pop[order(results)[1:(n_sim*tol)]]
    learn_prob_post <- params$learn_prob[order(results)[1:(n_sim*tol)]]
    learn_x_times_post <- params$learn_x_times[order(results)[1:(n_sim*tol)]]
    learn_x_pop_post <- params$learn_x_pop[order(results)[1:(n_sim*tol)]]
    n_top_post <- params$n_top[order(results)[1:(n_sim*tol)]]
    constraint_a_post <- params$constraint_a[order(results)[1:(n_sim*tol)]]
    constraint_b_post <- params$constraint_b[order(results)[1:(n_sim*tol)]]
    improve_rate_m_post <- params$improve_rate_m_post[order(results)[1:(n_sim*tol)]]
    improve_rate_sd_post <- params$improve_rate_sd_post[order(results)[1:(n_sim*tol)]]
    
    #get posteriors for each parameter for prior sampling
    innov_prob_post <- bde::bde(innov_prob_post, dataPointsCache = seq(min(innov_prob_post), max(innov_prob_post), length.out = dens_res), lower.limit = min(innov_prob_post), upper.limit = max(innov_prob_post), estimator = "betakernel", options = list(modified = TRUE, mbc = "jnl"))
    innov_x_times_post <- bde::bde(innov_x_times_post, dataPointsCache = seq(min(innov_x_times_post), max(innov_x_times_post), length.out = dens_res), lower.limit = min(innov_x_times_post), upper.limit = max(innov_x_times_post), estimator = "betakernel", options = list(modified = TRUE, mbc = "jnl"))
    innov_x_pop_post <- bde::bde(innov_x_pop_post, dataPointsCache = seq(min(innov_x_pop_post), max(innov_x_pop_post), length.out = dens_res), lower.limit = min(innov_x_pop_post), upper.limit = max(innov_x_pop_post), estimator = "betakernel", options = list(modified = TRUE, mbc = "jnl"))
    learn_prob_post <- bde::bde(learn_prob_post, dataPointsCache = seq(min(learn_prob_post), max(learn_prob_post), length.out = dens_res), lower.limit = min(learn_prob_post), upper.limit = max(learn_prob_post), estimator = "betakernel", options = list(modified = TRUE, mbc = "jnl"))
    learn_x_times_post <- bde::bde(learn_x_times_post, dataPointsCache = seq(min(learn_x_times_post), max(learn_x_times_post), length.out = dens_res), lower.limit = min(learn_x_times_post), upper.limit = max(learn_x_times_post), estimator = "betakernel", options = list(modified = TRUE, mbc = "jnl"))
    learn_x_pop_post <- bde::bde(learn_x_pop_post, dataPointsCache = seq(min(learn_x_pop_post), max(learn_x_pop_post), length.out = dens_res), lower.limit = min(learn_x_pop_post), upper.limit = max(learn_x_pop_post), estimator = "betakernel", options = list(modified = TRUE, mbc = "jnl"))
    n_top_post <- bde::bde(n_top_post, dataPointsCache = seq(min(n_top_post), max(n_top_post), length.out = dens_res), lower.limit = min(n_top_post), upper.limit = max(n_top_post), estimator = "betakernel", options = list(modified = TRUE, mbc = "jnl"))
    constraint_a_post <- bde::bde(constraint_a_post, dataPointsCache = seq(min(constraint_a_post), max(constraint_a_post), length.out = dens_res), lower.limit = min(constraint_a_post), upper.limit = max(constraint_a_post), estimator = "betakernel", options = list(modified = TRUE, mbc = "jnl"))
    constraint_b_post <- bde::bde(constraint_b_post, dataPointsCache = seq(min(constraint_b_post), max(constraint_b_post), length.out = dens_res), lower.limit = min(constraint_b_post), upper.limit = max(constraint_b_post), estimator = "betakernel", options = list(modified = TRUE, mbc = "jnl"))
    improve_rate_m_post <- bde::bde(improve_rate_m_post, dataPointsCache = seq(min(improve_rate_m_post), max(improve_rate_m_post), length.out = dens_res), lower.limit = min(improve_rate_m_post), upper.limit = max(improve_rate_m_post), estimator = "betakernel", options = list(modified = TRUE, mbc = "jnl"))
    improve_rate_sd_post <- bde::bde(improve_rate_sd_post, dataPointsCache = seq(min(improve_rate_sd_post), max(improve_rate_sd_post), length.out = dens_res), lower.limit = min(improve_rate_sd_post), upper.limit = max(improve_rate_sd_post), estimator = "betakernel", options = list(modified = TRUE, mbc = "jnl"))
    
    rm(list = c("params", "results"))
    
    #set new priors by sampling from posteriors
    priors <- data.frame(innov_prob = arb_sample(n_sim, innov_prob_post),
                         innov_x_times = arb_sample(n_sim, innov_x_times_post),
                         innov_x_pop = arb_sample(n_sim, innov_x_pop_post),
                         learn_prob = arb_sample(n_sim, learn_prob_post),
                         learn_x_times = arb_sample(n_sim, learn_x_times_post),
                         learn_x_pop = arb_sample(n_sim, learn_x_pop_post),
                         n_top = arb_sample(n_sim, n_top_post),
                         constraint_a = arb_sample(n_sim, constraint_a_post),
                         constraint_b = arb_sample(n_sim, constraint_b_post),
                         improve_rate_m = arb_sample(n_sim, improve_rate_m_post),
                         improve_rate_sd = arb_sample(n_sim, improve_rate_sd_post))
    
    rm(list = c("innov_prob_post", "innov_x_times_post", "innov_x_pop_post",
                "learn_prob_post", "learn_x_times_post", "learn_x_pop_post",
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
