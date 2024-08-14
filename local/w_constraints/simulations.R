#set working directory, load data, source code
setwd("/Users/masonyoungblood/Documents/Work/Summer_2021/Speed Climbing/SpeedClimbingABM")
load("data.RData")
grid <- read.csv("grid.csv")/1000
source("SpeedClimbingABM.R")

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
pop_data$gender <- data$gender[match(pop_data$ID, data$athlete)]
pop_data <- pop_data[, c(1, 5, 2, 3, 4)]

#get years
years <- sort(unique(data$year))

#get population sizes
n <- unlist(lapply(1:length(years), function(x){nrow(data[which(data$year == years[x]), ])}))

#number of simulations
n_sim <- 10000

#set random seed
set.seed(12345)

#set priors
priors <- data.frame(innov_prob = runif(n_sim, 0, 1),
                     learn_prob = runif(n_sim, 0, 1),
                     n_top = runif(n_sim, 0, 1),
                     improve_rate_m = rnorm(n_sim, 2, 1),
                     improve_rate_w = rnorm(n_sim, 2, 1),
                     constraint_a = truncnorm::rtruncnorm(n_sim, 0, Inf, 0, 1),
                     constraint_b = truncnorm::rtruncnorm(n_sim, 0, Inf, 0, 1),
                     innov_x_year = rnorm(n_sim, 0, 1),
                     innov_x_times = rnorm(n_sim, 0, 1),
                     innov_x_pop = rnorm(n_sim, 0, 1),
                     learn_x_year = rnorm(n_sim, 0, 1),
                     learn_x_times = rnorm(n_sim, 0, 1),
                     learn_x_pop = rnorm(n_sim, 0, 1))

#run simulations with everything varied
simulations <- do.call(rbind, parallel::mclapply(1:n_sim, function(x){
  unlist(SpeedClimbingABM(n = n, years = years, pop_data = pop_data, grid = grid, n_holds = 20,
                          beta_true_prob = 1,
                          innov_prob = priors$innov_prob[x],
                          learn_prob = priors$learn_prob[x],
                          n_top = priors$n_top[x],
                          improve_rate_m = priors$improve_rate_m[x],
                          improve_rate_w = priors$improve_rate_w[x],
                          innov_x_year = priors$innov_x_year[x],
                          innov_x_times = priors$innov_x_times[x],
                          innov_x_pop = priors$innov_x_pop[x],
                          learn_x_year = priors$learn_x_year[x],
                          learn_x_times = priors$learn_x_times[x],
                          learn_x_pop = priors$learn_x_pop[x],
                          sum_stats = FALSE, plot = FALSE)[[-1]])
}, mc.cores = 7))

#reset working directory to target
setwd("/Users/masonyoungblood/Documents/Work/Summer_2021/Speed Climbing/SpeedClimbingABM/local/w_constraints")

#save priors and simulations
priors_and_simulations <- list(priors, simulations)
save(priors_and_simulations, file = "priors_and_simulations.RData")
