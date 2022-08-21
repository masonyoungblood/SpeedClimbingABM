#set working directory, load data, source code
setwd(system("pwd", intern = T))
load("data.RData")
source("SpeedClimbingABM.R")

#euclidean distance function, where a and b are lists of vectors, with optional weighted (not default)
euclidean <- function(a, b, weighted = FALSE){
  divs <- c(length(a):1)
  if(weighted){
    return(sum(sapply(1:length(a), function(x){sqrt(sum((a[[x]]-b[[x]])^2)/divs[[x]])})))
  } else{
    return(sum(sapply(1:length(a), function(x){sqrt(sum((a[[x]]-b[[x]])^2))})))
  }
}

#random seed
set.seed(12345)

#subset data to only include men
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
SpeedClimbingABM_slurm <- function(innov_prob, learn_prob, n_top, adj_poss, improve_rate_m, improve_rate_sd, improve_min){
  temp <- SpeedClimbingABM(n = n, years = years, pop_data = pop_data, n_holds = 20,
                           beta_true_prob = 1, innov_prob = innov_prob, learn_prob = learn_prob, n_top = n_top, adj_poss = adj_poss, 
                           improve_rate_m = improve_rate_m, improve_rate_sd = improve_rate_sd, improve_min = improve_min,
                           sum_stats = FALSE, plot = FALSE)
  weighted_euclidean(obs_stats[-1], temp[-1])
}

#store required packages
pkgs <- unique(getParseData(parse("SpeedClimbingABM.R"))$text[getParseData(parse("SpeedClimbingABM.R"))$token == "SYMBOL_PACKAGE"])

#number of simulations
n_sim <- 100000000

#set priors
priors <- data.frame(innov_prob = truncnorm::rtruncnorm(n_sim, a = 0, b = 0.5, mean = 0, sd = 0.1),
                     learn_prob = truncnorm::rtruncnorm(n_sim, a = 0, b = 0.5, mean = 0, sd = 0.1),
                     n_top = runif(n_sim, min = 1, max = 34),
                     adj_poss = runif(n_sim, 1, 2),
                     improve_rate_m = runif(n_sim, min = 1, max = 4),
                     improve_rate_sd = truncnorm::rtruncnorm(n_sim, a = 0, mean = 0, sd = 2),
                     improve_min = runif(n_sim, min = 0, max = 0.5))

#save priors
save(priors, file = "prior_table.RData")

#run simulations
slurm <- rslurm::slurm_apply(SpeedClimbingABM_slurm, priors, jobname = "main",
                             nodes = 5, cpus_per_node = 48, pkgs = pkgs,
                             global_objects = objects(), slurm_options = list(mem = 0))
