#set working directory, load data, source code
setwd(system("pwd", intern = T))
load("data.RData")
grid <- read.csv("grid.csv")/1000
source("SpeedClimbingABM.R")

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

#number of simulations per round
n_sim <- 10000000

#set priors
priors <- data.frame(innov_prob = rbeta(n_sim, 1, 2),
                     learn_prob = rbeta(n_sim, 1, 2),
                     constraint_b = truncnorm::rtruncnorm(n_sim, a = 0, mean = 0, sd = 1),
                     improve_rate_m = runif(n_sim, 1, 3))

#wrap SpeedClimbingABM in a simpler function for slurm, that outputs the sum of the euclidean distances between the distributions in each timepoint
SpeedClimbingABM_slurm <- function(innov_prob, learn_prob, constraint_b, improve_rate_m){
  unlist(SpeedClimbingABM(n = n, years = years, pop_data = pop_data, grid = grid, n_holds = 20,
                          beta_true_prob = 1, innov_prob = innov_prob, learn_prob = learn_prob,
                          n_top = 0.1, max_dist = 1.645, constraint_b = constraint_b,
                          improve_rate_m = improve_rate_m, improve_min = 0.3427374, sum_stats = FALSE, plot = FALSE))
}

#store required packages
pkgs <- unique(getParseData(parse("SpeedClimbingABM.R"))$text[getParseData(parse("SpeedClimbingABM.R"))$token == "SYMBOL_PACKAGE"])

#run simulations
slurm <- rslurm::slurm_apply(SpeedClimbingABM_slurm, priors, jobname = "robust_simple",
                             nodes = 4, cpus_per_node = 45, pkgs = pkgs,
                             global_objects = objects(), slurm_options = list(mem = 0))
  