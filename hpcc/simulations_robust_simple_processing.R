#set working directory, load data, source code
setwd(system("pwd", intern = T))
load("data.RData")
grid <- read.csv("grid.csv")/1000
source("SpeedClimbingABM.R")

#euclidean distance function
euclidean <- function(obs_stats, sum_stats){
  return(sqrt(sum((obs_stats-sum_stats)^2)))
}

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

#load parameters
params <- readRDS("_rslurm_simulations_robust/params.RDS")

#load results and combine
results_0 <- readRDS(paste0("_rslurm_", i, "/results_0.RDS"))
results_1 <- readRDS(paste0("_rslurm_", i, "/results_1.RDS"))
results_2 <- readRDS(paste0("_rslurm_", i, "/results_2.RDS"))
results_3 <- readRDS(paste0("_rslurm_", i, "/results_3.RDS"))
results <- c(unlist(results_0), unlist(results_1), unlist(results_2), unlist(results_3))
rm(list = c("results_0", "results_1", "results_2", "results_3"))

#get known parameters
set.seed(12345)
known <- params[sample(nrow(params), 1), ]

#get simulated data from known parameter values
obs_stats <- SpeedClimbingABM(n = n, years = years, pop_data = pop_data, grid = grid, n_holds = 20,
                              beta_true_prob = 1, innov_prob = known$innov_prob, learn_prob = known$learn_prob,
                              n_top = 0.1, max_dist = 1.645, constraint_b = known$constraint_b,
                              improve_rate_m = known$improve_rate_m, improve_min = 0.3427374, sum_stats = FALSE, plot = FALSE)

#calculate distances
dists <- sapply(1:length(results), function(x){euclidean(unlist(obs_stats), results[[x]])})

#combine and save output
output <- list(params, known, dists)
save(output, file = "output.RData")
