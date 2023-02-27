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
params <- readRDS("_rslurm_robust_simple/params.RDS")

#get simulated data from known parameter values
obs_stats <- unlist(SpeedClimbingABM(n = n, years = years, pop_data = pop_data, grid = grid, n_holds = 20,
                                     beta_true_prob = 1, innov_prob = 0.2, learn_prob = 0.3,
                                     n_top = 0.1, max_dist = 1.645, constraint_b = 1,
                                     improve_rate_m = 2, improve_min = 0.3427374, sum_stats = FALSE, plot = FALSE))

#iterate through individual files to preserve memory
for(i in 1:4){
  temp_results <- readRDS(paste0("_rslurm_robust_simple/results_", i-1, ".RDS"))
  temp_dists <- unlist(lapply(1:length(temp_results), function(x){c(dist(rbind(obs_stats, temp_results[[x]])))}))
  assign(paste0("dists_", i), temp_dists)
  save(temp_dists, file = paste0("output/dists_", i, ".RData"))
  rm(list = c("temp_results", "temp_dists"))
  gc()
}

#combine and save output
output <- list(params, dists = c(dists_1, dists_2, dists_3, dists_4))
save(output, file = "output/output.RData")
