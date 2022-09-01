#set working directory, load data, source code
setwd(system("pwd", intern = T))
load("data.RData")
grid <- read.csv("grid.csv")/1000
source("SpeedClimbingABM.R")

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

#load main simulation data
n_m <- 500
temp <- paste0("simulations_m/_rslurm_", n_m)
params <- readRDS(paste0(temp, "/params.RDS"))
results <- c(unlist(readRDS(paste0(temp, "/results_0.RDS"))),
             unlist(readRDS(paste0(temp, "/results_1.RDS"))), 
             unlist(readRDS(paste0(temp, "/results_2.RDS"))),
             unlist(readRDS(paste0(temp, "/results_3.RDS"))), 
             unlist(readRDS(paste0(temp, "/results_4.RDS"))))
rm(temp)

#set number of simulations kept below tolerance level
tol <- 1000

#get posteriors to sample priors from
n_top_post <- density(params$n_top[order(results)[1:tol]], from = min(params$n_top[order(results)[1:tol]]), to = max(params$n_top[order(results)[1:tol]]))
constraint_post <- density(params$constraint[order(results)[1:tol]], from = min(params$constraint[order(results)[1:tol]]), to = max(params$constraint[order(results)[1:tol]]))
improve_rate_m_post <- density(params$improve_rate_m[order(results)[1:tol]], from = min(params$improve_rate_m[order(results)[1:tol]]), to = max(params$improve_rate_m[order(results)[1:tol]]))
improve_rate_sd_post <- density(params$improve_rate_sd[order(results)[1:tol]], from = min(params$improve_rate_sd[order(results)[1:tol]]), to = max(params$improve_rate_sd[order(results)[1:tol]]))
improve_min_post <- density(params$improve_min[order(results)[1:tol]], from = min(params$improve_min[order(results)[1:tol]]), to = max(params$improve_min[order(results)[1:tol]]))

#store number of simulations per pixel
n_sim <- 100

#store dimension of matrix to explore
dim <- 50

#generate priors
sub_priors <- data.frame(n_top = sample(n_top_post$x, n_sim, replace = TRUE, prob = n_top_post$y),
                         constraint = sample(constraint_post$x, n_sim, replace = TRUE, prob = constraint_post$y),
                         improve_rate_m = sample(improve_rate_m_post$x, n_sim, replace = TRUE, prob = improve_rate_m_post$y),
                         improve_rate_sd = sample(improve_rate_sd_post$x, n_sim, replace = TRUE, prob = improve_rate_sd_post$y),
                         improve_min = sample(improve_min_post$x, n_sim, replace = TRUE, prob = improve_min_post$y))

#generate innovation and learning values
innov_probs <- seq(from = 0, to = 1, length.out = dim)
learn_probs <- seq(from = 0, to = 1, length.out = dim)

#loop through innovation and learning values and output master priors
priors <- do.call(rbind, lapply(1:dim, function(x){cbind(innov_prob = innov_probs[x], do.call(rbind, lapply(1:dim, function(x){cbind(learn_prob = learn_probs[x], sub_priors)})))}))

#wrap SpeedClimbingABM in a simpler function for slurm, that outputs the sum of the euclidean distances between the distributions in each timepoint
SpeedClimbingABM_slurm <- function(innov_prob, learn_prob, n_top, constraint, improve_rate_m, improve_rate_sd, improve_min){
  temp <- SpeedClimbingABM(n = n, years = years, pop_data = pop_data, grid = grid, n_holds = 20,
                           beta_true_prob = 1, innov_prob = innov_prob, learn_prob = learn_prob, n_top = n_top, constraint = constraint, 
                           improve_rate_m = improve_rate_m, improve_rate_sd = improve_rate_sd, improve_min = improve_min,
                           sum_stats = FALSE, plot = FALSE)
  median(temp[[length(n)]])
}

#store required packages
pkgs <- unique(getParseData(parse("SpeedClimbingABM.R"))$text[getParseData(parse("SpeedClimbingABM.R"))$token == "SYMBOL_PACKAGE"])

#run simulations
slurm <- rslurm::slurm_apply(SpeedClimbingABM_slurm, priors, jobname = "opt_m",
                             nodes = 5, cpus_per_node = 48, pkgs = pkgs,
                             global_objects = objects(), slurm_options = list(mem = 0))
