#set working directory, load data, source code
setwd(system("pwd", intern = T))
load("data.RData")
grid <- read.csv("grid.csv")/1000
source("SpeedClimbingABM.R")

#random seed
set.seed(12345)

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

#load main simulation data
n_w <- 500
temp <- paste0("simulations_w/_rslurm_", n_w)
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
innov_x_times_post <- density(params$innov_x_times[order(results)[1:tol]], from = min(params$innov_x_times[order(results)[1:tol]]), to = max(params$innov_x_times[order(results)[1:tol]]), n = 2^12, bw = "SJ")
innov_x_pop_post <- density(params$innov_x_pop[order(results)[1:tol]], from = min(params$innov_x_pop[order(results)[1:tol]]), to = max(params$innov_x_pop[order(results)[1:tol]]), n = 2^12, bw = "SJ")
innov_x_year_post <- density(params$innov_x_year[order(results)[1:tol]], from = min(params$innov_x_year[order(results)[1:tol]]), to = max(params$innov_x_year[order(results)[1:tol]]), n = 2^12, bw = "SJ")
learn_x_times_post <- density(params$learn_x_times[order(results)[1:tol]], from = min(params$learn_x_times[order(results)[1:tol]]), to = max(params$learn_x_times[order(results)[1:tol]]), n = 2^12, bw = "SJ")
learn_x_pop_post <- density(params$learn_x_pop[order(results)[1:tol]], from = min(params$learn_x_pop[order(results)[1:tol]]), to = max(params$learn_x_pop[order(results)[1:tol]]), n = 2^12, bw = "SJ")
learn_x_year_post <- density(params$learn_x_year[order(results)[1:tol]], from = min(params$learn_x_year[order(results)[1:tol]]), to = max(params$learn_x_year[order(results)[1:tol]]), n = 2^12, bw = "SJ")
n_top_post <- density(params$n_top[order(results)[1:tol]], from = min(params$n_top[order(results)[1:tol]]), to = max(params$n_top[order(results)[1:tol]]), n = 2^12, bw = "SJ")
constraint_a_post <- density(params$constraint_a[order(results)[1:tol]], from = min(params$constraint_a[order(results)[1:tol]]), to = max(params$constraint_a[order(results)[1:tol]]), n = 2^12, bw = "SJ")
constraint_b_post <- density(params$constraint_b[order(results)[1:tol]], from = min(params$constraint_b[order(results)[1:tol]]), to = max(params$constraint_b[order(results)[1:tol]]), n = 2^12, bw = "SJ")
improve_rate_m_post <- density(params$improve_rate_m[order(results)[1:tol]], from = min(params$improve_rate_m[order(results)[1:tol]]), to = max(params$improve_rate_m[order(results)[1:tol]]), n = 2^12, bw = "SJ")
improve_rate_sd_post <- density(params$improve_rate_sd[order(results)[1:tol]], from = min(params$improve_rate_sd[order(results)[1:tol]]), to = max(params$improve_rate_sd[order(results)[1:tol]]), n = 2^12, bw = "SJ")

#store number of simulations per pixel
n_sim <- 200

#store dimension of matrix to explore
dim <- 50

#generate priors
sub_priors <- data.frame(innov_x_times = sample(innov_x_times_post$x, n_sim, replace = TRUE, prob = innov_x_times_post$y),
                         innov_x_pop = sample(innov_x_pop_post$x, n_sim, replace = TRUE, prob = innov_x_pop_post$y),
                         innov_x_year = sample(innov_x_year_post$x, n_sim, replace = TRUE, prob = innov_x_year_post$y),
                         learn_x_times = sample(learn_x_times_post$x, n_sim, replace = TRUE, prob = learn_x_times_post$y),
                         learn_x_pop = sample(learn_x_pop_post$x, n_sim, replace = TRUE, prob = learn_x_pop_post$y),
                         learn_x_year = sample(learn_x_year_post$x, n_sim, replace = TRUE, prob = learn_x_year_post$y),
                         n_top = sample(n_top_post$x, n_sim, replace = TRUE, prob = n_top_post$y),
                         constraint_a = sample(constraint_a_post$x, n_sim, replace = TRUE, prob = constraint_a_post$y),
                         constraint_b = sample(constraint_b_post$x, n_sim, replace = TRUE, prob = constraint_b_post$y),
                         improve_rate_m = sample(improve_rate_m_post$x, n_sim, replace = TRUE, prob = improve_rate_m_post$y),
                         improve_rate_sd = sample(improve_rate_sd_post$x, n_sim, replace = TRUE, prob = improve_rate_sd_post$y))

#generate innovation and learning values
innov_probs <- seq(from = 0, to = 1, length.out = dim)
learn_probs <- seq(from = 0, to = 1, length.out = dim)

#loop through innovation and learning values and output master priors
priors <- do.call(rbind, lapply(1:dim, function(x){cbind(innov_prob = innov_probs[x], do.call(rbind, lapply(1:dim, function(x){cbind(learn_prob = learn_probs[x], sub_priors)})))}))

#wrap SpeedClimbingABM in a simpler function for slurm, that outputs the sum of the euclidean distances between the distributions in each timepoint
SpeedClimbingABM_slurm <- function(innov_prob, learn_prob,
                                   innov_x_times, innov_x_pop, innov_x_year,
                                   learn_x_times, learn_x_pop, learn_x_year,
                                   n_top, max_dist, constraint_a, constraint_b, improve_rate_m, improve_rate_sd){
  temp <- SpeedClimbingABM(n = n, years = years, pop_data = pop_data, grid = grid, n_holds = 20,
                           beta_true_prob = 1, innov_prob = innov_prob, learn_prob = learn_prob, 
                           innov_x_times_post = innov_x_times_post, innov_x_pop_post = innov_x_pop_post, innov_x_year_post = innov_x_year_post,
                           learn_x_times_post = learn_x_times_post, learn_x_pop_post = learn_x_pop_post, learn_x_year_post = learn_x_year_post,
                           n_top = n_top, constraint_a = constraint_a, constraint_b = constraint_b, max_dist = 1.645,
                           improve_rate_m = improve_rate_m, improve_rate_sd = improve_rate_sd, improve_min = 0.4114153,
                           sum_stats = FALSE, plot = FALSE)
  temp[[length(n)]]
}

#store required packages
pkgs <- unique(getParseData(parse("SpeedClimbingABM.R"))$text[getParseData(parse("SpeedClimbingABM.R"))$token == "SYMBOL_PACKAGE"])

#run simulations
slurm <- rslurm::slurm_apply(SpeedClimbingABM_slurm, priors, jobname = "opt_w",
                             nodes = 5, cpus_per_node = 48, pkgs = pkgs,
                             global_objects = objects(), slurm_options = list(mem = 0))
