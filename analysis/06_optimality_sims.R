#set working directory, load data, source code
setwd("/Users/masonyoungblood/Documents/Work/Summer_2021/Speed Climbing/SpeedClimbingABM")
load("analysis/data_and_output/03_best_params/priors_and_simulations.RData")
load("climbing_times/best_climbing_times.RData"); data <- best_climbing_times; rm(best_climbing_times)
grid <- read.csv("grid.csv")/1000
source("SpeedClimbingABM.R")

#convert posteriors to usable table format
posterior_predictions <- jsonlite::read_json("analysis/data_and_output/04_best_params_inference/posterior_predictions.json")
posterior_predictions <- do.call(rbind, lapply(1:length(posterior_predictions), function(x){unlist(posterior_predictions[[x]])}))
posterior_predictions <- do.call(cbind, lapply(1:13, function(x){posterior_predictions[, x]*sd(priors_and_simulations[[1]][, x]) + mean(priors_and_simulations[[1]][, x])}))
posterior_predictions <- as.data.frame(posterior_predictions)
colnames(posterior_predictions) <- c("innov_prob", "learn_prob", "n_top", 
                                     "improve_min_m", "improve_min_w",
                                     "improve_rate_m", "improve_rate_w",
                                     "innov_x_year", "innov_x_times", "innov_x_pop",
                                     "learn_x_year", "learn_x_times", "learn_x_pop")

#prevent n_top from being < 0 or > 1
posterior_predictions$n_top[which(posterior_predictions$n_top < 0)] <- 0
posterior_predictions$n_top[which(posterior_predictions$n_top > 1)] <- 1

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
pop_data$time <- as.numeric(pop_data$time)
pop_data$gender <- toupper(data$gender[match(pop_data$ID, data$athlete)])
pop_data <- pop_data[, c(1, 5, 2, 3, 4)]

#get years
years <- sort(unique(data$year))

#get population sizes
n <- unlist(lapply(1:length(years), function(x){nrow(data[which(data$year == years[x]), ])}))

#number of simulations
#n_sim <- 1
n_sim <- 1000

#store dimension of matrix to explore
dim <- 11
#dim <- 20

#subset priors to include everything besides innov_prob and learn_prob
sub_priors <- posterior_predictions[1:n_sim, -c(1:2)]

#generate innovation and learning values, replacing 0 and 1 to avoid NAs in prob vector error
innov_probs <- seq(from = 0, to = 1, length.out = dim)
innov_probs[which(innov_probs == 0)] <- 0.00000001
innov_probs[which(innov_probs == 1)] <- 0.99999999
learn_probs <- seq(from = 0, to = 1, length.out = dim)
learn_probs[which(learn_probs == 0)] <- 0.00000001
learn_probs[which(learn_probs == 1)] <- 0.99999999

#loop through innovation and learning values and output master priors
priors <- do.call(rbind, lapply(1:dim, function(x){cbind(innov_prob = innov_probs[x], do.call(rbind, lapply(1:dim, function(x){cbind(learn_prob = learn_probs[x], sub_priors)})))}))

#run simulations with everything varied
optimality_sims <- do.call(rbind, parallel::mclapply(1:nrow(priors), function(x){
  temp <- SpeedClimbingABM(n = n, years = years, pop_data = pop_data, grid = grid, n_holds = 20,
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
                           sum_stats = FALSE, plot = FALSE, raw = TRUE)[[1]]
  temp_unique <- length(unique(sapply(1:nrow(temp), function(y){paste0(which(temp$beta[[y]]), collapse = " ")})))
  return(c(min(temp$current_record), median(temp$current_record), temp_unique))
}, mc.cores = 7))

#restructure and save simulations
colnames(optimality_sims) <- c("min", "median", "unique")
optimality_sims <- cbind(priors, optimality_sims)
save(optimality_sims, file = "analysis/data_and_output/06_optimality_sims/optimality_sims.RData")
