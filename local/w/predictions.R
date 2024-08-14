#set working directory, load data, source code
setwd("/Users/masonyoungblood/Documents/Work/Summer_2021/Speed Climbing/SpeedClimbingABM")
load("data.RData")
grid <- read.csv("grid.csv")/1000
source("SpeedClimbingABM.R")

#load packages
library(cito)
library(caret)

#logit and inverse logit functions with bounds
logit <- function(p, bounds){
  norm_p <- (p - bounds[1])/(bounds[2] - bounds[1])
  return(log(norm_p/(1-norm_p)))
}
inv_logit <- function(l, bounds){
  inv_l <- exp(l)/(1+exp(l))
  return((inv_l*(bounds[2] - bounds[1])) + bounds[1])
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

#load it
load("local/w/sum_stats.RData")

#https://topepo.github.io/caret/pre-processing.html
preprocessing <- preProcess(data.frame(sum_stats), method = c("center", "scale"))

#load models
load("local/w/models/innov_prob_model.RData")
innov_prob_model <- model
rm(model)
load("local/w/models/learn_prob_model.RData")
learn_prob_model <- model
rm(model)
load("local/w/models/n_top_model.RData")
n_top_model <- model
rm(model)
load("local/w/models/improve_rate_m_model.RData")
improve_rate_m_model <- model
rm(model)

#set number of iterations used for error, quantiles, etc.
n_iter <- 100

#pred intervals using mean +/- 1.96*RMSE: https://stats.stackexchange.com/questions/247879/using-mse-to-determine-prediction-intervals
#https://dtkaplan.github.io/SDS-book/mean-square-error.html#:~:text=The%20RMSE%20provides%20an%20operational,function%20output%20minus%20the%20RMSE.

# INNOVATION --------------------------------------------------------------

#get random values of parameter to simulate output from
known <- runif(n_iter, min = 0, max = 1)

#simulate output from known values
general_sims <- do.call(rbind, parallel::mclapply(1:n_iter, function(x){
  observed <- SpeedClimbingABM(n = n, years = years, pop_data = pop_data, grid = grid, n_holds = 20,
                               beta_true_prob = 1, max_dist = 1.645, improve_min = 0.3527184, 
                               innov_prob = known[x],
                               learn_prob = runif(1, min = 0, max = 1),
                               n_top = runif(1, min = 0, max = 1),
                               improve_rate_m = runif(1, min = 1, max = 3),
                               sum_stats = FALSE, plot = FALSE)
  return(predict(preprocessing, data.frame(t(matrix(unlist(observed[-1]))))))
}, mc.cores = 7))

#combine predicted and known values into a data frame
pred_and_known <- data.frame(pred = c(predict(innov_prob_model, general_sims, type = "response")), known = logit(known, c(0, 1)))

#get values of parameter to confirm model works for
explored_vals <- c(0.05, 0.2, 0.5, 0.8, 0.95)

#loop through these and get quantiles of the predicted values
known_quantiles <- cbind(value = explored_vals, do.call(rbind, lapply(1:length(explored_vals), function(x){
  #simulated output from known values
  explored_sims <- do.call(rbind, parallel::mclapply(1:n_iter, function(y){
    observed <- SpeedClimbingABM(n = n, years = years, pop_data = pop_data, grid = grid, n_holds = 20,
                                 beta_true_prob = 1, max_dist = 1.645, improve_min = 0.3527184, 
                                 innov_prob = explored_vals[x],
                                 learn_prob = runif(1, min = 0, max = 1),
                                 n_top = runif(1, min = 0, max = 1),
                                 improve_rate_m = runif(1, min = 1, max = 3),
                                 sum_stats = FALSE, plot = FALSE)
    return(predict(preprocessing, data.frame(t(matrix(unlist(observed[-1]))))))
  }, mc.cores = 7))
  
  #get quantiles from predictions of the model
  quantile(inv_logit(c(predict(innov_prob_model, explored_sims, type = "response")), c(0, 1)), probs = c(0.025, 0.1, 0.5, 0.9, 0.975))
})))

#return rmse, rsquared, and quantiles
innov_prob_gof <- list(rmse = sqrt(mean((pred_and_known$known - pred_and_known$pred)^2)),
                       rsquared = summary(lm(pred ~ known, pred_and_known))$r.squared,
                       quantiles = known_quantiles)

#remove temporary objects
rm(list = c("explored_vals", "known_quantiles", "pred_and_known", "general_sims", "known"))

# LEARNING ----------------------------------------------------------------

#get random values of parameter to simulate output from
known <- runif(n_iter, min = 0, max = 1)

#simulate output from known values
general_sims <- do.call(rbind, parallel::mclapply(1:n_iter, function(x){
  observed <- SpeedClimbingABM(n = n, years = years, pop_data = pop_data, grid = grid, n_holds = 20,
                               beta_true_prob = 1, max_dist = 1.645, improve_min = 0.3527184, 
                               innov_prob = runif(1, min = 0, max = 1),
                               learn_prob = known[x],
                               n_top = runif(1, min = 0, max = 1),
                               improve_rate_m = runif(1, min = 1, max = 3),
                               sum_stats = FALSE, plot = FALSE)
  return(predict(preprocessing, data.frame(t(matrix(unlist(observed[-1]))))))
}, mc.cores = 7))

#combine predicted and known values into a data frame
pred_and_known <- data.frame(pred = c(predict(learn_prob_model, general_sims, type = "response")), known = logit(known, c(0, 1)))

#get values of parameter to confirm model works for
explored_vals <- c(0.05, 0.2, 0.5, 0.8, 0.95)

#loop through these and get quantiles of the predicted values
known_quantiles <- cbind(value = explored_vals, do.call(rbind, lapply(1:length(explored_vals), function(x){
  #simulated output from known values
  explored_sims <- do.call(rbind, parallel::mclapply(1:n_iter, function(y){
    observed <- SpeedClimbingABM(n = n, years = years, pop_data = pop_data, grid = grid, n_holds = 20,
                                 beta_true_prob = 1, max_dist = 1.645, improve_min = 0.3527184, 
                                 innov_prob = runif(1, min = 0, max = 1),
                                 learn_prob = explored_vals[x],
                                 n_top = runif(1, min = 0, max = 1),
                                 improve_rate_m = runif(1, min = 1, max = 3),
                                 sum_stats = FALSE, plot = FALSE)
    return(predict(preprocessing, data.frame(t(matrix(unlist(observed[-1]))))))
  }, mc.cores = 7))
  
  #get quantiles from predictions of the model
  quantile(inv_logit(c(predict(learn_prob_model, explored_sims, type = "response")), c(0, 1)), probs = c(0.025, 0.1, 0.5, 0.9, 0.975))
})))

#return rmse, rsquared, and quantiles
learn_prob_gof <- list(rmse = sqrt(mean((pred_and_known$known - pred_and_known$pred)^2)),
                       rsquared = summary(lm(pred ~ known, pred_and_known))$r.squared,
                       quantiles = known_quantiles)

#remove temporary objects
rm(list = c("explored_vals", "known_quantiles", "pred_and_known", "general_sims", "known"))

# NUMBER TOP --------------------------------------------------------------

#get random values of parameter to simulate output from
known <- runif(n_iter, min = 0, max = 1)

#simulate output from known values
general_sims <- do.call(rbind, parallel::mclapply(1:n_iter, function(x){
  observed <- SpeedClimbingABM(n = n, years = years, pop_data = pop_data, grid = grid, n_holds = 20,
                               beta_true_prob = 1, max_dist = 1.645, improve_min = 0.3527184, 
                               innov_prob = runif(1, min = 0, max = 1),
                               learn_prob = runif(1, min = 0, max = 1),
                               n_top = known[x],
                               improve_rate_m = runif(1, min = 1, max = 3),
                               sum_stats = FALSE, plot = FALSE)
  return(predict(preprocessing, data.frame(t(matrix(unlist(observed[-1]))))))
}, mc.cores = 7))

#combine predicted and known values into a data frame
pred_and_known <- data.frame(pred = c(predict(n_top_model, general_sims, type = "response")), known = logit(known, c(0, 1)))

#get values of parameter to confirm model works for
explored_vals <- c(0.05, 0.2, 0.5, 0.8, 0.95)

#loop through these and get quantiles of the predicted values
known_quantiles <- cbind(value = explored_vals, do.call(rbind, lapply(1:length(explored_vals), function(x){
  #simulated output from known values
  explored_sims <- do.call(rbind, parallel::mclapply(1:n_iter, function(y){
    observed <- SpeedClimbingABM(n = n, years = years, pop_data = pop_data, grid = grid, n_holds = 20,
                                 beta_true_prob = 1, max_dist = 1.645, improve_min = 0.3527184, 
                                 innov_prob = runif(1, min = 0, max = 1),
                                 learn_prob = runif(1, min = 0, max = 1),
                                 n_top = explored_vals[x],
                                 improve_rate_m = runif(1, min = 1, max = 3),
                                 sum_stats = FALSE, plot = FALSE)
    return(predict(preprocessing, data.frame(t(matrix(unlist(observed[-1]))))))
  }, mc.cores = 7))
  
  #get quantiles from predictions of the model
  quantile(inv_logit(c(predict(n_top_model, explored_sims, type = "response")), c(0, 1)), probs = c(0.025, 0.1, 0.5, 0.9, 0.975))
})))

#return rmse, rsquared, and quantiles
n_top_gof <- list(rmse = sqrt(mean((pred_and_known$known - pred_and_known$pred)^2)),
                  rsquared = summary(lm(pred ~ known, pred_and_known))$r.squared,
                  quantiles = known_quantiles)

#remove temporary objects
rm(list = c("explored_vals", "known_quantiles", "pred_and_known", "general_sims", "known"))

# IMPROVEMENT RATE --------------------------------------------------------

#get random values of parameter to simulate output from
known <- runif(n_iter, min = 1, max = 3)

#simulate output from known values
general_sims <- do.call(rbind, parallel::mclapply(1:n_iter, function(x){
  observed <- SpeedClimbingABM(n = n, years = years, pop_data = pop_data, grid = grid, n_holds = 20,
                               beta_true_prob = 1, max_dist = 1.645, improve_min = 0.3527184, 
                               innov_prob = runif(1, min = 0, max = 1),
                               learn_prob = runif(1, min = 0, max = 1),
                               n_top = runif(1, min = 0, max = 1),
                               improve_rate_m = known[x],
                               sum_stats = FALSE, plot = FALSE)
  return(predict(preprocessing, data.frame(t(matrix(unlist(observed[-1]))))))
}, mc.cores = 7))

#combine predicted and known values into a data frame
pred_and_known <- data.frame(pred = c(predict(improve_rate_m_model, general_sims, type = "response")), known = logit(known, c(1, 3)))

#get values of parameter to confirm model works for
explored_vals <- c(1.1, 1.5, 2, 2.5, 2.9)

#loop through these and get quantiles of the predicted values
known_quantiles <- cbind(value = explored_vals, do.call(rbind, lapply(1:length(explored_vals), function(x){
  #simulated output from known values
  explored_sims <- do.call(rbind, parallel::mclapply(1:n_iter, function(y){
    observed <- SpeedClimbingABM(n = n, years = years, pop_data = pop_data, grid = grid, n_holds = 20,
                                 beta_true_prob = 1, max_dist = 1.645, improve_min = 0.3527184, 
                                 innov_prob = runif(1, min = 0, max = 1),
                                 learn_prob = runif(1, min = 0, max = 1),
                                 n_top = runif(1, min = 0, max = 1),
                                 improve_rate_m = explored_vals[x],
                                 sum_stats = FALSE, plot = FALSE)
    return(predict(preprocessing, data.frame(t(matrix(unlist(observed[-1]))))))
  }, mc.cores = 7))
  
  #get quantiles from predictions of the model
  quantile(inv_logit(c(predict(improve_rate_m_model, explored_sims, type = "response")), c(1, 3)), probs = c(0.025, 0.1, 0.5, 0.9, 0.975))
})))

#return rmse, rsquared, and quantiles
improve_rate_m_gof <- list(rmse = sqrt(mean((pred_and_known$known - pred_and_known$pred)^2)),
                           rsquared = summary(lm(pred ~ known, pred_and_known))$r.squared,
                           quantiles = known_quantiles)

#remove temporary objects
rm(list = c("explored_vals", "known_quantiles", "pred_and_known", "general_sims", "known"))

# GOODNESS OF FIT ---------------------------------------------------------

#combine goodness-of-fit data for all four parameters and save (rmse is on logit scale not original scale)
gof <- list(innov_prob_gof, learn_prob_gof, n_top_gof, improve_rate_m_gof)
save(gof, file = "local/w/gof.RData")

# REAL PREDICTION ---------------------------------------------------------

#get observed summary statistics
obs_stats <- lapply(years, function(x){sort(data$time[which(data$year == x)])})
obs_stats <- predict(preprocessing, (data.frame(t(matrix(unlist(obs_stats[-1]))))))

#set seed
set.seed(12345)

#get 1000 predictions for each value
innov_prob_pred <- sapply(1:1000, function(x){inv_logit(predict(innov_prob_model, obs_stats, type = "response"), c(0, 1))})
learn_prob_pred <- sapply(1:1000, function(x){inv_logit(predict(learn_prob_model, obs_stats, type = "response"), c(0, 1))})
n_top_pred <- sapply(1:1000, function(x){inv_logit(predict(n_top_model, obs_stats, type = "response"), c(0, 1))})
improve_rate_m_pred <- sapply(1:1000, function(x){inv_logit(predict(improve_rate_m_model, obs_stats, type = "response"), c(1, 3))})

#create and save data frame of predictions
predictions <- data.frame(innov_prob = innov_prob_pred, learn_prob = learn_prob_pred, n_top = n_top_pred, improve_rate_m = improve_rate_m_pred)
save(predictions, file = "local/w/predictions.RData")

