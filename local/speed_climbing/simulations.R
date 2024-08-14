#set working directory, load data, source code
#setwd("/Users/masonyoungblood/Documents/Work/Summer_2021/Speed Climbing/SpeedClimbingABM")
setwd("~/speed_climbing")
load("data.RData")
grid <- read.csv("grid.csv")/1000
source("SpeedClimbingABM.R")

#load packages
library(caret)
library(gbm)
library(ranger)
library(glmnet)
library(cito)
#library(xgboost)
library(doParallel)

#create custom model specification for cito (torch) in caret
cito_fit <- function(x, y, wts, param, lev, last, weight, classProbs, ...){
  dat <- if(is.data.frame(x)) x else as.data.frame(x, stringsAsFactors = TRUE)
  dat$.outcome <- y
  if(param$regular){
    model <- cito::dnn(.outcome ~ ., 
                       data = dat,
                       lr = param$lr,
                       dropout = param$dropout,
                       hidden = c(1000, 1000),
                       loss = "mse", optimizer = "adam", activation = "relu",
                       alpha = 0.2, lambda = 0.001,
                       device = "cpu",
                       plot = FALSE,
                       ...)
  }
  if(!param$regular){
    model <- cito::dnn(.outcome ~ ., 
                       data = dat,
                       lr = param$lr,
                       dropout = param$dropout,
                       hidden = c(1000, 1000),
                       loss = "mse", optimizer = "adam", activation = "relu",
                       device = "cpu",
                       plot = FALSE,
                       ...)
  }
  class(model) <- "citodnn"
  return(model)
}
cito_predict <- function(modelFit, newdata, preProc = NULL, submodels = NULL){
  out <- as.numeric(cito:::predict.citodnn(modelFit, newdata, type = "response", device = "cpu"))
}
cito_dnn <- list(
  library = c("cito", "torch", "coro", "checkmate", "gridExtra", "parabar", "abind", "progress", "cli"),
  type = "Regression",
  parameters = data.frame(parameter = c("lr", "dropout", "regular"), 
                           class = c("numeric", "numeric", "boolean"), 
                           label = c("Learning Rate", "Dropout", "Regularization"))
)
cito_dnn$fit <- cito_fit
cito_dnn$predict <- cito_predict
cito_dnn$grid <- function(){}
cito_dnn$prob <- function(){}

#logit and inverse logit functions with bounds
logit <- function(p, bounds){
  norm_p <- (p - bounds[1])/(bounds[2] - bounds[1])
  return(log(norm_p/(1-norm_p)))
}
inv_logit <- function(l, bounds){
  inv_l <- exp(l)/(1+exp(l))
  return((inv_l*(bounds[2] - bounds[1])) + bounds[1])
}

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

# #number of simulations
# n_sim <- 100000
# 
# #set random seed
# set.seed(12345)
# 
# #set priors
# priors <- data.frame(innov_prob = runif(n_sim, 0, 1),
#                      learn_prob = runif(n_sim, 0, 1),
#                      n_top = runif(n_sim, 0, 1),
#                      improve_rate_m = runif(n_sim, 1, 3),
#                      improve_rate_w = runif(n_sim, 1, 3))
# 
# #run simulations with everything varied
# simulations <- parallel::mclapply(1:n_sim, function(x){
#   SpeedClimbingABM(n = n, years = years, pop_data = pop_data, grid = grid, n_holds = 20,
#                    beta_true_prob = 1,
#                    innov_prob = priors$innov_prob[x],
#                    learn_prob = priors$learn_prob[x],
#                    n_top = priors$n_top[x],
#                    improve_rate_m = priors$improve_rate_m[x],
#                    improve_rate_w = priors$improve_rate_w[x],
#                    sum_stats = FALSE, plot = FALSE)
# }, mc.cores = 95)
# 
# #save priors and simulations
# priors_and_simulations <- list(priors, simulations)
# save(priors_and_simulations, file = "priors_and_simulations.RData")

#load data
load("priors_and_simulations.RData")

#load into separate objects
priors <- priors_and_simulations[[1]]
simulations <- priors_and_simulations[[2]]

#set number of simulations
n_sim <- nrow(priors)

# #expand into summary statistics object for prediction
# sum_stats <- do.call(rbind, parallel::mclapply(1:n_sim, function(y){unlist(simulations[[y]][-1])}, mc.cores = 95))
# 
# #save it
# save(sum_stats, file = "sum_stats.RData")

#remove original simulations object and do garbage collection
rm(simulations)
gc()

#load it
load("sum_stats.RData")

#https://topepo.github.io/caret/pre-processing.html
preprocessing <- preProcess(data.frame(sum_stats), method = c("center", "scale"))
preprocessed_data <- predict(preprocessing, data.frame(sum_stats))

#remove original sum_stats object and do garbage collection
rm(sum_stats)
gc()

#set random seed
set.seed(12345)

#set parameters
n_tune <- 50
n_rounds <- 5
n_valid <- 1

# #register parallel backend for caret
# #https://stackoverflow.com/questions/65716884/how-to-prevent-rstudio-from-crashing
# #https://stackoverflow.com/questions/6543999/why-is-caret-train-taking-up-so-much-memory
# cl <- makeCluster(2)
# registerDoParallel(cl)
# 
# #create output object
# tuning <- list()
# 
# #run four forms of ABC for each parameter
# for(x in 1:5){
#   #construct data
#   if(x == 1){model_data <- cbind(param = logit(priors$innov_prob, c(0, 1)), preprocessed_data)}
#   if(x == 2){model_data <- cbind(param = logit(priors$learn_prob, c(0, 1)), preprocessed_data)}
#   if(x == 3){model_data <- cbind(param = logit(priors$n_top, c(0, 1)), preprocessed_data)}
#   if(x == 4){model_data <- cbind(param = logit(priors$improve_rate_m, c(1, 3)), preprocessed_data)}
#   if(x == 5){model_data <- cbind(param = logit(priors$improve_rate_w, c(1, 3)), preprocessed_data)}
#   
#   #hyperparameter tuning for gradient boosting
#   #this is the only not allowed to parallelize, because it increase computation time to do so
#   xgb_grid <- data.frame(nrounds = n_rounds,
#                          max_depth = sample(1:10, replace = TRUE, size = n_tune),
#                          eta = runif(n_tune, min = 0.001, max = 0.6),
#                          gamma = runif(n_tune, min = 0.0, max = 10.0),
#                          subsample = runif(n_tune, min = 0.25, max = 1.00),
#                          colsample_bytree = runif(n_tune, min = 0.30, max = 0.70),
#                          rate_drop = runif(n_tune, min = 0.01, max = 0.50),
#                          skip_drop = runif(n_tune, min = 0.05, max = 0.95),
#                          min_child_weight = sample(0:20, size = n_tune, replace = TRUE)) #modified from https://github.com/topepo/caret/blob/master/models/files/xgbDART.R
#   xgb_caret_model <- train(param ~ ., data = model_data, method = "xgbDART",
#                            trControl = trainControl(number = n_valid, verboseIter = FALSE, allowParallel = FALSE, trim = TRUE, returnData = FALSE),
#                            tuneGrid = xgb_grid, verbosity = 0)
# 
#   #print update
#   cat(paste0("Done with xgboost for parameter ", x, "...\n\n"))
#   
#   # #hyperparameter tuning for random forest
#   # rf_grid <- data.frame(min.node.size = sample(5:100, size = n_tune, replace = TRUE),
#   #                       mtry = sample(1:ncol(model_data[, -1]), size = n_tune, replace = TRUE),
#   #                       splitrule = sample(c("variance", "extratrees", "maxstat"), size = n_tune, replace = TRUE)) #modified from https://github.com/topepo/caret/blob/master/models/files/ranger.R
#   # rf_caret_model <- train(param ~ ., data = model_data, method = "ranger",
#   #                         num.trees = n_rounds,
#   #                         trControl = trainControl(number = n_valid, verboseIter = FALSE, allowParallel = TRUE, trim = TRUE, returnData = FALSE),
#   #                         tuneGrid = rf_grid)
#   # 
#   # #print update
#   # cat(paste0("Done with random forest for parameter ", x, "...\n\n"))
#   # 
#   # #hyperparameter tuning for glm with elastic net
#   # glm_caret_model <- train(param ~ ., data = model_data, method = "glmnet",
#   #                          trControl = trainControl(number = n_valid, search = "random", verboseIter = FALSE, allowParallel = TRUE, trim = TRUE, returnData = FALSE),
#   #                          tuneLength = n_tune)
#   # 
#   # #print update
#   # cat(paste0("Done with glmnet for parameter ", x, "...\n\n"))
#   
#   #hyperparameter tuning for deep learning
#   dnn_grid <- data.frame(lr = sample(c(0.01, 0.001, 0.0001, 0.00001), n_tune, replace = TRUE),
#                          dropout = runif(n_tune, min = 0, max = 0.9),
#                          regular = sample(c(TRUE, FALSE), n_tune, replace = TRUE))
#   dnn_caret_model <- train(y = model_data[, 1],
#                            x = model_data[, -1],
#                            method = cito_dnn,
#                            epochs = n_rounds,
#                            trControl = trainControl(number = n_valid, verboseIter = FALSE, allowParallel = TRUE, trim = TRUE, returnData = FALSE),
#                            tuneGrid = dnn_grid,
#                            verbose = FALSE)
#   
#   #print update
#   cat(paste0("Done with deep neural network for parameter ", x, "...\n\n"))
#   
#   #combine everything into an object
#   tuning[[x]] <- list(xgb = xgb_caret_model$results,
#                       # glm = glm_caret_model$results,
#                       # rf = rf_caret_model$results,
#                       dnn = dnn_caret_model$results)
#   
#   #remove all objects
#   rm(list = c("xgb_caret_model", "xgb_grid", 
#               # "rf_caret_model", "rf_grid", 
#               # "glm_caret_model", 
#               "dnn_caret_model", "dnn_grid"))
# 
#   #print update
#   cat(paste0("Done modeling parameter ", x, ": ", colnames(priors)[x], "\n\n"))
#   
#   #garbage collection
#   gc()
# }
# 
# #turn off parallelization
# stopCluster(cl)
# 
# #save output
# save(tuning, file = "tuning.RData")

#load data
load("tuning.RData")

#set index for best model (deep neural network)
y <- 2

#set random seed
set.seed(12345)

#set working directory to scratch directory
setwd("/scratch")

#set start time
time <- Sys.time()

#run four forms of approximate bayesian computation for each parameter
for(x in 1:5){
  #construct data
  if(x == 1){model_data <- cbind(param = logit(priors$innov_prob, c(0, 1)), preprocessed_data)}
  if(x == 2){model_data <- cbind(param = logit(priors$learn_prob, c(0, 1)), preprocessed_data)}
  if(x == 3){model_data <- cbind(param = logit(priors$n_top, c(0, 1)), preprocessed_data)}
  if(x == 4){model_data <- cbind(param = logit(priors$improve_rate_m, c(1, 3)), preprocessed_data)}
  if(x == 5){model_data <- cbind(param = logit(priors$improve_rate_w, c(1, 3)), preprocessed_data)}
  
  #get index for best parameter combination
  best_ind <- which.max(tuning[[x]][[y]]$Rsquared)
  
  #run best deep learning
  if(tuning[[x]][[y]]$regular[best_ind]){
    model <- dnn(param ~ ., data = model_data,
                 lr = tuning[[x]][[y]]$lr[best_ind],
                 dropout = tuning[[x]][[y]]$dropout[best_ind],
                 alpha = 0.2, lambda = 0.001,
                 epochs = 100,
                 hidden = c(1000, 1000),
                 loss = "mse", optimizer = "adam", activation = "relu",
                 device = "cpu", plot = FALSE)
  }
  if(!tuning[[x]][[y]]$regular[best_ind]){
    model <- dnn(param ~ ., data = model_data,
                 lr = tuning[[x]][[y]]$lr[best_ind],
                 dropout = tuning[[x]][[y]]$dropout[best_ind],
                 epochs = 100,
                 hidden = c(1000, 1000),
                 loss = "mse", optimizer = "adam", activation = "relu",
                 device = "cpu", plot = FALSE)
  }
  
  #save model
  save(model, file = paste0("models/", colnames(priors)[x], "_model.RData"))
  
  #print update
  cat(paste0("Done modeling parameter ", x, ": ", colnames(priors)[x], "\n\n"))
  
  #remove objects and do garbage collection
  rm(model)
  gc()
}

#print run time
Sys.time() - time
