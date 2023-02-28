setwd(system("pwd", intern = T))
library(abcrf)
load("output/output.RData")

data <- cbind(output[[1]], output[[2]])
names(data)[ncol(data)] <- "dist"

sampsize <- 0.5

innov_tree <- regAbcrf(formula = innov_prob ~ dist, data = data, sampsize = nrow(data)*sampsize, mtry = 1, paral = TRUE, ncores = 20)
innov_pred <- predict(innov_tree, obs = data.frame(results = 0), training = data, sampsize = nrow(data)*sampsize, paral = TRUE, ncores = 20)
rm(innov_tree)

learn_tree <- regAbcrf(formula = learn_prob ~ dist, data = data, sampsize = nrow(data)*sampsize, mtry = 1, paral = TRUE, ncores = 20)
learn_pred <- predict(learn_tree, obs = data.frame(results = 0), training = data, sampsize = nrow(data)*sampsize, paral = TRUE, ncores = 20)
rm(learn_tree)

constraint_tree <- regAbcrf(formula = constraint_b ~ dist, data = data, sampsize = nrow(data)*sampsize, mtry = 1, paral = TRUE, ncores = 20)
constraint_pred <- predict(constraint_tree, obs = data.frame(results = 0), training = data, sampsize = nrow(data)*sampsize, paral = TRUE, ncores = 20)
rm(constraint_tree)

improve_tree <- regAbcrf(formula = improve_rate_m ~ dist, data = data, sampsize = nrow(data)*sampsize, mtry = 1, paral = TRUE, ncores = 20)
improve_pred <- predict(improve_tree, obs = data.frame(results = 0), training = data, sampsize = nrow(data)*sampsize, paral = TRUE, ncores = 20)
rm(improve_tree)

rf_preds <- list(innov_pred, learn_pred, constraint_pred, improve_pred)
save(rf_preds, file = "rf_preds.RData")
