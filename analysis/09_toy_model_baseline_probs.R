
#this script is a toy model to demonstrate that our approach to modeling
#interactions leads to the intended behavior

#set working directory
setwd("~/Documents/Work/Summer_2021/Speed Climbing/SpeedClimbingABM/other/submissions/proc_b/revision_3")

#set random seed
set.seed(123)

#logit functions
logit <- function(p){return(log(p/(1-p)))}
inv_logit <- function(l){return(exp(l)/(1+exp(l)))}

#setup
n <- 1000 #pop size per simulation
k <- 1000 #number of simulations
x <- rnorm(n, 0, 1) #standardized covariate
mu <- 0.1 #parameter in the model
betas <- rnorm(k, 0, 1) #distribution of effects sampled from priors

#simulate
#for each of the k beta values, get the median value of mu across all n agents
median_mus <- sapply(1:k, function(y){
  mu_inds <- sapply(1:n, function(z){
    inv_logit(logit(mu) + betas[y]*x[z])
  })
  median(mu_inds)
})

#plot density of median mu values
plot(density(median_mus), main = "Distribution of median mu values", xlim = c(0.09, 0.11))

#save as png
png("toy_model.png", width = 5, height = 4, units = "in", res = 300)
plot(density(median_mus), main = "Distribution of median mu values", xlim = c(0.09, 0.11))
dev.off()
