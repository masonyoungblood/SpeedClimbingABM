Speed Climbing ABM
================
Mason Youngblood

<img src="https://helios-i.mashable.com/imagery/articles/03XKxYMwkNiZJGiHUu2FWnw/hero-image.fill.size_1248x702.v1628089817.png" width="70%" style="display: block; margin: auto;" />

This agent-based model simulates a dynamic population of professional
speed climbers, and incorporates parameters for athletic improvement,
innovation of “beta” (or route sequence), and copying of other climbers’
beta.

## Athletic Improvement

The biggest challenge in developing this model was figuring out how to
simulation athletic improvement. I ended up assigning each climber a
beta distribution that changes over time and controls the amount of time
spent on each hold. In order to control the shape of the beta
distribution, I used the [following
function](https://stats.stackexchange.com/questions/12232/calculating-the-parameters-of-a-beta-distribution-using-the-mean-and-variance)
that converts a specified mean and variance into alpha and beta
parameters.

``` r
#function for calculating beta parameters
est_beta_params <- function(mean, var){
  #calculate alpha from mean and var
  alpha <- (((1-mean)/var)-(1/mean))*(mean^2)
  
  #calculate beta from alpha and mean
  beta <- alpha*((1/mean)-1)
  
  #return alpha and beta
  return(params = list(alpha = alpha, beta = beta))
}
```

With this function, we can easily construct beta distributions for each
climber. Here is a plotted exampe of what this looks like. This starting
distribution is heavily skewed towards 1, but as the distribution shifts
to the left the amount of time spent on each hold, and thus the overall
climbing time, will get lower and lower.

``` r
#generate range
range <- seq(0.01, 0.99, length.out = 1000)

#calculate parameters of beta dist
mean <- 0.99
var <- 0.001
beta_params <- est_beta_params(mean, var)

#construct and normalize dist
dist <- dbeta(range, beta_params$alpha, beta_params$beta)
dist <- (length(range)*dist)/sum(dist)

#plot it
par(mar = c(4, 4, 1, 1))
plot(range, dist, type = "l", xlab = "P", ylab = "Density")
```

![](README_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

In order to simulate improvement in this distribution, we will construct
an exponential distribution with a max bound drawn from the current
distribution, and randomly draw the new mean of the distribution from
it. Below is an example of an exponential distribution (in red) drawn
with a max bound drawn from a beta distribution (in blue).

![](README_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

If we iterate this out, we can see that the distribution shifts to the
left over time. In this case, the rate of the exponential curves
controls the rate of athletic improvement, where larger values mean that
the distributions shifts more slowly over time. For now, let’s set the
rate to 100. Note that, when calculating the parameters of the new beta
distribution, the variance is multiplied by the peak of the current
distribution to the power of 2, in order to “hit the brakes” on
directional improvement.

``` r
#set rate of exp curves when sampling
rate <- 100

#begin plot
par(mar = c(4, 4, 1, 1))
plot(range, dist, type = "l", xlab = "P", ylab = "Density")

#iterate
for(i in 1:20){
  #sample the max bound of exp distribution from the current distribution
  sampled_max <- sample(range, 1, prob = dist)
  
  #sample the exp distribution
  exp <- range[which(range <= sampled_max)]^rate
  
  #draw random point sample from exp distribution
  point_sample <- sample(range[which(range <= sampled_max)], 1, prob = exp)
  
  #create new dist with point sample as the mean
  beta_params <- est_beta_params(point_sample, var*(range[which.max(dist)]^2))
  dist <- dbeta(range, beta_params$alpha, beta_params$beta)
  
  #normalize it
  dist <- (length(range)*dist)/sum(dist)
  
  #plot new line
  lines(range, dist, type = "l", col = i+1)
  
  #remove temporary objects
  rm(list = c("sampled_max", "exp", "point_sample", "beta_params"))
}
```

![](README_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

Now let’s set the rate to 20 and see what happens.

![](README_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

As you can see a lower rate parameter for the exponential distribution
causes athletic improvement to happen at a faster rate over time, but
the rate of improvement slows down over time because how how we are
handling the variance. If we assign each climber a different rate
parameter for the exponential distributions then we can model variation
in climbers’ athletic ability.

## Beta & Reference Times

Now, how do we use these beta distributions in combination with hold
sequences? First, let’s look at a diagram of the standardized speed
wall.

<img src="https://walltopia.com/images/climbing_walls/speedwall-sketch.png" width="50%" style="display: block; margin: auto;" />

As you can see, the standardized speed wall has a total of 20 hand holds
and 11 foot holds. Speed climbers often “smear” their feet on the wall
or use hand holds for feet, so we will only be modeling the time spent
on hand holds. Each climber will be initialized with a vector of their
beta, or a TRUE/FALSE for whether they use each hand hold in the route,
along with a vector of reference times for the amount of time spent on
each hold. For now the beta vectors will start out as all TRUE, so that
all climbers start out using every hold in the route. The reference
times will be drawn from truncated normal distributions with a lower
bound at 0, a mean of the initial mean speed (let’s say 20 seconds)
divided by the number of holds on the route, and a standard deviation of
the mean divided by a parameter that controls the initial variation in
times across holds. Here is an example of how the beta and reference
time vectors are intialized.

``` r
#set number of holds
n_holds <- 20

#set initial mean speed
init_mean <- 20

#set probability of initial beta holds at 1 (all holds on the route)
beta_true_prob <- 1

#set parameter controlling the SD of ref times
sd_divider <- 10

#initialize starting beta
beta <- sample(c(TRUE, FALSE), n_holds, prob = c(beta_true_prob, 1-beta_true_prob), replace = TRUE)

#initialize ref times
ref_times <- truncnorm::rtruncnorm(n_holds, a = 0, mean = init_mean/n_holds, sd = (init_mean/n_holds)/sd_divider)

#print the vectors
beta
```

    ##  [1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ## [16] TRUE TRUE TRUE TRUE TRUE

``` r
ref_times
```

    ##  [1] 1.1220552 0.9902822 0.9048368 0.9206747 1.0805348 1.1248050 0.8082014
    ##  [8] 1.0841870 0.8132262 1.0091081 0.9177684 1.0258905 1.0037559 1.0447402
    ## [15] 0.9787296 0.8168301 0.9025306 1.1172947 1.0858090 0.9697394

In each timestep, each climber has the same `innov_prob` probability of
innovation. Innovations are changes to the beta of the route, in this
case flipping one of the booleans of the `beta` vector from TRUE to
FALSE (or vice versa). Not all boolean flips are possible. The parameter
`adj_poss` controls the number of successive FALSE booleans that are
allowed in the model (to simulate the fact that skipping three holds in
a row and still being able to complete the route is *extremely*
unlikely). When a boolean flip occurs (i.e. when a hold is added to or
dropped from the beta) the `ref_times` of the adjacent holds are
resampled, since the amount of time spent on the adjacent holds is
dependent on the presence of the added/dropped hold. When each climber
innovates, they compare the overall route time for their current beta
with the innovated beta, and if the innovated beta is better then they
overwrite their current beta. See `SpeedClimbingABM.R` for details.

## Social Learning

In each timestep, each climber also has the same `learn_prob`
probability of copying the beta of another climber. When a climber
copies another climber, they only have access to the top `n_top` fastest
climbers in the current timestep. They “try out” all of the betas of the
fastest climbers that are different from their own beta, and if one of
them yields a faster time then they overwrite their current `beta` and
`ref_times` with those of better beta. See `SpeedClimbingABM.R` for
details.

## Dynamic Population Size

We also allow the population size to grow over time. For now we just
have an initial population size of `n` climbers, a growth rate of
`growth_rate`, and a probability of leaving the sport of `leave_prob`.
Basically, at the end of each timestep a random subset of climbers leave
the sport, and new climbers enter the sport to bring the population up
to the size specified by `growth_rate`. The probability of a climber
leaving the sport is a function of how far their current record is from
the current best as well as how long they’ve been in the sport, so older
climbers with slower times are more likely to leave. Specifically it is
their `current_record` divided by the minimum `current_record` in the
population multiplied by `age`, or the number of timesteps that they
have been in the sport. See `SpeedClimbingABM.R` for details.

## Test Run

The ABM functions are in the `SpeedClimbingABM.R` file. Please refer to
the functions in this file for details.

``` r
source("SpeedClimbingABM.R")
```

Let’s do a test run of this model with a population size of 200 over 15
timesteps with a growth rate of 5%, with a 0.2 probability of climbers
leaving the sport during each timestep, an initial climbing time of 18
seconds across 20 holds where the climbers start out using every hold, a
0.2 probability of learning from the top 50 climbers, and a 0.1
probability of innovation where no more than 3 adjacent holds can be
skipped.

``` r
output <- SpeedClimbingABM(n = 200, t = 15, growth_rate = 0.05, leave_prob = 0.2, init_mean = 18,
                           n_holds = 20, beta_true_prob = 1, learn_prob = 0.2, n_top = 50,
                           innov_prob = 0.1, adj_poss = 3, rate_sd = 5, sd_divider = 10)
```

![](README_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

In the above plot, each distribution (from right to left) is the
distribution of `current_records` for climbers in the population in each
timestep.
