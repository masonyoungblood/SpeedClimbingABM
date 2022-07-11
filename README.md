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
simulate athletic improvement. I ended up assigning each climber a
truncated normal distribution from 0 to 1 that changes over time and
controls the amount of time spent on each hold. Here is a plotted
example of what this looks like. The example starting distribution has a
mean of 1 and a standard deviation of 0.1, but as the distribution
shifts to the left the amount of time spent on each hold, and thus the
overall climbing time, will get lower and lower.

``` r
#store distribution statistics (mean and sd)
dist_stats <- c(1, 0.1)

#generate range
range <- seq(0, 1, by = 0.001)

#construct and normalize dist so max is 1
dist <- truncnorm::dtruncnorm(range, a = 0, b = 1, mean = dist_stats[1], sd = dist_stats[2])
dist <- dist/max(dist)

#plot it
par(mar = c(4, 4, 1, 1))
plot(range, dist, type = "l", xlab = "P", ylab = "Density", xlim = c(0, 1))
```

![](README_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

In order to simulate improvement in this distribution, we will randomly
draw a point sample and use it as the new mean of the normal
distribution in the next timestep. In order to “hit the brakes” on
improvement, the standard deviation of each new distribution will be the
`init_sd` multiplied by the new mean raised to the power of `brakes`.
For this first example we’ll use `init_sd` = 0.1 and `brakes` = 2.

``` r
#set seed
set.seed(123)

#set initial sd and brakes values
init_sd <- 0.1
brakes <- 2

#store distribution statistics (mean and sd)
dist_stats <- c(1, init_sd)

#begin plot
par(mar = c(4, 4, 1, 1))
plot(range, dist, type = "l", xlab = "P", ylab = "Density", xlim = c(0, 1))

#colors for plotting after each timepoint
colors <- rainbow((20-1)*1.25) #times 1.2 so it doesn't loop back around

#iterate
for(i in 1:20){
  #replace mean with point sample multiplied by improve
  dist_stats[1] <- truncnorm::rtruncnorm(1, a = 0, b = 1, mean = dist_stats[1], sd = dist_stats[2])
  
  #replace sd with initial sd multiplied by the new mean to the power of brakes
  dist_stats[2] <- init_sd*(dist_stats[1]^brakes)
  
  #create new dist to plot normalize it
  dist <- truncnorm::dtruncnorm(range, a = 0, b = 1, mean = dist_stats[1], sd = dist_stats[2])
  dist <- dist/max(dist)
  
  #plot new line
  lines(range, dist, type = "l", col = colors[i])
}
```

![](README_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

Now let’s set increase `init_sd` to 0.4 and `brakes` to 4 and see what
happens.

![](README_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

As you can see, we can get a wide range of outcomes by varying these two
simple parameters.

## Beta & Reference Times

Now, how do we use these distributions in combination with hold
sequences? First, let’s look at a diagram of the standardized speed
wall.

<img src="https://www.kindpng.com/picc/m/592-5929143_speed-climbing-wall-sketch-speed-climbing-route-map.png" width="65%" style="display: block; margin: auto;" />

As you can see, the standardized speed wall has a total of 20 hand holds
and 11 foot holds. Speed climbers often “smear” their feet on the wall
or use hand holds for feet, so we will only be modeling the time spent
on hand holds. Each climber will be initialized with a vector of their
beta, or a TRUE/FALSE for whether they use each hand hold in the route,
along with a vector of sequence ratios. The sequence ratios will be
drawn from a truncated normal distribution with a lower bound at 0, a
mean of 1, and a standard deviation parameter that controls the initial
variation in times across holds. The actual amount of time spent on each
hold, then, will be these sequence ratios multiplied by an initial
climbing time (let’s say 18 seconds) divided by the number of holds. For
now the beta vectors will start out as all TRUE, so that all climbers
start out using every hold in the route. Here is an example of how the
beta and sequence ratio vectors are initialized.

``` r
#set number of holds
n_holds <- 20

#set initial mean speed
init_time <- 18

#set probability of initial beta holds at 1 (all holds on the route)
beta_true_prob <- 1

#set parameter controlling the SD of sequence ratios
sd_multiplier <- 0.5

#initialize starting beta
beta <- sample(c(TRUE, FALSE), n_holds, prob = c(beta_true_prob, 1-beta_true_prob), replace = TRUE)

#initialize sequence ratios
seq_ratios <- truncnorm::rtruncnorm(n_holds, a = 0, mean = 1, sd = sd_multiplier)

#print the beta and climbing time vectors
beta
```

    ##  [1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ## [16] TRUE TRUE TRUE TRUE TRUE

``` r
(init_time/n_holds)*seq_ratios
```

    ##  [1] 1.3028065 1.2951601 1.2697115 1.2098881 1.1492629 0.8721397 0.7623168
    ##  [8] 0.7287880 0.5873819 0.8064372 0.3305716 1.8760302 1.4435829 0.3946011
    ## [15] 0.7187018 0.6900051 1.2509843 0.8624839 1.0139933 0.8871540

This `sd_multiplier` value of 0.5 generates a distribution of times per
hold that is similar to the example distribution in [Reveret et
al. (2020)](https://www.frontiersin.org/articles/10.3389/fpsyg.2020.02188/full)
(see below), and is close to the variation in times per hold observed in
lead climbing [(Seifert et al.,
2020)](https://www.tandfonline.com/doi/full/10.1080/14763141.2020.1830161).
A more recent study that estimated times per hold for two high-level
climbers in the 2019 IFSC World Cup suggests that an `sd_multiplier`
value of 0.33 may be more appropriate [(Pandurevic et al.,
2022)](https://www.mdpi.com/1424-8220/22/6/2251), but in exploratory
analyses the posterior for this parameter converged to 0.5. We will make
a simplifying assumption and use 0.5 for all of our simulations.

![](README_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

In each timestep, each climber has a certain `innov_prob` probability of
innovation. Innovations are changes to the beta of the route, in this
case flipping one of the booleans of the `beta` vector from TRUE to
FALSE (or vice versa). Not all boolean flips are possible. The parameter
`adj_poss` controls the number of successive FALSE booleans that are
allowed in the model (to simulate the fact that skipping three holds in
a row and still being able to complete the route is *extremely*
unlikely). When a boolean flip occurs (i.e. when a hold is added to or
dropped from the beta) the `seq_ratios` of the adjacent holds are
resampled (based on the mean initial climbing times of the population),
since the amount of time spent on the adjacent holds is dependent on the
presence of the added/dropped hold. When each climber innovates, they
compare the overall route time for their current beta with the innovated
beta, and if the innovated beta is better then they overwrite their
current beta. See `SpeedClimbingABM.R` for details.

## Social Learning

In each timestep, each climber also has a certain `learn_prob`
probability of copying the beta of another climber. When a climber
copies another climber, they only have access to the top `n_top` fastest
climbers in the current timestep. They “try out” all of the betas of the
fastest climbers that are different from their own beta, and if one of
them yields a faster time then they overwrite their current `beta` and
`seq_ratios` with those of better beta. See `SpeedClimbingABM.R` for
details.

## Dynamic Population Size

We also allow the population size to grow over time. For now we just
have a dynamic population size of `n` climbers and the probability
`leave_prob` of leaving the sport in each timestep. Basically, at the
end of each timestep a random subset of climbers leave the sport, and
new climbers enter the sport to bring the population up to the size
specified by `n`. The probability of a climber leaving the sport is a
function of how far their current record is from the current best among
climbers of the same gender as well as how long they’ve been in the
sport, so older climbers with slower times are more likely to leave.
Specifically it is their `current_record` divided by the minimum
`current_record` in the population multiplied by `age`, or the number of
timesteps that they have been in the sport. Each new climber that enters
the sport inherits information from two randomly sampled climbers from
that timestep: `seq_ratios` and `beta` come from one, and `dist` comes
from another. See `SpeedClimbingABM.R` for details.

## Interaction Effects

For our study, we want to add some interaction effects between
innovation, social learning, ranking, and population size. More
specifically, we want to innovation and social learning change with
climbers’ current rankings and with population size. To do this, we will
adjust the mean innovation and social learning rates of the population
with the following two equations for innovation and social learning,
respectively:

*μ* = *μ*<sub>*a**v**g*</sub> + (*μ*<sub>*a**v**g*</sub> \* *t*<sub>*n**o**r**m*</sub> \* *ϕ*) + (*μ*<sub>*a**v**g*</sub> \* *p*<sub>*n**o**r**m*</sub> \* *ω*)

*λ* = *λ*<sub>*a**v**g*</sub> + (*λ*<sub>*a**v**g*</sub> \* *t*<sub>*n**o**r**m*</sub> \* *ϵ*) + (*λ*<sub>*a**v**g*</sub> \* *p*<sub>*n**o**r**m*</sub> \* *σ*)

In each case there are three parameters in play. For innovation (*μ*)
there is the population average (*μ*<sub>*a**v**g*</sub>), the effect of
a climber’s record time (*ϕ*), and the effect of population size (*ω*).
The effects of time and population size are computed as follows: (1) all
values are collected and normalize to be between -0.5 and 0.5 (to
prevent *μ* from taking a value that is negative or greater than 1), (2)
the normalized values are multiplied by number between -1 and 1, where
the sign and absolute value control the direction and strength of the
effect, respectively, and (3) the values are the multiplied by the
innovation rate. Social learning is computed identically.

## Test Run

The ABM functions are in the `SpeedClimbingABM.R` file. Please refer to
the functions in this file for details.

``` r
source("SpeedClimbingABM.R")
```

Let’s do a test run of this model based on the observed data from the
women in the population.

``` r
load("data.RData")
data <- data[which(data$gender == "W"), ]
```

First we need the observed population sizes, leave probabilities, and
initial climbing times.

``` r
#get years
years <- sort(unique(data$year))

#population sizes
n <- unlist(lapply(1:length(years), function(x){nrow(data[which(data$year == years[x]), ])}))

#leave probabilities
leave_prob <- unlist(lapply(1:(length(years)-1), function(x){
  temp_a <- data[which(data$year == years[x]), ]
  temp_b <- data[which(data$year == years[x+1]), ]
  length(which(!(temp_a$athlete %in% temp_b$athlete)))/nrow(temp_a)
}))

#initial climbing times
init_times <- data[which(data$year == years[1]), ]$time

#get initial sd
init_sd <- sd(data$time[which(data$year == sort(unique(data$year))[1])]/mean(data$time[which(data$year == sort(unique(data$year))[1])]))
```

Now let’s run the model with an initial population size of 37, a 0.2
probability of learning from the top 20 climbers, and a 0.2 probability
of innovation where no more than 2 adjacent holds can be skipped. `bw`
and `ylim`, control the density bandwidth and y-axis limit in the plot,
respectively. `sum_stats` controls whether summary statistics are
calculated from the output, `plot` controls whether the output is
plotted, and `bw` and `ylim` control the density bandwidth and y-axis
limit of the plot, respectively.

``` r
#set seed
set.seed(1234)

#store starting time
start <- Sys.time()

#run model
output <- SpeedClimbingABM(n = n, leave_prob = leave_prob, init_times = init_times,
                           n_holds = 20, beta_true_prob = 1, learn_prob = 0.4, n_top = 20,
                           innov_prob = 0.2, adj_poss = 2, init_sd = 0.15, brakes = 4,
                           sum_stats = TRUE, plot = TRUE, bw = 0.6, ylim = 0.4)
```

![](README_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

``` r
#print run time
Sys.time() - start
```

    ## Time difference of 0.280817 secs

In the above plot, each distribution (from right to left) is the
distribution of `current_records` for climbers in the population in each
timestep. The output of this ABM a table of the summary statistics
(quantiles) from each timestep.

``` r
output
```

    ##              0%       10%       20%       30%      40%      50%      60%
    ##  [1,] 23.010000 24.152000 26.850000 30.834000 34.50600 38.38000 47.34800
    ##  [2,] 12.761826 20.023355 23.793427 26.518366 32.88718 38.12152 39.77302
    ##  [3,] 13.332129 16.102238 18.217846 20.026528 20.39239 22.17648 24.07452
    ##  [4,] 10.494364 13.488622 14.078820 14.839543 17.72227 17.98869 18.65948
    ##  [5,] 10.129161 11.875197 13.279061 13.859285 14.28961 14.91859 16.99825
    ##  [6,]  9.412327 10.737420 12.248010 13.613367 14.33061 14.75506 15.24896
    ##  [7,]  8.742520 10.330741 11.563589 12.618762 13.35801 13.79717 15.90748
    ##  [8,]  9.154391 10.753670 11.544902 12.224039 12.62664 13.34815 14.27216
    ##  [9,]  9.002533 10.616144 11.661566 12.261403 12.75051 13.03875 14.58703
    ## [10,]  8.849515 10.031542 10.878508 11.619035 12.07584 13.04801 13.41002
    ## [11,]  7.346235  9.358459 10.258824 10.827286 11.46269 12.33103 13.03682
    ## [12,]  6.999091  8.611335  9.335621  9.734487 10.07129 10.85267 11.47201
    ## [13,]  6.746023  8.396509  8.917587  9.690593 10.08004 10.46992 10.91139
    ##            70%      80%      90%     100%
    ##  [1,] 50.29400 53.99400 62.84000 86.39000
    ##  [2,] 41.68019 45.42380 52.22495 55.25121
    ##  [3,] 25.18877 31.36218 37.77011 54.06845
    ##  [4,] 23.79976 24.79142 27.17880 39.05654
    ##  [5,] 18.30585 20.39896 21.80457 29.96750
    ##  [6,] 17.25045 18.98560 21.45520 23.28590
    ##  [7,] 16.71615 18.22624 19.77208 21.24571
    ##  [8,] 15.60521 16.77040 18.29740 20.52947
    ##  [9,] 15.01069 16.74055 18.08729 22.98771
    ## [10,] 14.99778 16.59659 18.30505 23.69740
    ## [11,] 13.84865 15.02037 17.23025 22.32250
    ## [12,] 12.29134 12.88204 14.23525 21.69915
    ## [13,] 11.84014 12.78893 13.65159 19.37199

Let’s do another run of the model, but let’s add interaction effects so
that climbers with slower times are more likely to learn
(`learn_x_times` = 1), climbers with faster times are more likely to
innovate (`innov_x_times` = -1), and both learning and innovation become
more common as population size increases (`learn_x_pop` and
`innov_x_pop` = 1).

``` r
#set seed
set.seed(123)

#store starting time
start <- Sys.time()

#run model
output <- SpeedClimbingABM(n = n, leave_prob = leave_prob, init_times = init_times,
                           n_holds = 20, beta_true_prob = 1, learn_prob = 0, n_top = 20,
                           innov_prob = 0.2, adj_poss = 2, init_sd = 0.15, brakes = 4,
                           learn_x_times = 1, innov_x_times = -1, learn_x_pop = 1, innov_x_pop = 1,
                           sum_stats = TRUE, plot = TRUE, bw = 0.6, ylim = 0.4)
```

![](README_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

``` r
#print run time
Sys.time() - start
```

    ## Time difference of 0.085361 secs
