Speed Climbing ABM
================
Mason Youngblood

<img src="https://helios-i.mashable.com/imagery/articles/03XKxYMwkNiZJGiHUu2FWnw/hero-image.fill.size_1248x702.v1628089817.png" width="70%" style="display: block; margin: auto;" />

This agent-based model simulates a dynamic population of professional
speed climbers, and incorporates parameters for athletic improvement,
innovation of “beta” (or route sequence), and copying of other climbers’
beta.

## Athletic Improvement

Athletic improvement is simulated using a bounded exponential function
controlled by three parameters: `rate_m`, `rate_sd`, and `min`. `rate_m`
controls the mean exponential rate across all climbers, `rate_sd`
controls the standard deviation of the exponential rate to introduce
variation in ability, and `min` controls the asymptotic lower bound. The
result of this function is an `athletic_improvement` index for each
climber over time (x-axis). Below is an example with a `rate_m` of 2, a
`rate_sd` of 0.2, and a `min` of 0.4.

``` r
#bounded exponential function
bounded_exp <- function(x, rate, min){
  return((1-min)*(rate/rate^x)+min)
}

#generate values
x <- 1:12
rates <- truncnorm::rtruncnorm(100, a = 1, mean = 2, sd = 0.2)
y <- sapply(1:length(rates), function(h){bounded_exp(x, rates[h], 0.4)})

#plot
par(mar = c(4, 4, 1, 1))
matplot(x, y, type = "l", xlab = "Timestep", ylab = "Athletic Improvement", ylim = c(0, 1), col = scales::alpha("black", 0.2), lty = 1)
```

![](README_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

In each timestep, a climbers current record comes from multiplying the
time per handhold by this `athletic_improvement` index.

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

    ##  [1] 1.1710974 0.2207594 1.4080416 0.2674713 0.9134322 0.8957299 0.5969995
    ##  [8] 0.7342470 0.6648977 0.8517578 0.7224867 1.2618807 0.8806159 0.5630386
    ## [15] 0.6412850 1.1541968 0.6629331 0.9719648 0.7169040 1.3648290

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

![](README_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

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
```

Now let’s run the model with an initial population size of 37, a 0.2
probability of learning from the top 20 climbers, and a 0.2 probability
of innovation where no more than 2 adjacent holds can be skipped. The
`improve_rate_m` of athletic improvement will be 2, the
`improve_rate_sd` will be 0.2, and the `min` improvement possible will
be 0.4. `bw` and `ylim`, control the density bandwidth and y-axis limit
in the plot, respectively. `sum_stats` controls whether summary
statistics are calculated from the output, `plot` controls whether the
output is plotted, and `bw` and `ylim` control the density bandwidth and
y-axis limit of the plot, respectively.

``` r
#set seed
#set.seed(1234)

#store starting time
start <- Sys.time()

#run model
output <- SpeedClimbingABM(n = n, leave_prob = leave_prob, init_times = init_times,
                           n_holds = 20, beta_true_prob = 1, learn_prob = 0.2, n_top = 20,
                           innov_prob = 0.2, adj_poss = 2, improve_rate_m = 2, improve_rate_sd = 0.2, improve_min = 0.4,
                           sum_stats = TRUE, plot = TRUE, bw = 0.6, ylim = 0.4)
```

![](README_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

``` r
#print run time
Sys.time() - start
```

    ## Time difference of 0.2062211 secs

In the above plot, each distribution (from right to left) is the
distribution of `current_records` for climbers in the population in each
timestep. The output of this ABM a table of the summary statistics
(quantiles) from each timestep.

``` r
output
```

    ##              0%       10%       20%       30%       40%       50%       60%
    ##  [1,] 23.010000 24.152000 26.850000 30.834000 34.506000 38.380000 47.348000
    ##  [2,] 12.317890 14.875258 18.679820 19.697300 21.596174 22.835628 24.150158
    ##  [3,]  9.602480 11.740445 13.993279 15.923207 16.882375 17.459200 18.091515
    ##  [4,]  9.048284  9.421418  9.987964 12.500537 13.097471 13.562483 14.118041
    ##  [5,]  7.274149  8.068610  8.532725  8.763266  9.122720  9.357439 11.827960
    ##  [6,]  6.526504  7.010012  7.530442  7.770917  7.904872  8.349489  8.690358
    ##  [7,]  6.341333  6.574879  6.697899  6.811523  6.999505  7.424195  8.092159
    ##  [8,]  6.101528  6.476853  6.588419  6.704977  6.797921  7.286760  7.399286
    ##  [9,]  5.756849  5.974566  6.417760  6.516267  6.762173  7.173239  7.332569
    ## [10,]  5.237934  5.767548  6.397673  6.495518  7.064164  7.314365  7.322926
    ## [11,]  4.446071  4.677941  4.962160  4.980171  6.024904  6.533448  6.789394
    ## [12,]  4.259187  4.668821  4.790820  4.957376  4.963649  5.225461  6.015705
    ## [13,]  4.008788  4.257248  4.445312  4.667748  4.955430  4.964028  5.728012
    ##             70%       80%       90%     100%
    ##  [1,] 50.294000 53.994000 62.840000 86.39000
    ##  [2,] 25.108856 30.691648 41.519839 55.98966
    ##  [3,] 18.547126 19.575823 25.709494 41.33017
    ##  [4,] 15.039039 15.760290 17.815326 22.40885
    ##  [5,] 13.033629 13.806240 14.451229 15.67023
    ##  [6,] 11.854670 12.440098 13.318386 14.32805
    ##  [7,]  8.443771 11.471929 12.102746 13.19669
    ##  [8,]  8.032530  9.619175 11.078388 12.10552
    ##  [9,]  7.347140  8.189635 11.382328 12.02559
    ## [10,]  7.344286  7.446167 11.287757 11.49439
    ## [11,]  7.312783  7.328007  8.141557 11.40845
    ## [12,]  6.569899  7.086993  7.318022 11.32464
    ## [13,]  6.017425  6.712825  7.307341 11.27231

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
                           n_holds = 20, beta_true_prob = 1, learn_prob = 0.2, n_top = 20,
                           innov_prob = 0.2, adj_poss = 2, improve_rate_m = 2, improve_rate_sd = 0.2, improve_min = 0.4,
                           learn_x_times = 1, innov_x_times = -1, learn_x_pop = 1, innov_x_pop = 1,
                           sum_stats = TRUE, plot = TRUE, bw = 0.6, ylim = 0.4)
```

![](README_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

``` r
#print run time
Sys.time() - start
```

    ## Time difference of 0.1122799 secs
