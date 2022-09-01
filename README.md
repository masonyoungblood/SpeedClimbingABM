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
rates <- truncnorm::rtruncnorm(100, a = 1, mean = 2.5, sd = 0.3)
y <- sapply(1:length(rates), function(h){bounded_exp(x, rates[h], 0.07)})

#plot
par(mar = c(4, 4, 1, 1))
matplot(x, y, type = "l", xlab = "Timestep", ylab = "Athletic Improvement", ylim = c(0, 1), col = scales::alpha("black", 0.2), lty = 1)
```

![](README_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

In each timestep, a climbers current record comes from multiplying the
time per handhold by this `athletic_improvement` index. After the first
timestep, new climbers entering the population will have their
distribution of `athletic_improvement` indices divided by the
`athletic_improvement` value from the timestep that they enter. In this
way, a climber who enters the population with a lower climbing time will
follow the same general improvement trajectory as everyone else in the
population.

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
reference climbing time per hold (initially their starting climbing time
divided by the number of holds). For now the beta vectors will start out
as all TRUE, so that all climbers start out using every hold in the
route. Here is an example of how the beta and sequence ratio vectors are
initialized.

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

    ##  [1] 1.2016715 0.3635775 1.0146782 0.8979425 1.2744094 0.4522616 0.8356682
    ##  [8] 1.0342582 0.9435056 0.7262213 1.6958151 1.4004775 1.0410451 0.4203595
    ## [15] 0.9260607 0.3129652 0.4928775 0.8009349 1.3722188 1.0508054

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
FALSE. Not all boolean flips are possible. The parameter `adj_poss`
controls the number of successive FALSE booleans that are allowed in the
model. We will assume `adj_poss` = 2, because (to our knowledge) nobody
has successfully skipped three or more holds in speed climbing.

Additionally, we will include a parameter `max_dist` that controls the
maximum travel distance (Euclidean distance in meters) that a climber
can travel after skipping a hold. The prior for this parameter will go
from the minimum distance between holds that are one hold apart to the
maximum distance between holds that are two (`adj_poss`) holds apart.
The first hold can never be skipped, because every climber is required
to start on it. The Euclidean distance for the last hold will be the
distance between the second-to-last (or third-to-last) hold and the
buzzer.

Below are the distributions of distances between holds that are one and
two holds apart …

![](README_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-6-2.png)<!-- -->

Among the available holds within `max_dist`, the probability of
innovating a hold will be weighted by the inverse of the distance
between the adjacent holds, raised to the power of a parameter
`constraint`. In other words, if skipping a hold means traversing a
larger distance, then that hold will be less likely to be skipped.
`constraint` controls the strength of this physical constraint.

When a boolean flip occurs (i.e. when a hold is added to or dropped from
the beta) the `seq_ratios` of the adjacent holds are resampled (based on
the mean initial climbing times of the population), since the amount of
time spent on the adjacent holds is dependent on the presence of the
added/dropped hold. When each climber innovates, they compare the
overall route time for their current beta with the innovated beta, and
if the innovated beta is better then they overwrite their current beta.
See `SpeedClimbingABM.R` for details.

## Social Learning

In each timestep, each climber also has a certain `learn_prob`
probability of copying the beta of another climber. When a climber
copies another climber, they only have access to the top `n_top` fastest
climbers in the current timestep. They “try out” all of the betas of the
fastest climbers that are different from their own beta, and if one of
them yields a faster time then they overwrite their current `beta` and
`seq_ratios` with those of better beta. See `SpeedClimbingABM.R` for
details.

## Population Size & Turnover

In the original version of the model we did not use all available
information about when specific climbers entered the sport, left the
sport, etc. This version was far too stochastic to use for inference, so
we are now incorporating this information explicitly into the model.
First let’s take a look at the data:

    ##       athlete gender   time year
    ##    1:    1474      W 25.170 2007
    ##    2:      42      W 23.450 2007
    ##    3:      71      W 50.060 2007
    ##    4:      46      W 23.260 2007
    ##    5:    4735      W 48.950 2007
    ##   ---                           
    ## 2007:    4206      W 15.058 2019
    ## 2008:    2651      W 15.494 2019
    ## 2009:    4215      W 18.339 2019
    ## 2010:    2636      W 20.085 2019
    ## 2011:    2625      W 21.796 2019

Each row is an climber in a particular year with their current record.
When a climber first enters the sport we will initialize them with their
current record in that year, and then we will simulate innovation,
learning, and athletic improvement from that baseline until they
eventually leave the sport. In the cases when climbers’ have gaps in
their careers we will treat them as separate cascades by re-initializing
their current record when they re-enter the sport and simulating change
from there. Climbers who only appear in the dataset once be initialized
with their current record in the year they appear and will disappear in
the next year with no simulated improvement (see below).

![](README_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

To do this, let’s create a new `pop_data` object that includes all of
the information we need. Each row corresponds to a continuous sequence
of years: the first column is the climbers’ ID, the second is their
starting year, the third is their end year, and the fourth is their best
time in the starting year.

    ##         ID start  end   time
    ##    1: 1474  2007 2009 25.170
    ##    2:   42  2007 2007 23.450
    ##    3:   71  2007 2008 50.060
    ##    4:   46  2007 2009 23.260
    ##    5:   46  2011 2011 11.770
    ##   ---                       
    ## 1063: 4206  2019 2019 15.058
    ## 1064: 2651  2019 2019 15.494
    ## 1065: 4215  2019 2019 18.339
    ## 1066: 2636  2019 2019 20.085
    ## 1067: 2625  2019 2019 21.796

After the first timestep, each new climber that joins the population
will be initialized with a random `seq_ratios` and `beta` drawn from an
existing climber. See `SpeedClimbingABM.R` for details.

## Interaction Effects

For our study, we want to add some interaction effects between
innovation, social learning, ranking, and population size. More
specifically, we want to innovation and social learning change with
climbers’ current rankings and with population size. To do this, we will
adjust the mean innovation and social learning rates of the population
with the following two equations for innovation and social learning,
respectively:

*l**o**g**i**t*(*μ*) = *l**o**g**i**t*(*μ*<sub>*a**v**g*</sub>) + *t*<sub>*s**c**a**l**e*</sub> \* *ϕ* + *p*<sub>*s**c**a**l**e*</sub> \* *ω*

*l**o**g**i**t*(*λ*) = *l**o**g**i**t*(*λ*<sub>*a**v**g*</sub>) + *t*<sub>*s**c**a**l**e*</sub> \* *ϵ* + *p*<sub>*s**c**a**l**e*</sub> \* *σ*

The calculations are done on the logit scale so that everything remains
between 0 and 1. In each case there are three parameters in play. For
innovation (*μ*) there is the population average
(*μ*<sub>*a**v**g*</sub>), the effect of a climber’s record time (*ϕ*),
and the effect of population size (*ω*). The effects of time and
population size are computed as follows: (1) all values are scaled and
centered with a mean of 0 and standard deviation of 1, (2) the scale
values are multiplied by number between -1 and 1, where the sign and
absolute value control the direction and strength of the effect,
respectively, and (3) the values are the multiplied by the
logit-transformed innovation rate. Social learning is computed
identically.

## Test Run

The ABM functions are in the `SpeedClimbingABM.R` file. Please refer to
the functions in this file for details.

``` r
source("SpeedClimbingABM.R")
```

Let’s do a test run of this model based on the observed data from the
men in the population.

``` r
load("data.RData")
pop_data <- pop_data[which(data$gender[match(pop_data$ID, data$athlete)] == "M"), ]
data <- data[which(data$gender == "M"), ]
```

First we need the observed population sizes, leave probabilities, and
initial climbing times.

``` r
#get years
years <- sort(unique(data$year))

#population sizes
n <- unlist(lapply(1:length(years), function(x){nrow(data[which(data$year == years[x]), ])}))
```

Now let’s run the model with an initial population size of 53, a 0.2
probability of learning from the top 20 climbers, a 0.2 probability of
innovation, `max_dist` = 2, and `constraint` = 1. The `improve_rate_m`
of athletic improvement will be 2, the `improve_rate_sd` will be 0.2,
and the `min` improvement possible will be 0.4. `bw` and `ylim`, control
the density bandwidth and y-axis limit in the plot, respectively.
`sum_stats` controls whether summary statistics are calculated from the
output, `plot` controls whether the output is plotted, and `bw` and
`ylim` control the density bandwidth and y-axis limit of the plot,
respectively.

``` r
#store starting time
start <- Sys.time()

#run model
output <- SpeedClimbingABM(n = n, years = years, pop_data = pop_data, grid = grid,
                           n_holds = 20, beta_true_prob = 1, learn_prob = 0.2, n_top = 20,
                           innov_prob = 0.2, max_dist = 2, constraint = 1, improve_rate_m = 2, improve_rate_sd = 0.2, improve_min = 0.4,
                           sum_stats = TRUE, plot = TRUE, bw = 0.6, ylim = 0.4)
```

![](README_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

``` r
#print run time
Sys.time() - start
```

    ## Time difference of 0.3574409 secs

In the above plot, each distribution (from right to left) is the
distribution of `current_records` for climbers in the population in each
timestep. The output of this ABM a table of the summary statistics
(quantiles) from each timestep.

``` r
output
```

    ##              0%       10%       20%       30%       40%       50%       60%
    ##  [1,] 15.400000 18.256000 20.408000 21.188000 24.432000 28.640000 29.712000
    ##  [2,]  8.080000 11.212208 13.006475 14.087226 15.603768 20.951223 24.921408
    ##  [3,]  5.004672  9.922364 10.685502 11.247696 11.614398 13.203548 15.065357
    ##  [4,]  4.565625  7.884936  8.773138  9.237135  9.880000 11.505000 12.690000
    ##  [5,]  4.197945  7.030000  7.335178  7.835000  8.229885  8.590000  9.030000
    ##  [6,]  4.011660  6.343744  7.120000  7.272186  7.641014  8.171768  8.610000
    ##  [7,]  5.401533  6.077252  6.475737  6.980000  7.234000  7.788274  8.534687
    ##  [8,]  4.790783  5.780095  6.372643  6.968657  7.561200  7.758000  8.562361
    ##  [9,]  4.609093  5.713902  6.349991  6.646000  6.915274  7.410000  8.092000
    ## [10,]  4.413766  5.771984  6.356000  6.659423  7.104361  7.640000  8.064000
    ## [11,]  4.407839  5.648615  6.146090  6.568001  7.032888  7.510000  8.124978
    ## [12,]  4.404906  5.619296  6.460000  6.817600  7.318000  7.711564  8.176606
    ## [13,]  4.380829  5.519397  5.958194  6.570256  7.049352  7.525000  7.974382
    ##             70%       80%      90%     100%
    ##  [1,] 32.472000 39.548000 43.25400 62.06000
    ##  [2,] 28.557000 30.598000 33.13300 48.08000
    ##  [3,] 17.331187 21.293987 27.09493 33.23000
    ##  [4,] 16.642057 22.180000 24.34453 31.90000
    ##  [5,] 10.754062 12.844932 18.39930 24.74039
    ##  [6,] 10.134444 11.530000 14.71543 23.37186
    ##  [7,] 11.686735 13.546000 17.35353 24.79000
    ##  [8,]  8.974000  9.786211 11.86200 22.74628
    ##  [9,]  8.427145  9.434648 10.64612 18.15000
    ## [10,]  8.558119 10.001742 11.95917 23.37000
    ## [11,]  8.635860  9.745827 10.87693 16.59000
    ## [12,]  8.852000  9.496408 10.52140 19.83000
    ## [13,]  8.468800  9.178912 10.28040 13.30600

Let’s do another run of the model, but let’s add interaction effects so
that climbers with slower times are more likely to learn
(`learn_x_times` = 1), climbers with faster times are more likely to
innovate (`innov_x_times` = -1), and both learning and innovation become
more common as population size increases (`learn_x_pop` and
`innov_x_pop` = 1).

``` r
#store starting time
start <- Sys.time()

#run model
output <- SpeedClimbingABM(n = n, years = years, pop_data = pop_data, grid = grid,
                           n_holds = 20, beta_true_prob = 1, learn_prob = 0.2, n_top = 20,
                           innov_prob = 0.2, max_dist = 2, constraint = 1, improve_rate_m = 2, improve_rate_sd = 0.2, improve_min = 0.4,
                           learn_x_times = 1, innov_x_times = -1, learn_x_pop = 1, innov_x_pop = 1,
                           sum_stats = TRUE, plot = TRUE, bw = 0.6, ylim = 0.4)
```

![](README_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

``` r
#print run time
Sys.time() - start
```

    ## Time difference of 0.4235671 secs

``` r
output
```

    ##              0%       10%       20%       30%       40%       50%       60%
    ##  [1,] 15.400000 18.256000 20.408000 21.188000 24.432000 28.640000 29.712000
    ##  [2,]  8.080000 11.070635 12.884395 14.486272 16.128637 20.618114 24.783051
    ##  [3,]  5.524309  8.498236  9.564546 11.008008 11.880000 13.036565 15.098352
    ##  [4,]  4.745187  7.503882  8.649010  9.216619 10.260000 11.346142 12.410000
    ##  [5,]  3.959606  6.877707  7.290000  7.745000  8.263424  8.590000  9.030000
    ##  [6,]  3.830641  5.344099  6.270000  7.188000  7.535021  7.946475  8.610000
    ##  [7,]  4.016897  5.931408  6.498952  6.980000  7.407757  7.800000  8.541838
    ##  [8,]  3.824089  5.003719  6.195525  6.928800  7.426266  7.741000  8.292817
    ##  [9,]  3.674447  4.481168  5.406298  6.216059  6.658000  7.380000  8.049877
    ## [10,]  3.321633  4.587125  5.834472  6.640419  7.230000  7.645000  8.028000
    ## [11,]  3.482513  4.506797  5.290784  5.946015  6.529663  7.185000  7.738000
    ## [12,]  2.713401  4.583468  5.416347  6.314000  6.816800  7.435000  7.838000
    ## [13,]  2.620836  3.710568  4.464717  4.961161  5.377142  5.916462  6.604074
    ##             70%       80%      90%     100%
    ##  [1,] 32.472000 39.548000 43.25400 62.06000
    ##  [2,] 28.862000 30.598000 33.13300 56.35012
    ##  [3,] 19.921562 21.905738 25.52600 46.32093
    ##  [4,] 17.445147 21.570000 23.75463 31.90000
    ##  [5,] 10.404615 13.125095 17.62395 20.69978
    ##  [6,]  9.115920 10.403819 12.34282 17.20412
    ##  [7,]  9.399614 13.546000 16.79027 24.79000
    ##  [8,]  8.850100  9.403000 11.79280 17.89650
    ##  [9,]  8.327514  8.712000 10.11400 18.15000
    ## [10,]  8.604000  9.648102 11.82700 23.37000
    ## [11,]  8.439186  9.342000 10.10109 16.59000
    ## [12,]  8.496400  9.180403 10.25600 19.83000
    ## [13,]  7.381000  8.109200  9.33540 13.30600
