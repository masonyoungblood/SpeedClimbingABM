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

    ##  [1] 1.3410881 1.6738316 1.8399749 0.9471451 0.4868738 1.2733414 1.1009485
    ##  [8] 1.5056046 1.2812829 0.5276838 0.6359598 0.1962652 0.7314883 1.2037818
    ## [15] 1.4770855 1.6936720 0.8822396 1.2475899 0.7959805 0.7323823

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

Among all possible boolean flips, we will weight their probability by
the inverse of the Euclidean distance (in meters) between the adjacent
holds, raised to the power of a parameter `constraint`. In other words,
if skipping a hold means traversing a large distance, then that hold
will be less likely to be skipped. `constraint` controls the strength of
this physical constraint. The probability of skipping a hold will be a
normalized version of the following:

*P*(*s**k**i**p**p**i**n**g*) = (1/*d**i**s**t*)<sup>*c**o**n**s**t**r**a**i**n**t*</sup>

The first two holds can never be skipped, because every climber is
required to start on them. The Euclidean distance for the last hold will
be the distance between the second-to-last hold and the buzzer. Below is
a density plot showing the (yet to be normalized) probabilities of
skipping each hold under `constraint` = {0.5, 1, 2} (blue, black, and
red, respectively).

![](README_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

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
innovation, and `constraint` = 2. The `improve_rate_m` of athletic
improvement will be 2, the `improve_rate_sd` will be 0.2, and the `min`
improvement possible will be 0.4. `bw` and `ylim`, control the density
bandwidth and y-axis limit in the plot, respectively. `sum_stats`
controls whether summary statistics are calculated from the output,
`plot` controls whether the output is plotted, and `bw` and `ylim`
control the density bandwidth and y-axis limit of the plot,
respectively.

``` r
#store starting time
start <- Sys.time()

#run model
output <- SpeedClimbingABM(n = n, years = years, pop_data = pop_data, grid = grid,
                           n_holds = 20, beta_true_prob = 1, learn_prob = 0.2, n_top = 20,
                           innov_prob = 0.2, constraint = 2, improve_rate_m = 2, improve_rate_sd = 0.2, improve_min = 0.4,
                           sum_stats = TRUE, plot = TRUE, bw = 0.6, ylim = 0.4)
```

![](README_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

``` r
#print run time
Sys.time() - start
```

    ## Time difference of 0.222872 secs

In the above plot, each distribution (from right to left) is the
distribution of `current_records` for climbers in the population in each
timestep. The output of this ABM a table of the summary statistics
(quantiles) from each timestep.

``` r
output
```

    ##              0%       10%       20%       30%       40%       50%       60%
    ##  [1,] 15.400000 18.256000 20.408000 21.188000 24.432000 28.640000 29.712000
    ##  [2,]  8.080000 10.970601 12.238180 13.486824 16.354903 20.672619 24.314000
    ##  [3,]  5.379413  7.188737  9.213149 10.202282 11.874527 12.745415 13.939666
    ##  [4,]  5.182645  7.525000  8.506954  8.982154  9.690000 11.390000 12.370000
    ##  [5,]  4.829812  6.773774  7.209507  7.502846  8.000270  8.495000  8.756732
    ##  [6,]  4.668466  6.195000  6.843343  7.258000  7.620000  7.968129  8.758111
    ##  [7,]  4.790386  6.388059  6.837891  7.181000  7.563798  7.800000  8.846330
    ##  [8,]  5.091298  6.061813  6.430000  6.982200  7.550516  7.753000  8.554000
    ##  [9,]  5.068381  6.138628  6.327384  6.681232  7.036209  7.410000  7.918714
    ## [10,]  4.642779  5.630434  6.356000  6.799000  7.224000  7.605000  7.950000
    ## [11,]  4.505981  5.406971  6.009795  6.484111  6.869982  7.401279  8.128000
    ## [12,]  4.419796  5.461823  6.240868  6.696000  7.318000  7.691094  8.113927
    ## [13,]  4.264803  5.119903  5.735800  6.199191  6.744700  7.158822  7.643400
    ##             70%       80%      90%     100%
    ##  [1,] 32.472000 39.548000 43.25400 62.06000
    ##  [2,] 28.862000 31.218000 33.57986 48.08000
    ##  [3,] 17.466524 20.992011 26.70123 33.23000
    ##  [4,] 15.940227 21.570000 24.23000 31.90000
    ##  [5,] 10.483111 14.836487 18.47500 25.06290
    ##  [6,]  9.939494 11.530000 14.93640 23.80832
    ##  [7,] 10.934367 13.610000 17.26600 24.79000
    ##  [8,]  8.996000  9.802295 11.81320 24.51685
    ##  [9,]  8.410000  9.415867 11.26200 18.15000
    ## [10,]  8.654000 10.150000 12.54600 23.37000
    ## [11,]  8.635000  9.544000 10.86386 16.59000
    ## [12,]  8.760000  9.239800 10.44400 19.83000
    ## [13,]  8.091800  8.786233  9.74960 14.01663

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
                           innov_prob = 0.2, constraint = 2, improve_rate_m = 2, improve_rate_sd = 0.2, improve_min = 0.4,
                           learn_x_times = 1, innov_x_times = -1, learn_x_pop = 1, innov_x_pop = 1,
                           sum_stats = TRUE, plot = TRUE, bw = 0.6, ylim = 0.4)
```

![](README_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

``` r
#print run time
Sys.time() - start
```

    ## Time difference of 0.266813 secs

``` r
output
```

    ##              0%       10%       20%       30%       40%       50%       60%
    ##  [1,] 15.400000 18.256000 20.408000 21.188000 24.432000 28.640000 29.712000
    ##  [2,]  8.080000 10.987142 12.946643 14.194191 16.412899 20.810919 25.378000
    ##  [3,]  5.174844  7.532564  9.373042 10.892661 11.630945 12.650603 14.945497
    ##  [4,]  5.421619  7.308746  8.174333  9.045000  9.710990 11.390000 12.410000
    ##  [5,]  5.220905  6.537623  7.061955  7.325000  7.820000  8.520000  9.188215
    ##  [6,]  5.014731  5.418833  6.162348  6.666148  7.256000  7.690000  8.819247
    ##  [7,]  4.540514  5.208982  6.308512  6.980000  7.313723  7.812868  8.664763
    ##  [8,]  4.390889  5.084521  6.002863  6.928800  7.412009  7.741000  8.398352
    ##  [9,]  4.501665  5.036194  5.626756  5.951686  6.773039  7.600811  8.220000
    ## [10,]  4.496466  5.069759  5.991568  6.605000  7.224000  7.625000  7.994668
    ## [11,]  4.097029  5.025533  5.721323  6.083507  6.562987  7.218126  7.771216
    ## [12,]  3.712273  4.800847  5.562263  6.314000  6.680000  7.280000  7.840000
    ## [13,]  3.665448  4.486276  5.025545  5.569550  5.941799  6.250691  6.576176
    ##             70%       80%      90%     100%
    ##  [1,] 32.472000 39.548000 43.25400 62.06000
    ##  [2,] 28.862000 30.598000 33.01574 48.08000
    ##  [3,] 19.215290 20.402669 25.64417 33.23000
    ##  [4,] 16.016298 20.530000 23.29000 31.90000
    ##  [5,] 10.116920 13.312550 16.12117 18.52000
    ##  [6,]  9.238217 11.530000 13.60085 15.29905
    ##  [7,]  9.575216 13.546000 16.49700 24.79000
    ##  [8,]  8.888400  9.498400 11.79280 18.80665
    ##  [9,]  8.502000  8.790256 11.01986 18.15000
    ## [10,]  8.604000 10.018000 12.05100 23.37000
    ## [11,]  8.587000  9.504604 10.68966 16.59000
    ## [12,]  8.496400  9.206000 10.28300 19.83000
    ## [13,]  7.266805  8.109200  9.38600 13.30600
