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

    ##  [1] 0.8663235 1.0316852 0.3996837 0.6682919 0.8531268 0.3107566 1.6109695
    ##  [8] 1.3270237 0.1902318 0.9461157 1.0656838 0.7622323 0.3475074 0.3846847
    ## [15] 0.7474412 0.3104971 1.2578201 0.8464915 0.7891963 0.2259564

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

The first hold can never be skipped, because every climber is required
to start on it The Euclidean distance for the last hold will be the
distance between the second-to-last hold and the buzzer. `constraint` =
0 equalizes the probabilities, whereas `constraint` = 2 leads to a 10x
difference in the probability of innovating holds with the closest
vs. furthest adjacent holds. Below is a density plot showing the (yet to
be normalized) probabilities of skipping each hold under `constraint` =
{0.5, 1, 2} (blue, black, and red, respectively).

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
innovation, and `constraint` = 1. The `improve_rate_m` of athletic
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
                           innov_prob = 0.2, constraint = 1, improve_rate_m = 2, improve_rate_sd = 0.2, improve_min = 0.4,
                           sum_stats = TRUE, plot = TRUE, bw = 0.6, ylim = 0.4)
```

![](README_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

``` r
#print run time
Sys.time() - start
```

    ## Time difference of 0.3552511 secs

In the above plot, each distribution (from right to left) is the
distribution of `current_records` for climbers in the population in each
timestep. The output of this ABM a table of the summary statistics
(quantiles) from each timestep.

``` r
output
```

    ##              0%       10%       20%       30%       40%       50%       60%
    ##  [1,] 15.400000 18.256000 20.408000 21.188000 24.432000 28.640000 29.712000
    ##  [2,]  8.080000 10.951000 12.549569 13.851195 16.418657 19.633757 25.160000
    ##  [3,]  5.724069  7.273228  8.621227 10.249998 11.240000 13.219168 15.046140
    ##  [4,]  5.291820  7.296886  8.276121  9.045000  9.690000 11.390000 12.410000
    ##  [5,]  4.828431  6.569069  7.022928  7.325000  7.820000  8.475000  8.690000
    ##  [6,]  4.578932  6.090390  6.467761  7.050147  7.520000  7.947431  8.610000
    ##  [7,]  5.027576  5.968948  6.440297  7.089000  7.424000  7.805504  8.532677
    ##  [8,]  4.995140  5.826290  6.381385  6.932617  7.312000  7.741000  8.554000
    ##  [9,]  4.810681  5.820714  6.235133  6.646000  7.102216  7.593189  7.904107
    ## [10,]  4.702173  5.803929  6.478000  6.806599  7.230000  7.654218  8.006178
    ## [11,]  4.631086  5.767669  6.256013  6.605912  6.949718  7.442200  8.181147
    ## [12,]  4.620815  5.920000  6.402000  6.730837  7.285829  7.660000  8.181159
    ## [13,]  4.614992  5.482803  6.050427  6.594244  6.958668  7.364674  7.873232
    ##             70%       80%      90%     100%
    ##  [1,] 32.472000 39.548000 43.25400 62.06000
    ##  [2,] 28.862000 30.598000 33.13300 48.08000
    ##  [3,] 17.504887 23.894206 26.11259 33.23000
    ##  [4,] 16.979656 22.471647 24.23000 31.90000
    ##  [5,] 10.952123 14.900000 18.10706 21.56951
    ##  [6,] 10.550316 11.824156 16.20297 20.45811
    ##  [7,] 10.764312 13.546000 17.26600 24.79000
    ##  [8,]  8.974000  9.778831 11.86200 19.92865
    ##  [9,]  8.324000  9.421237 11.26200 18.15000
    ## [10,]  8.604000 10.150000 12.05100 23.37000
    ## [11,]  8.689236  9.770639 11.00600 16.59000
    ## [12,]  8.772000  9.300400 10.44400 19.83000
    ## [13,]  8.343955  8.880738 10.04744 14.02597

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
                           innov_prob = 0.2, constraint = 1, improve_rate_m = 2, improve_rate_sd = 0.2, improve_min = 0.4,
                           learn_x_times = 1, innov_x_times = -1, learn_x_pop = 1, innov_x_pop = 1,
                           sum_stats = TRUE, plot = TRUE, bw = 0.6, ylim = 0.4)
```

![](README_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

``` r
#print run time
Sys.time() - start
```

    ## Time difference of 0.4214411 secs

``` r
output
```

    ##              0%       10%       20%       30%       40%       50%       60%
    ##  [1,] 15.400000 18.256000 20.408000 21.188000 24.432000 28.640000 29.712000
    ##  [2,]  8.080000 11.204119 12.480000 14.173174 16.848000 21.267532 25.378000
    ##  [3,]  5.672369  8.982705  9.717226 10.746943 11.838823 14.009765 16.283206
    ##  [4,]  5.494154  7.790000  8.505359  9.130000 10.328315 11.505000 12.816501
    ##  [5,]  4.472480  6.772648  7.230000  7.719526  8.166333  8.590000  9.317773
    ##  [6,]  4.294320  5.997802  6.811376  7.258000  7.520000  8.312918  8.829000
    ##  [7,]  4.791437  6.226965  6.689193  7.098220  7.296322  7.800000  8.780000
    ##  [8,]  4.368213  6.015159  6.622849  7.007596  7.462000  7.741000  8.473214
    ##  [9,]  4.319913  5.934286  6.265114  6.646000  7.092659  7.465291  8.078863
    ## [10,]  4.693393  5.904901  6.403796  6.671904  7.224000  7.605000  7.844000
    ## [11,]  4.227967  5.693922  6.019117  6.467862  6.762540  7.510000  8.124955
    ## [12,]  3.742245  5.084202  5.920000  6.488492  6.890000  7.504000  7.966000
    ## [13,]  3.675545  4.513036  5.346854  5.784922  6.331800  6.837637  7.352410
    ##             70%       80%      90%     100%
    ##  [1,] 32.472000 39.548000 43.25400 62.06000
    ##  [2,] 28.557000 30.598000 33.13300 48.08000
    ##  [3,] 17.860299 22.135771 25.52600 33.23000
    ##  [4,] 15.897763 21.570000 23.69500 31.90000
    ##  [5,] 11.300717 14.545354 19.42523 22.69265
    ##  [6,] 10.242858 12.247586 16.74478 19.24696
    ##  [7,] 10.772664 13.546000 17.30500 24.79000
    ##  [8,]  8.996000  9.734000 11.86200 19.38308
    ##  [9,]  8.418375  9.370154 11.38325 18.15000
    ## [10,]  8.500907 10.018000 12.54600 23.37000
    ## [11,]  8.822877  9.572000 11.00600 16.59000
    ## [12,]  8.530000  9.206000 10.32436 19.83000
    ## [13,]  7.911299  8.438150  9.66340 13.30600
