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

    ##  [1] 0.65229230 0.85017774 0.58124922 1.25796814 1.18137493 0.63833767
    ##  [7] 1.25367652 1.65471092 0.01098751 1.03753524 0.39950981 0.33821582
    ## [13] 1.49153421 0.89728872 1.55949420 0.60290737 1.64239926 0.29313278
    ## [19] 0.86920733 0.93018777

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
FALSE. Not all boolean flips are possible. We will include a parameter
`max_dist` that controls the maximum travel distance (Euclidean distance
in meters) that a climber can travel after skipping a hold. The prior
for this parameter will go from the minimum distance between holds that
are one hold apart to the maximum distance between holds that are two
holds apart To our knowledge, nobody has successfully skipped three or
more holds in speed climbing, so this is an appropriately wide prior.
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

    ## Time difference of 0.1917951 secs

In the above plot, each distribution (from right to left) is the
distribution of `current_records` for climbers in the population in each
timestep. The output of this ABM a table of the summary statistics
(quantiles) from each timestep.

``` r
output
```

    ##              0%       10%       20%       30%       40%       50%       60%
    ##  [1,] 15.400000 18.256000 20.408000 21.188000 24.432000 28.640000 29.712000
    ##  [2,]  8.080000 11.713221 13.230382 14.825509 16.848000 20.681732 24.314000
    ##  [3,]  6.137332  8.464931  9.977822 11.096318 11.840509 13.474007 16.231523
    ##  [4,]  5.332223  7.670721  8.480000  9.265000 10.313183 11.505000 12.690000
    ##  [5,]  4.957141  7.058353  7.430000  7.847544  8.460000  8.705862  9.299056
    ##  [6,]  4.782399  6.834947  7.166333  7.363209  7.766841  8.243902  8.829000
    ##  [7,]  4.641932  6.657534  6.884983  7.181000  7.636000  8.219211  8.822027
    ##  [8,]  4.034412  6.138800  6.672423  7.028000  7.619031  7.849000  8.634020
    ##  [9,]  3.924371  5.272891  5.832917  6.544156  6.864190  7.380000  8.030181
    ## [10,]  3.917408  5.019389  5.627358  6.562000  7.224000  7.640000  8.032009
    ## [11,]  3.446072  4.691355  5.563178  6.186535  6.643469  7.310000  7.881432
    ## [12,]  3.589976  5.110404  5.957216  6.519785  6.826800  7.600000  8.016130
    ## [13,]  3.586636  4.726747  5.386970  5.838980  6.256800  6.792751  7.371800
    ##             70%       80%       90%     100%
    ##  [1,] 32.472000 39.548000 43.254000 62.06000
    ##  [2,] 28.862000 31.218000 34.151000 48.08000
    ##  [3,] 20.206797 24.435371 25.634489 33.23000
    ##  [4,] 17.612979 22.129083 23.862902 31.90000
    ##  [5,] 10.715717 14.900000 19.936331 22.66092
    ##  [6,] 10.380000 12.830000 18.586290 21.46823
    ##  [7,] 11.749000 13.546000 18.425050 24.79000
    ##  [8,]  9.026400  9.734000 11.862000 24.97120
    ##  [9,]  8.396859  8.960445 10.675733 18.26773
    ## [10,]  8.604000 10.150000 12.051000 23.37000
    ## [11,]  8.558815  9.640000 10.841409 16.59000
    ## [12,]  8.772000  9.413400 10.503000 19.83000
    ## [13,]  8.008068  8.788600  9.975361 13.59109

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

    ## Time difference of 0.217845 secs

``` r
output
```

    ##              0%       10%       20%       30%       40%       50%       60%
    ##  [1,] 15.400000 18.256000 20.408000 21.188000 24.432000 28.640000 29.712000
    ##  [2,]  8.080000 11.562981 12.479813 14.172841 16.168967 20.296796 25.160000
    ##  [3,]  6.043035  8.277218  9.775976 11.075495 11.700462 12.779405 15.373980
    ##  [4,]  5.228193  7.599383  8.469579  9.130000 10.104529 11.469959 12.690000
    ##  [5,]  4.830081  6.933944  7.357581  7.835561  8.490000  8.739531  9.508784
    ##  [6,]  4.635572  6.263508  7.120000  7.323120  7.620000  8.382609  8.960000
    ##  [7,]  5.236650  6.380862  6.938383  7.228073  7.424000  7.944887  8.633391
    ##  [8,]  5.195750  6.263138  6.875348  7.229287  7.684000  7.958395  8.506301
    ##  [9,]  5.684568  5.948669  6.545060  6.870749  7.153807  7.900000  8.220000
    ## [10,]  5.476944  6.338154  6.640346  6.961427  7.242036  7.645000  8.064000
    ## [11,]  4.853170  6.039679  6.480162  6.835541  7.278000  7.633688  8.275006
    ## [12,]  3.632268  5.294434  6.023574  6.580369  7.089442  7.620000  8.014000
    ## [13,]  3.050733  5.012867  5.620358  6.109421  6.489192  6.860000  7.555770
    ##             70%       80%      90%     100%
    ##  [1,] 32.472000 39.548000 43.25400 62.06000
    ##  [2,] 28.889246 30.598000 33.13300 48.08000
    ##  [3,] 18.662874 23.207538 25.52600 33.23000
    ##  [4,] 17.203950 21.710000 23.69500 31.90000
    ##  [5,] 10.900073 14.400619 19.41311 21.56242
    ##  [6,] 10.318678 12.359520 16.12290 19.46955
    ##  [7,] 11.714826 13.539249 17.48478 24.79000
    ##  [8,]  8.993197  9.734000 12.24469 21.18766
    ##  [9,]  8.502000  9.110170 11.21308 18.15000
    ## [10,]  8.604000 10.150000 12.54600 23.37000
    ## [11,]  9.158000  9.782000 11.96911 18.74039
    ## [12,]  8.772000  9.239800 10.50300 19.83000
    ## [13,]  8.117585  8.675284  9.66654 13.30600
