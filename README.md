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

# #generate values
# x <- 1:12
# rates <- truncnorm::rtruncnorm(100, a = 1, mean = 2, sd = 0.2)
# y <- sapply(1:length(rates), function(h){bounded_exp(x, rates[h], 0.4)})

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

    ##  [1] 0.7020366 0.7505131 1.5601573 0.4376044 0.7571112 0.8110028 0.1830058
    ##  [8] 0.7243347 0.9740226 2.0312848 1.6003573 0.8397763 1.2809386 1.3197038
    ## [15] 1.1861912 1.1418648 0.8981191 1.3700834 1.4813591 1.6863201

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

![](README_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

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
#store starting time
start <- Sys.time()

#run model
output <- SpeedClimbingABM(n = n, years = years, pop_data = pop_data,
                           n_holds = 20, beta_true_prob = 1, learn_prob = 0.2, n_top = 20,
                           innov_prob = 0.2, adj_poss = 2, improve_rate_m = 2, improve_rate_sd = 0.2, improve_min = 0.4,
                           sum_stats = TRUE, plot = TRUE, bw = 0.6, ylim = 0.4)
```

![](README_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

``` r
#print run time
Sys.time() - start
```

    ## Time difference of 0.3081088 secs

In the above plot, each distribution (from right to left) is the
distribution of `current_records` for climbers in the population in each
timestep. The output of this ABM a table of the summary statistics
(quantiles) from each timestep.

``` r
output
```

    ##              0%       10%       20%       30%       40%       50%       60%
    ##  [1,] 15.400000 18.256000 20.408000 21.188000 24.432000 28.640000 29.712000
    ##  [2,]  8.040347 10.815672 12.103209 13.581734 15.939330 21.280104 24.554723
    ##  [3,]  5.043516  7.754359  9.100716 10.227195 11.880000 13.702529 16.552613
    ##  [4,]  4.353751  7.330000  8.280000  9.045000  9.880000 11.505000 12.850000
    ##  [5,]  4.006596  6.770278  7.100000  7.376533  7.980000  8.475000  8.690000
    ##  [6,]  3.831875  5.581637  6.270000  7.095559  7.358773  7.817781  8.362237
    ##  [7,]  4.948169  5.755661  6.426130  6.888000  7.234000  7.629122  8.308000
    ##  [8,]  4.798627  5.652986  6.275270  6.766000  7.295516  7.700000  8.138765
    ##  [9,]  4.742585  5.527332  5.803689  6.242308  6.600000  7.138911  7.643753
    ## [10,]  4.710918  5.538128  6.197870  6.565593  7.105216  7.595000  7.844000
    ## [11,]  4.693024  5.524281  6.007117  6.333105  6.679472  7.275000  7.942715
    ## [12,]  4.526985  5.523528  6.267066  6.587757  7.108947  7.620000  7.993465
    ## [13,]  4.322934  5.442188  5.892151  6.313096  6.690013  7.137684  7.609802
    ##             70%       80%       90%     100%
    ##  [1,] 32.472000 39.548000 43.254000 62.06000
    ##  [2,] 28.557000 30.598000 32.705593 48.08000
    ##  [3,] 18.278795 20.932386 25.378146 33.23000
    ##  [4,] 15.736858 20.530000 23.290000 31.90000
    ##  [5,] 10.454543 13.858764 17.315281 20.81474
    ##  [6,]  9.262631 11.576707 14.261686 17.84559
    ##  [7,]  9.937324 13.482000 17.266000 24.79000
    ##  [8,]  8.974000  9.734000 11.792800 20.35791
    ##  [9,]  8.290000  8.817964 10.933019 18.15000
    ## [10,]  8.484000 10.018000 12.428890 23.37000
    ## [11,]  8.587000  9.593633 10.783000 16.59000
    ## [12,]  8.760000  9.300400 10.444000 19.83000
    ## [13,]  8.115099  8.787356  9.986058 13.30600

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
output <- SpeedClimbingABM(n = n, years = years, pop_data = pop_data,
                           n_holds = 20, beta_true_prob = 1, learn_prob = 0.2, n_top = 20,
                           innov_prob = 0.2, adj_poss = 2, improve_rate_m = 2, improve_rate_sd = 0.2, improve_min = 0.4,
                           learn_x_times = 1, innov_x_times = -1, learn_x_pop = 1, innov_x_pop = 1,
                           sum_stats = TRUE, plot = TRUE, bw = 0.6, ylim = 0.4)
```

![](README_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

``` r
#print run time
Sys.time() - start
```

    ## Time difference of 0.3436639 secs

``` r
output
```

    ##              0%       10%       20%       30%       40%       50%       60%
    ##  [1,] 15.400000 18.256000 20.408000 21.188000 24.432000 28.640000 29.712000
    ##  [2,]  8.080000 10.742886 12.480000 13.893208 15.533201 20.965456 25.160000
    ##  [3,]  6.386811  8.583561  9.841354 10.870492 11.375709 12.839691 15.285648
    ##  [4,]  5.577313  7.790000  8.676870  9.381616  9.880000 11.390000 12.690000
    ##  [5,]  5.218444  6.978805  7.360000  7.780580  8.500000  8.729471  9.052693
    ##  [6,]  5.059349  6.383908  7.120000  7.286414  7.933110  8.317257  8.829000
    ##  [7,]  5.563659  6.449403  6.882000  7.226000  7.636000  8.029496  8.553660
    ##  [8,]  5.188753  6.154628  6.650903  7.182901  7.708200  7.876217  8.554000
    ##  [9,]  5.161489  6.051033  6.526585  6.909945  7.233582  7.820720  8.032846
    ## [10,]  5.111258  6.043986  6.478000  6.862000  7.238015  7.646579  8.028000
    ## [11,]  3.934889  5.661430  6.249632  6.681000  7.135740  7.553466  8.094000
    ## [12,]  3.407496  5.108179  5.920000  6.526444  6.938784  7.504000  7.940341
    ## [13,]  3.404797  4.754109  5.389326  5.938797  6.478037  6.966250  7.459717
    ##             70%       80%       90%     100%
    ##  [1,] 32.472000 39.548000 43.254000 62.06000
    ##  [2,] 28.557000 30.598000 33.133000 48.08000
    ##  [3,] 18.686123 22.914965 25.526000 33.23000
    ##  [4,] 16.286312 22.068577 24.230000 31.90000
    ##  [5,] 10.605973 15.150350 18.937788 22.67335
    ##  [6,] 10.121246 11.601399 16.505278 21.39943
    ##  [7,] 10.753404 13.610000 17.266000 24.79000
    ##  [8,]  8.976150  9.734000 11.862000 22.71448
    ##  [9,]  8.502000  8.926733 11.262000 18.15000
    ## [10,]  8.654000 10.150000 12.355230 23.37000
    ## [11,]  8.854000  9.594059 11.055349 16.59000
    ## [12,]  8.668000  9.248506 10.322000 19.83000
    ## [13,]  7.981297  8.439097  9.527263 13.30600
