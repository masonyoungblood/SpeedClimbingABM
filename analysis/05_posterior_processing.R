#set working directory, load data, source code
setwd("/Users/masonyoungblood/Documents/Work/Summer_2021/Speed Climbing/SpeedClimbingABM")
load("analysis/data_and_output/03_best_params/priors_and_simulations.RData")
load("climbing_times/best_climbing_times.RData"); data <- best_climbing_times; rm(best_climbing_times)
grid <- read.csv("grid.csv")/1000
source("SpeedClimbingABM.R")

#convert posteriors to usable table format
posterior_predictions <- jsonlite::read_json("analysis/data_and_output/04_best_params_inference/posterior_predictions.json")
posterior_predictions <- do.call(rbind, lapply(1:length(posterior_predictions), function(x){unlist(posterior_predictions[[x]])}))
posterior_predictions <- do.call(cbind, lapply(1:13, function(x){posterior_predictions[, x]*sd(priors_and_simulations[[1]][, x]) + mean(priors_and_simulations[[1]][, x])}))
posterior_predictions <- as.data.frame(posterior_predictions)
colnames(posterior_predictions) <- c("innov_prob", "learn_prob", "n_top", 
                                     "improve_min_m", "improve_min_w",
                                     "improve_rate_m", "improve_rate_w",
                                     "innov_x_year", "innov_x_times", "innov_x_pop",
                                     "learn_x_year", "learn_x_times", "learn_x_pop")

#prevent n_top from being < 0 or > 1
posterior_predictions$n_top[which(posterior_predictions$n_top < 0)] <- 0
posterior_predictions$n_top[which(posterior_predictions$n_top > 1)] <- 1

#get posteriors for athletic improvement for manuscript (not plotted)
quantile(posterior_predictions$improve_min_m, probs = c(0.025, 0.5, 0.975))
quantile(posterior_predictions$improve_min_w, probs = c(0.025, 0.5, 0.975))
quantile(posterior_predictions$improve_rate_m, probs = c(0.025, 0.5, 0.975))
quantile(posterior_predictions$improve_rate_w, probs = c(0.025, 0.5, 0.975))

# POSTERIOR SIMULATIONS ---------------------------------------------------

#get all unique climbers
uniq_climbers <- unique(data$athlete)

#separate continuous sequences of years
seqs <- lapply(uniq_climbers, function(x){split(data$year[which(data$athlete == x)], cumsum(seq_along(data$year[which(data$athlete == x)]) %in% (which(diff(data$year[which(data$athlete == x)]) > 1) + 1)))})

#for each unique climber, iterate through their sequences, and and extract their ID, start year, end year, and time in start year (separate row per sequence)
pop_data <- data.table::data.table(do.call(rbind, lapply(1:length(uniq_climbers), function(i){
  t(sapply(1:length(seqs[[i]]), function(j){
    c(uniq_climbers[i], min(unlist(seqs[[i]][j])), max(unlist(seqs[[i]][j])), data$time[which(data$athlete == uniq_climbers[i] & data$year == min(unlist(seqs[[i]][j])))])
  }))
})))
colnames(pop_data) <- c("ID", "start", "end", "time")
pop_data$time <- as.numeric(pop_data$time)
pop_data$gender <- toupper(data$gender[match(pop_data$ID, data$athlete)])
pop_data <- pop_data[, c(1, 5, 2, 3, 4)]

#get years
years <- sort(unique(data$year))

#get population sizes
n <- unlist(lapply(1:length(years), function(x){nrow(data[which(data$year == years[x]), ])}))

#number of simulations
n_sim <- 10000

#run simulations with everything varied
posterior_simulations <- parallel::mclapply(1:n_sim, function(x){
  SpeedClimbingABM(n = n, years = years, pop_data = pop_data, grid = grid, n_holds = 20,
                   beta_true_prob = 1,
                   innov_prob = posterior_predictions$innov_prob[x],
                   learn_prob = posterior_predictions$learn_prob[x],
                   n_top = posterior_predictions$n_top[x],
                   improve_min_m = posterior_predictions$improve_min_m[x],
                   improve_min_w = posterior_predictions$improve_min_w[x],
                   improve_rate_m = posterior_predictions$improve_rate_m[x],
                   improve_rate_w = posterior_predictions$improve_rate_w[x],
                   innov_x_year = posterior_predictions$innov_x_year[x],
                   innov_x_times = posterior_predictions$innov_x_times[x],
                   innov_x_pop = posterior_predictions$innov_x_pop[x],
                   learn_x_year = posterior_predictions$learn_x_year[x],
                   learn_x_times = posterior_predictions$learn_x_times[x],
                   learn_x_pop = posterior_predictions$learn_x_pop[x],
                   sum_stats = FALSE, plot = FALSE, raw = TRUE)
}, mc.cores = 7)

#save simulations
save(posterior_simulations, file = "analysis/data_and_output/05_posterior_processing/posterior_simulations.RData")

# POSTERIOR DISTRIBUTIONS -------------------------------------------------

#load libraries for plotting
library(ggplot2)
library(cowplot)

#create function for constructing plots
plot_constructor <- function(post, tol = 1000, xlim, bound = 0.001, two_axes = TRUE, xlabel = NULL, ylabel = NULL, xlabs = TRUE, ybuffer = 0.02, small = FALSE, font_size = 6, color){
  post <- density(post, from = xlim[1], to = xlim[2], bw = "SJ", kernel = "gaussian")
  plot_data <- data.frame(x = c(post$x), y = c(post$y))
  if(length(which(post$y < bound) > 0)){plot_data$y[which(plot_data$y < bound)] <- NA}
  temp <- ggplot(plot_data, aes(x, y)) + geom_line(aes(y = y), color = color) + 
    geom_segment(x = plot_data$x[which.max(plot_data$y)], xend = plot_data$x[which.max(plot_data$y)],
                 y = plot_data$y[which.max(plot_data$y)], yend = 0, color = color, linetype = "solid") + 
    xlim(xlim[1], xlim[2]) + xlab(xlabel) + ylab("Density") + theme_linedraw(base_size = font_size) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, max(plot_data$y, na.rm = TRUE) + max(plot_data$y, na.rm = TRUE)*ybuffer)) + 
    theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), text = element_text(family = "Avenir Next"))
  if(small){
    temp <- temp + theme(axis.text.y = element_blank()) + geom_vline(xintercept = 0, linetype = "dotted", size = 0.5) + ylab(ylabel)
  } 
  if(!xlabs){
    temp <- temp + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank())
  }
  return(temp)
}

#construct subpanels for plots
bounds <- c(-5, 5)
innov_prob <- plot_constructor(posterior_predictions$innov_prob, xlim = c(0, 1), xlabel = "Innovation Probability", color = "red")
innov_x_times <- plot_constructor(posterior_predictions$innov_x_times, xlim = bounds, small = TRUE, ylabel = "Current\nRecord", color = "red")
innov_x_pop <- plot_constructor(posterior_predictions$innov_x_pop, xlim = bounds, small = TRUE, ylabel = "Pop.\nSize", xlabel = "Effect on Innovation", color = "red")
innov_x_year <- plot_constructor(posterior_predictions$innov_x_year, xlim = bounds, small = TRUE, ylabel = "Time\nTrend", color = "red")
learn_prob <- plot_constructor(posterior_predictions$learn_prob, xlim = c(0, 1), xlabel = "Learning Probability", color = "blue")
learn_x_times <- plot_constructor(posterior_predictions$learn_x_times, xlim = bounds, small = TRUE, ylabel = "Current\nRecord", color = "blue")
learn_x_pop <- plot_constructor(posterior_predictions$learn_x_pop, xlim = bounds, small = TRUE, ylabel = "Pop.\nSize", xlabel = "Effect on Learning", color = "blue")
learn_x_year <- plot_constructor(posterior_predictions$learn_x_year, xlim = bounds, small = TRUE, ylabel = "Time\nTrend", color = "blue")
#improve_rate_m <- plot_constructor(posterior_predictions$improve_rate_m, xlim = c(1, 3), xlabel = "Rate of Athletic Improvement (Men)")
#improve_rate_w <- plot_constructor(posterior_predictions$improve_rate_w, xlim = c(1, 3), xlabel = "Rate of Athletic Improvement (Women)")

#set theme for cowplot
theme_set(theme_cowplot(font_family = "Avenir Next")) 

#combined and save panels
png("analysis/data_and_output/05_posterior_processing/posterior_distributions.png", units = "in", width = 7.2, height = 3, res = 1000)
space <- plot_grid(NULL)
innov_interactions <- plot_grid(innov_x_year, innov_x_times, innov_x_pop, ncol = 1, rel_heights = c(1, 1, 1.2))
first_row <- plot_grid(space, innov_prob, innov_interactions, nrow = 1, rel_widths = c(0.02, 1, 1))
learn_interactions <- plot_grid(learn_x_year, learn_x_times, learn_x_pop, ncol = 1, rel_heights = c(1, 1, 1.2))
second_row <- plot_grid(space, learn_prob, learn_interactions, nrow = 1, rel_widths = c(0.02, 1, 1))
#third_row <- plot_grid(space, improve_rate_m, improve_rate_w, nrow = 1, rel_widths = c(0.02, 1, 1))
plot_grid(first_row, space, second_row, ncol = 1, labels = c("A", "", "B"), rel_heights = c(0.6, 0.02, 0.6))
#plot_grid(first_row, space, second_row, space, third_row, ncol = 1, labels = c("A", "", "B", "", "C"), rel_heights = c(0.6, 0.02, 0.6, 0.02, 0.6))
dev.off()

# ROUTES AND HOLDS --------------------------------------------------------

#load libraries and data
library(ggnetwork)
library(network)
library(cowplot)
library(ggplot2)
#extrafont::loadfonts(quiet = TRUE)
load("analysis/data_and_output/05_posterior_processing/posterior_simulations.RData")
load("data.RData")
grid <- read.csv("grid.csv")/1000

#get proportion of attempted innovation and learning attempts that actually happened
actual_innovation <- mean(sapply(1:length(posterior_simulations), function(x){posterior_simulations[[x]]$act_innov/sum(sapply(2007:2019, function(y){length(unique(data$athlete[which(data$year == y)]))}))}), xlab = "", main = "")
actual_learning <- mean(sapply(1:length(posterior_simulations), function(x){posterior_simulations[[x]]$act_learn/sum(sapply(2007:2019, function(y){length(unique(data$athlete[which(data$year == y)]))}))}), xlab = "", main = "")

#compute diversity and evenness over time for each posterior simulation
diversity_data <- do.call(rbind, lapply(1:length(posterior_simulations), function(x){
  #compute hill numbers
  temp <- do.call(rbind, lapply(1:length(posterior_simulations[[x]]$n_unique_routes), function(y){
    c(hillR::hill_taxa(posterior_simulations[[x]]$n_unique_routes[[y]], q = 0),
      hillR::hill_taxa(posterior_simulations[[x]]$n_unique_routes[[y]], q = 1),
      hillR::hill_taxa(posterior_simulations[[x]]$n_unique_routes[[y]], q = 2))
  }))
  
  #evenness is computed as (diversity - 1)/(richness - 1), and cannot be computed for q = 0
  #as per 10.1126/science.abl7655, which uses E3 from 10.1002/ecy.2852
  temp <- cbind(x, 1:length(posterior_simulations[[x]]$n_unique_routes), temp, (temp[, 2] - 1)/(temp[, 1] - 1), (temp[, 3] - 1)/(temp[, 1] - 1))
  temp <- as.data.frame(temp)
  colnames(temp) <- c("sim", "t", "q0", "q1", "q2", "even_q1", "even_q2")
  return(temp)
}))

#create evenness plot
even_med <- data.frame(x = 1:max(diversity_data$t), y = sapply(1:max(diversity_data$t), function(x){median(diversity_data$even_q2[which(diversity_data$t == x)])}))
even_plot <- ggplot() + 
  geom_line(data = diversity_data, aes(x = t, y = even_q2, group = sim), color = scales::alpha("black", 0.01)) + 
  geom_line(data = even_med, aes(x = x, y = y), color = "red") + 
  ylab(expression(Evenness~(D^{"q=2"}~-~1)/(D^{"q=0"}~-~1))) + 
  scale_x_continuous(name = "Year", breaks = 1:max(diversity_data$t), labels = (2019-max(diversity_data$t)+1):2019, limits = c(1, max(diversity_data$t)), expand = c(0, 0)) + 
  theme_linedraw(base_size = 6) + 
  theme(text = element_text(family = "Avenir Next"), plot.margin = margin(5.5, 11, 5.5, 5.5))

#create diversity plot
div_med <- data.frame(x = 1:max(diversity_data$t), y = sapply(1:max(diversity_data$t), function(x){median(diversity_data$q2[which(diversity_data$t == x)])}))
div_plot <- ggplot() + 
  geom_line(data = diversity_data, aes(x = t, y = q2, group = sim), color = scales::alpha("black", 0.01)) + 
  geom_line(data = div_med, aes(x = x, y = y), color = "red") + 
  ylab(expression(Shannon~Diversity~(D^{"q=2"}))) + 
  scale_x_continuous(name = "Year", breaks = 1:max(diversity_data$t), labels = (2019-max(diversity_data$t)+1):2019, limits = c(1, max(diversity_data$t)), expand = c(0, 0)) + 
  theme_linedraw(base_size = 6) + 
  theme(text = element_text(family = "Avenir Next"), plot.margin = margin(5.5, 11, 5.5, 5.5))

#create richness plot
rich_med <- data.frame(x = 1:max(diversity_data$t), y = sapply(1:max(diversity_data$t), function(x){median(diversity_data$q0[which(diversity_data$t == x)])}))
rich_plot <- ggplot() + 
  geom_line(data = diversity_data, aes(x = t, y = q0, group = sim), color = scales::alpha("black", 0.01)) + 
  geom_line(data = rich_med, aes(x = x, y = y), color = "red") + 
  ylab(expression(Richness~(D^{"q=0"}))) + 
  scale_x_continuous(name = "Year", breaks = 1:max(diversity_data$t), labels = (2019-max(diversity_data$t)+1):2019, limits = c(1, max(diversity_data$t)), expand = c(0, 0)) + 
  theme_linedraw(base_size = 6) + 
  theme(text = element_text(family = "Avenir Next"), plot.margin = margin(5.5, 11, 5.5, 5.5))

#set theme for cowplot
theme_set(theme_cowplot(font_family = "Avenir Next")) 

#export richness and diversity plots
png("analysis/data_and_output/05_posterior_processing/diversity.png", units = "in", width = 7.2, height = 2.6, res = 1000)
plot_grid(rich_plot, div_plot, labels = c("A", "B"))
dev.off()

#simplify posterior simulations to remove actual innovation and learning
posterior_simulations <- lapply(1:length(posterior_simulations), function(x){posterior_simulations[[x]][[1]]})

#get all routes from all posterior simulations
all_routes <- parallel::mclapply(1:length(posterior_simulations), function(x){
  sapply(1:length(posterior_simulations[[x]]$beta), function(y){
    paste0(which(posterior_simulations[[x]]$beta[[y]]), collapse = " ")
  })
}, mc.cores = 7)
all_routes <- as.data.frame(sort(table(unlist(all_routes)), decreasing = TRUE))

#construct frequency table for histogram
freq_table <- all_routes
freq_table[, 1] <- sapply(1:nrow(freq_table), function(x){length(strsplit(as.character(freq_table$Var1[x]), " ")[[1]])})
freq_table <- aggregate(Freq ~ Var1, freq_table, FUN = sum)
freq_table[, 1] <- 20 - freq_table[, 1]
colnames(freq_table) <- c("skipped", "freq")

#compress them into an aggregated data frame
all_routes <- do.call(rbind, lapply(1:nrow(all_routes), function(x){
  temp <- c(as.numeric(strsplit(as.character(all_routes[x, 1]), " ")[[1]]), 21)
  data.frame(from = temp[1:(length(temp)-1)], to = temp[2:length(temp)], weight = all_routes[x, 2])
}))
all_routes <- aggregate(weight ~ from + to, all_routes, FUN = sum)

#get the frequency at which each hold is use and which hold is skipped
skipped <- sapply(1:21, function(x){
  temp <- all_routes$weight[which(all_routes$from == x - 1 & all_routes$to == x + 1)]
  if(length(temp) == 0){temp <- 0}
  return(temp)
})
#used <- unlist(parallel::mclapply(1:length(posterior_simulations), function(x){unlist(lapply(1:length(posterior_simulations[[x]]$beta), function(y){which(posterior_simulations[[x]]$beta[[y]])}))}, mc.cores = parallel::detectCores()-1))
#used <- as.numeric(table(used))
#used <- c(used, max(used))

#construct network for plotting and add attributes
net <- network(all_routes, directed = FALSE)
set.vertex.attribute(net, "shape", c(rep(19, 20), 17))
set.vertex.attribute(net, attrname = "size", value = skipped)
#set.vertex.attribute(net, attrname = "size", value = used/max(used))
net <- ggnetwork(net, layout = as.matrix(grid))

#create plot of transition densities
hold_plot <- ggplot(net, aes(x = x, y = y, xend = xend, yend = yend)) + 
  geom_edges(aes(linewidth = weight), curvature = -0.02) +
  geom_nodes(aes(size = size), color = "blue") + 
  scale_size(trans = "log", range = c(0, 4), guide = NULL) + 
  #scale_size(range = c(0, 4), guide = NULL) + 
  scale_linewidth(range = c(0, 1.2), guide = NULL) + 
  xlim(-0.05, 1.05) + 
  theme_void() + 
  scale_y_continuous(breaks = sort(unique(net$y))[-21], labels = c(1:20), expand = c(0, 0), limits = c(-0.05, 1.05))

#create empty plot to store the axis alone
hold_plot_axis <- ggplot(net, aes(x = x, y = y, xend = xend, yend = yend)) + 
  theme_linedraw(base_size = 6) + 
  theme(text = element_text(family = "Avenir Next")) + 
  xlim(-0.05, 1.05) + scale_y_continuous(breaks = sort(unique(net$y))[-21], labels = c(1:20), expand = c(0, 0), limits = c(-0.05, 1.05))

#create bar plot of number of holds skipped
bar_plot <- ggplot(freq_table, aes(fill = factor(1), x = skipped, y = freq)) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = c("limegreen"), name = "") + 
  geom_vline(xintercept = 5, linetype = "dashed", size = 0.5) + 
  theme_linedraw(base_size = 6) + 
  xlab("# Holds Skipped") + ylab("Count") + 
  scale_x_continuous(labels = scales::number_format(accuracy = 1), breaks = c(0:9)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA), breaks = c(0, 2e5, 4e5, 6e5), labels = c("0", "0.2M", "0.4M", "0.6M")) + 
  theme(strip.text.x = element_blank(), legend.position = "none", text = element_text(family = "Avenir Next"))

#construct combined object of priors, posteriors, and observed data
temp_prior <- do.call(rbind, parallel::mclapply(1:10000, function(x){
  temp <- density(priors_and_simulations[[2]][[x]][[12]], bw = "SJ", kernel = "gaussian")
  data.frame(x = temp$x, y = temp$y, a = "a", b = x)
}, mc.cores = 7))
temp_posterior_a <- do.call(rbind, parallel::mclapply(1:2000, function(x){
  temp <- density(posterior_simulations[[x]]$current_record, bw = "SJ", kernel = "gaussian")
  data.frame(x = temp$x, y = temp$y, a = "b", b = x)
}, mc.cores = 7))
temp_posterior_b <- temp_posterior_a
temp_posterior_b$a <- "c"
temp_observed <- data.frame(x = density(data$time[which(data$year == "2019")], bw = "SJ", kernel = "gaussian")$x,
                            y = density(data$time[which(data$year == "2019")], bw = "SJ", kernel = "gaussian")$y,
                            a = "d", b = 1)
curve_data <- rbind(temp_prior, temp_posterior_a, temp_posterior_b, temp_observed)

#create plot that includes priors, posteriors, and observed data
curve_plot <- ggplot(curve_data, aes(x = x, y = y, group = paste0(a, b), color = as.factor(a))) + 
  geom_line() + xlim(0, 25) +
  scale_linewidth_manual(values = c(5, 0.5, 0.5, 1)) + 
  scale_linetype_manual(values = c("solid", "solid", "solid", "dashed")) + 
  scale_color_manual(values = c(scales::alpha("gray88", 1), scales::alpha("white", 0.02), scales::alpha("red", 0.01), "black")) +
  theme_linedraw(base_size = 6) + theme(legend.position = "none", text = element_text(family = "Avenir Next")) + xlab("Climbing Times (s)") + ylab("Density") + 
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01))

#process things for cowplot
theme_set(theme_cowplot(font_family = "Avenir Next")) 
plot_a <- plot_grid(plot_grid(NULL), curve_plot, nrow = 1, rel_widths = c(0.06, 1))
plot_a_b <- plot_grid(plot_a, plot_grid(plot_grid(NULL), bar_plot, rel_widths = c(0.06, 1)), nrow = 2, rel_heights = c(1.5, 1), labels = c("A", "B"))
plot_c <- plot_grid(get_y_axis(hold_plot_axis), hold_plot, nrow = 1, rel_widths = c(0.2, 1, 0.2, 1))

#export
png("analysis/data_and_output/05_posterior_processing/times_holds_routes.png", units = "in", width = 3.5, height = 3, res = 1000)
plot_grid(plot_a_b, plot_c, nrow = 1, labels = c("", "C"), rel_widths = c(2, 0.8))
dev.off()

# CLIMBING TIMES ----------------------------------------------------------

#load libraries for plotting
library(ggplot2)
library(cowplot)

#make copy of the data for the plot of climbing times
times_plot_data <- data
times_plot_data$gender <- ifelse(times_plot_data$gender == "M", "Men", "Women")

#create plot of climbing times for men and women
times_plot <- ggplot(data = times_plot_data, aes(x = year, y = time, color = gender)) + geom_jitter(width = 0.25, height = 0, size = 0.2) + 
  geom_smooth(span = 1, se = FALSE, size = 0.5) + 
  xlab("Year") + ylab("Time (s)") + theme_linedraw(base_size = 6) + 
  scale_x_continuous(breaks = c(2007, 2010, 2013, 2016, 2019)) + scale_y_continuous(expand = c(0, 0), limits = c(0, 90)) + 
  scale_color_manual(values = c("blue", "red")) + theme(legend.title = element_blank(), 
                                                        legend.position = c(0.825, 0.845), 
                                                        legend.background = element_blank(),
                                                        legend.key=element_blank(),
                                                        panel.grid.minor.x = element_blank(),
                                                        text = element_text(family = "Avenir Next"))

#create speed wall plot with adjustable buffers
theme_set(theme_cowplot(font_family = "Avenir Next")) 
speed_wall <- plot_grid(ggplot(NULL) + theme_void(), 
                        plot_grid(ggplot(NULL) + theme_void(), 
                                  ggdraw() + draw_image("speed_wall.png"), 
                                  ggplot(NULL) + theme_void(), ncol = 1, rel_heights = c(0.1, 1, 0.07)), 
                        ggplot(NULL) + theme_void(), 
                        nrow = 1, rel_widths = c(0.05, 1, 0.05))

#export
png("analysis/data_and_output/05_posterior_processing/climbing_times.png", width = 5, height = 3, units = "in", res = 1000)
plot_grid(speed_wall, times_plot, labels = c("A", "B"), rel_widths = c(0.4, 1))
dev.off()
