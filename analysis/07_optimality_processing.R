#set working directory and load data and libraries for plotting
setwd("/Users/masonyoungblood/Documents/Work/Summer_2021/Speed Climbing/SpeedClimbingABM")
load("analysis/data_and_output/06_optimality_sims/optimality_sims.RData")
load("analysis/data_and_output/03_best_params/priors_and_simulations.RData")
posterior_predictions <- jsonlite::read_json("analysis/data_and_output/04_best_params_inference/posterior_predictions.json")
library(ggplot2)
library(cowplot)

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

#get coordinates for innovation and learning for optimality plots
coords <- c(median(posterior_predictions$innov_prob), median(posterior_predictions$learn_prob))

#aggregate values across combinations of innovation and learning probability
plot_data <- cbind(aggregate(min ~ innov_prob + learn_prob, optimality_sims, FUN = median),
                   aggregate(median ~ innov_prob + learn_prob, optimality_sims, FUN = median)$median,
                   aggregate(unique ~ innov_prob + learn_prob, optimality_sims, FUN = median)$unique)
colnames(plot_data) <- c("innov_prob", "learn_prob", "min", "median", "unique")

#create heatplot of minimum times
min_plot <- ggplot(plot_data, aes(x = innov_prob, y = learn_prob, fill = min)) + geom_tile() + xlab("Innovation Probability") + ylab("Copying Probability") + 
  scale_fill_gradientn(name = "Time (s)", colours = c("red", "white", "blue"), values = scales::rescale(c(min(plot_data$min), median(plot_data$min), max(plot_data$min)))) +
  ggtitle(bquote(bold("Minimum Time"))) + theme_linedraw(base_size = 6) + 
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) + 
  geom_point(aes(x = coords[1], y = coords[2]), size = 3) +
  theme(text = element_text(family = "Avenir Next"))

#create heatplot of median times
med_plot <- ggplot(plot_data, aes(x = innov_prob, y = learn_prob, fill = median)) + geom_tile() + xlab("Innovation Probability") + ylab("Copying Probability") + 
  scale_fill_gradientn(name = "Time (s)", colours = c("red", "white", "blue"), values = scales::rescale(c(min(plot_data$median), median(plot_data$median), max(plot_data$median)))) +
  ggtitle(bquote(bold("Median Time"))) + theme_linedraw(base_size = 6) + 
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) + 
  geom_point(aes(x = coords[1], y = coords[2]), size = 3) +
  theme(text = element_text(family = "Avenir Next"))

#create heatplot of unique number of routes
unique_plot <- ggplot(plot_data, aes(x = innov_prob, y = learn_prob, fill = unique)) + geom_tile() + xlab("Innovation Probability") + ylab("Copying Probability") + 
  scale_fill_gradientn(name = "# Routes", colours = c("red", "white", "blue"), values = scales::rescale(c(min(plot_data$unique), median(plot_data$unique), max(plot_data$unique)))) +
  ggtitle(bquote(bold("Unique Routes"))) + theme_linedraw(base_size = 6) + 
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) + 
  geom_point(aes(x = coords[1], y = coords[2]), size = 3) +
  theme(text = element_text(family = "Avenir Next"))

#export
png("analysis/data_and_output/07_optimality_processing/optimality_heatmaps.png", units = "in", width = 7.2, height = 1.9, res = 1000)
plot_grid(min_plot + theme(legend.position = "none"),
          get_legend(min_plot),
          med_plot + theme(legend.position = "none"),
          get_legend(med_plot),
          unique_plot + theme(legend.position = "none"),
          get_legend(unique_plot),
          nrow = 1, rel_widths = c(1, 0.3, 1, 0.3, 1, 0.3))
dev.off()
