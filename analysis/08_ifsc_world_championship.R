
#this script statistically models and then plots the manually-transcribed 
#climbing sequences from ifsc events

# MODELING ----------------------------------------------------------------

#load libraries
library(brms)
library(mice)
library(miceadds)
library(parallel)
library(posterior)
library(car)
library(lme4)
library(showtext)

#set font family
font = "Myriad Pro Condensed"
font_add(font, regular = "/Users/masonyoungblood/Library/Fonts/MyriadPro-Cond.otf")
showtext_auto()

#set working directory
setwd("/Users/masonyoungblood/Documents/Work/Summer_2021/Speed Climbing/SpeedClimbingABM")

#create function to extract diagnostics from brm_multiple
diagnostics <- function(model, num_imputed){
  #get draws for each imputed dataset
  draws <- as_draws_array(model)
  draws_per_dat <- lapply(1:num_imputed, \(i) subset_draws(draws, chain = i))
  
  #get diagnostics for each dataset
  lapply(draws_per_dat, summarise_draws, default_convergence_measures())
}

#based on vignette
#https://cran.r-project.org/web/packages/brms/vignettes/brms_missings.html

#load in data
data <- rbind(
  cbind(year = 2012, read.csv("climbing_times/ifsc_world_championship/fra_2012.csv")),
  cbind(year = 2014, read.csv("climbing_times/ifsc_world_championship/spa_2014.csv")),
  cbind(year = 2016, read.csv("climbing_times/ifsc_world_championship/fra_2016.csv")),
  cbind(year = 2018, read.csv("climbing_times/ifsc_world_championship/aus_2018.csv")),
  cbind(year = 2019, read.csv("climbing_times/ifsc_world_championship/jpn_2019.csv"))
)
ranks_2016 <- read.csv("climbing_times/ifsc_world_championship/fra_2016_ranks.csv")

#remove athletes with NA values in the route (e.g. from falls), or have text in the time
data <- data[-which(is.na(rowMeans(data[, grep("hold", colnames(data))]))), ]
data <- data[-which(is.na(as.numeric(data$first_final_s))), ]

#remove sources
data <- data[, -which(colnames(data) %in% c("source_a", "source_b", "source_c"))]

#restructure for analysis
data$first_final_s <- as.numeric(data$first_final_s)
data$gender <- factor(data$gender)
data$name <- factor(data$name)
data$birthdate <- as.Date(data$birthdate, "%d %b %y")
data$age_bio <- as.numeric(as.Date("2019-08-17")-data$birthdate)
data$age_comp <- as.numeric(as.Date("2019-08-17")-as.Date(as.character(data$active_since), "%Y"))
data <- data[, -which(colnames(data) %in% c("active_since", "birthdate"))]

#get proportion of complete data for height and weight
length(which(!is.na(data$age_bio)))/nrow(data)
length(which(!is.na(data$age_comp)))/nrow(data)
length(which(!is.na(data$height_cm)))/nrow(data)
length(which(!is.na(data$weight_kg)))/nrow(data)

#set parameters for model fitting
num_imputed <- 10
chains <- 8
iters <- 5000

#generate imputed datasets
imp_data <- mice(data, m = num_imputed, method = "rf")
reza_imp_data <- datlist2mids(
  #iterate through imputed datasets
  lapply(lapply(1:num_imputed, function(x){complete(imp_data, x)}), function(data){
    #get people who are present both before and after the reza, excluding the man himself
    post_reza_ppl <- unique(data$name[which(data$year %in% c(2018, 2019))])
    pre_reza_ppl <- unique(data$name[which(data$year %in% c(2016))])
    target_ppl <- post_reza_ppl[which(post_reza_ppl %in% pre_reza_ppl)]
    target_ppl <- target_ppl[-which(target_ppl == "Reza Alipour")]
    
    #return processed data frame with variables
    data.frame(name = target_ppl,
               hold_4 = sapply(target_ppl, function(x){min(data$hold_4[which(data$name == x & data$year %in% c(2018, 2019))])}),
               pre_rank = ranks_2016$rank[match(target_ppl, ranks_2016$name)],
               gender = data$gender[match(target_ppl, data$name)],
               age_bio = data$age_bio[match(target_ppl, data$name)],
               age_comp = data$age_comp[match(target_ppl, data$name)],
               height_cm = data$height[match(target_ppl, data$name)],
               weight_kg = data$weight[match(target_ppl, data$name)])
  })
)

#vif exploration
#vif(lmer(scale(first_final_s) ~ scale(age_bio) + scale(age_comp) + scale(height_cm) + scale(weight_kg) + (1|gender) + (1|name) + (1|year), data = complete(imp_data)))

#set priors for model of time
prior <- get_prior(scale(first_final_s) ~ scale(age_bio) + scale(age_comp) + scale(height_cm) + scale(weight_kg) + (1|gender) + (1|name) + (1|year), 
                   data = data, family = gaussian())
prior$prior[which(prior$class == "b")] <- "normal(0, 2)"

#run model of time
time_model <- brm_multiple(scale(first_final_s) ~ scale(age_bio) + scale(age_comp) + scale(height_cm) + scale(weight_kg) + (1|gender) + (1|name) + (1|year), 
                           data = imp_data, family = gaussian(), 
                           prior = prior, iter = iters, chains = chains, cores = chains)

#check which holds have variation across climbers
apply(data[, grep("hold", colnames(data))], 2, sd)

#set priors for model of hold use
#initially planned on cauchy(0, 1), tighter than the gelman (2008) recommendation of (0, 2.5) because we scale to 1 rather than 0.5
#https://projecteuclid.org/journals/annals-of-applied-statistics/volume-2/issue-4/A-weakly-informative-default-prior-distribution-for-logistic-and-other/10.1214/08-AOAS191.full?tab=ArticleLink
#now using normal(0, 1), mentioned by gelman in these forums
#https://statmodeling.stat.columbia.edu/2015/11/01/cauchy-priors-for-logistic-regression-coefficients/
#https://discourse.mc-stan.org/t/how-to-scale-your-priors-by-number-of-predictors-in-logistic-regression/33263/19
#also note that gender was included as a random effect, because it doesn't predict any holds when fixed effect
prior <- get_prior(mvbind(hold_2, hold_4, hold_5, hold_9, hold_10, hold_14, hold_15, hold_16) ~ scale(age_bio) + scale(age_comp) + scale(height_cm) + scale(weight_kg) + gender + (1|name) + (1|year), 
                   data = data, family = bernoulli())
prior$prior[which(prior$class == "b")] <- "normal(0, 1)"

#run model of hold use
hold_model <- brm_multiple(mvbind(hold_2, hold_4, hold_5, hold_9, hold_10, hold_14, hold_15, hold_16) ~ scale(age_bio) + scale(age_comp) + scale(height_cm) + scale(weight_kg) + gender + (1|name) + (1|year),
                           data = imp_data, family = bernoulli(), 
                           prior = prior, iter = iters, chains = chains, cores = chains)

#set priors for reza model
prior <- get_prior(hold_4 ~ scale(age_bio) + scale(age_comp) + scale(height_cm) + scale(weight_kg) + gender + scale(pre_rank), 
                   data = complete(reza_imp_data, 1), family = bernoulli())
prior$prior[which(prior$class == "b")] <- "normal(0, 1)"

#run reza model
reza_model <- brm_multiple(hold_4 ~ scale(age_bio) + scale(age_comp) + scale(height_cm) + scale(weight_kg) + gender + scale(pre_rank), 
                           data = reza_imp_data, family = bernoulli(), 
                           prior = prior, iter = iters, chains = chains, cores = chains)

#get diagnostics for both models
time_diagnostics <- diagnostics(time_model, num_imputed)
hold_diagnostics <- diagnostics(hold_model, num_imputed)
reza_diagnostics <- diagnostics(reza_model, num_imputed)

#create table of time results
time_table <- round(summary(time_model)$fixed[-1, c(1:4)], 3)
rownames(time_table) <- c("Age", "Career", "Height", "Weight")
time_table$Sig <- ""
time_table$Sig[which(time_table$`l-95% CI`*time_table$`u-95% CI` > 0)] <- "*"
time_table

#create table of hold results
hold_table <- round(summary(hold_model)$fixed[-c(1:8), c(1:4)], 3)
rownames(hold_table) <- paste0(rep(c("Hold 2", "Hold 4", "Hold 5", "Hold 9", "Hold 10", "Hold 14", "Hold 15", "Hold 16"), each = 5), ": ", rep(c("Age", "Career", "Height", "Weight", "Gender (W)"), 8))
hold_table$Sig <- ""
hold_table$Sig[which(hold_table$`l-95% CI`*hold_table$`u-95% CI` > 0)] <- "*"
hold_table

#create table of reza results
reza_table <- round(summary(reza_model)$fixed[-1, c(1:4)], 3)
rownames(reza_table) <- c("Age", "Career", "Height", "Weight", "Gender", "Past Rank")
reza_table$Sig <- ""
reza_table$Sig[which(reza_table$`l-95% CI`*reza_table$`u-95% CI` > 0)] <- "*"
reza_table

#create table of time diagnostics
time_diag_table <- data.frame(do.call(rbind, lapply(time_diagnostics, function(x){x$rhat})))
colnames(time_diag_table) <- time_diagnostics[[1]]$variable
time_diag_table <- time_diag_table[, grep("b_", colnames(time_diag_table))]
time_diag_table <- time_diag_table[, -grep("_Intercept", colnames(time_diag_table))]
colnames(time_diag_table) <- c("Age", "Career", "Height", "Weight")
time_diag_table <- round(time_diag_table, 3)

#create table of hold diagnostics
hold_diag_table <- data.frame(do.call(rbind, lapply(hold_diagnostics, function(x){x$rhat})))
colnames(hold_diag_table) <- hold_diagnostics[[1]]$variable
hold_diag_table <- hold_diag_table[, grep("b_", colnames(hold_diag_table))]
hold_diag_table <- hold_diag_table[, -grep("_Intercept", colnames(hold_diag_table))]
colnames(hold_diag_table) <- paste0(rep(c("Hold 2", "Hold 4", "Hold 5", "Hold 9", "Hold 10", "Hold 14", "Hold 15", "Hold 16"), each = 5), ": ", rep(c("Age", "Career", "Height", "Weight", "Gender (W)"), 8))
hold_diag_table <- round(hold_diag_table, 3)

#create table of reza diagnostics
reza_diag_table <- data.frame(do.call(rbind, lapply(reza_diagnostics, function(x){x$rhat})))
colnames(reza_diag_table) <- reza_diagnostics[[1]]$variable
reza_diag_table <- reza_diag_table[, grep("b_", colnames(reza_diag_table))]
reza_diag_table <- reza_diag_table[, -grep("_Intercept", colnames(reza_diag_table))]
colnames(reza_diag_table) <- c("Age", "Career", "Height", "Weight", "Gender", "Past Rank")
reza_diag_table <- round(reza_diag_table, 3)

#save models
models <- list(data = imp_data, 
               time = list(model = time_table, diagnostics = time_diag_table), 
               hold = list(model = hold_table, diagnostics = hold_diag_table), 
               reza = list(model = reza_table, diagnostics = reza_diag_table))
save(models, file = "analysis/data_and_output/08_ifsc_world_championship/ifsc_models.RData")

# UNSEEN SPECIES ----------------------------------------------------------

#set working directory
setwd("/Users/masonyoungblood/Documents/Speed Climbing/SpeedClimbingABM")

#load libraries
library(iNEXT)
library(ggplot2)
library(cowplot)

#load in data
data <- rbind(
  cbind(year = 2012, read.csv("climbing_times/ifsc_world_championship/fra_2012.csv")),
  cbind(year = 2014, read.csv("climbing_times/ifsc_world_championship/spa_2014.csv")),
  cbind(year = 2016, read.csv("climbing_times/ifsc_world_championship/fra_2016.csv")),
  cbind(year = 2018, read.csv("climbing_times/ifsc_world_championship/aus_2018.csv")),
  cbind(year = 2019, read.csv("climbing_times/ifsc_world_championship/jpn_2019.csv"))
)

#remove athletes with NA values in the route (e.g. from falls), or have text in the time
data <- data[-which(is.na(rowMeans(data[, grep("hold", colnames(data))]))), ]
data <- data[-which(is.na(as.numeric(data$first_final_s))), ]

#remove sources
data <- data[, -which(colnames(data) %in% c("source_a", "source_b", "source_c"))]

#restructure for analysis
data$first_final_s <- as.numeric(data$first_final_s)
data$gender <- factor(data$gender)
data$name <- factor(data$name)
data$birthdate <- as.Date(data$birthdate, "%d %b %y")
data$age_bio <- as.numeric(as.Date("2019-08-17")-data$birthdate)
data$age_comp <- as.numeric(as.Date("2019-08-17")-as.Date(as.character(data$active_since), "%Y"))
data <- data[, -which(colnames(data) %in% c("active_since", "birthdate"))]

#get estimate for asymptotic species richness
holds_2019 <- data[which(data$year == 2019), which(grepl("hold", colnames(data)))]
route_freqs <- as.numeric(sort(table(sapply(1:nrow(holds_2019), function(x){paste0(holds_2019[x, ], collapse = "")})), decreasing = TRUE))
iNEXT(route_freqs, nboot = 1000)$AsyEst[1, ]

# PLOTTING ----------------------------------------------------------------

#load libraries
library(tidyverse)
library(ggplot2)
library(patchwork)
library(cowplot)
library(showtext)

#set font family
font = "Myriad Pro Condensed"
font_add(font, regular = "/Users/masonyoungblood/Library/Fonts/MyriadPro-Cond.otf")
showtext_auto()

#set working directory
setwd("/Users/masonyoungblood/Documents/Work/Summer_2021/Speed Climbing/SpeedClimbingABM")

#load in data
data <- rbind(
  cbind(year = 2012, read.csv("climbing_times/ifsc_world_championship/fra_2012.csv")),
  cbind(year = 2014, read.csv("climbing_times/ifsc_world_championship/spa_2014.csv")),
  cbind(year = 2016, read.csv("climbing_times/ifsc_world_championship/fra_2016.csv")),
  cbind(year = 2018, read.csv("climbing_times/ifsc_world_championship/aus_2018.csv")),
  cbind(year = 2019, read.csv("climbing_times/ifsc_world_championship/jpn_2019.csv"))
)

#remove athletes with NA values in the route (e.g. from falls), or have text in the time
data <- data[-which(is.na(rowMeans(data[, grep("hold", colnames(data))]))), ]
data <- data[-which(is.na(as.numeric(data$first_final_s))), ]

#remove sources
data <- data[, -which(colnames(data) %in% c("source_a", "source_b", "source_c"))]

#restructure for analysis
data$first_final_s <- as.numeric(data$first_final_s)
data$gender <- factor(data$gender)
data$name <- factor(data$name)
data$birthdate <- as.Date(data$birthdate, "%d %b %y")
data$age_bio <- as.numeric(as.Date("2019-08-17")-data$birthdate)
data$age_comp <- as.numeric(as.Date("2019-08-17")-as.Date(as.character(data$active_since), "%Y"))
data <- data[, -which(colnames(data) %in% c("active_since", "birthdate"))]

#find athletes who appear in at least 4 years
frequent_athletes <- data %>%
  group_by(name) %>%
  summarize(years_competed = n_distinct(year)) %>%
  filter(years_competed >= 4) %>%
  pull(name)

#create a complete dataset with all years for each athlete
all_years <- sort(unique(data$year))
all_holds <- 1:20

#create a complete grid of all athlete-year-hold combinations
complete_grid <- expand.grid(
  name = frequent_athletes,
  year = all_years,
  hold = all_holds,
  stringsAsFactors = FALSE
)

#join with actual data
evolution_matrix <- complete_grid %>%
  left_join(
    data %>%
      filter(name %in% frequent_athletes) %>%
      select(name, year, gender, first_final_s, starts_with("hold_")) %>%
      pivot_longer(cols=starts_with("hold_"), 
                   names_to="hold_col", 
                   values_to="used") %>%
      mutate(
        hold = as.numeric(str_replace(hold_col, "hold_", "")),
        name = as.character(name)
      ) %>%
      select(-hold_col),
    by = c("name", "year", "hold")
  ) %>%
  #fill in gender for missing years (it's constant per athlete)
  group_by(name) %>%
  fill(gender, .direction = "downup") %>%
  ungroup() %>%
  #create a data status flag
  mutate(
    data_status = case_when(
      is.na(first_final_s) ~ "missing",
      !is.na(used) ~ "present",
      TRUE ~ "error"
    )
  ) %>%
  #extract last name for sorting
  mutate(
    #extract last name by taking the last word in the name
    last_name = sapply(name, function(x) {
      parts <- strsplit(x, " ")[[1]]
      return(parts[length(parts)])
    })
  ) %>%
  #order by last name
  arrange(last_name, year, hold)

#create a factor of names ordered by last name
ordered_names <- evolution_matrix %>%
  select(name, last_name) %>%
  distinct() %>%
  arrange(last_name) %>%
  pull(name)

#create a separate plot for each athlete
plot_list <- list()
for(athlete in ordered_names) {
  athlete_data <- evolution_matrix %>% filter(name == athlete)
  athlete_gender <- unique(athlete_data$gender)[1]
  
  p <- ggplot(athlete_data, aes(y = hold, x = factor(year))) +
    #add grey background for missing years
    geom_tile(data = subset(athlete_data, data_status == "missing"),
              fill = "grey90", color = NA) +
    #add colored tiles for present data
    geom_tile(data = subset(athlete_data, data_status == "present"),
              aes(fill = factor(used)), color = NA) +
    #custom fill colors based on gender
    scale_fill_manual(values = if (athlete_gender == "m") {
      c("0" = "white", "1" = "blue")
    } else {
      c("0" = "white", "1" = "red")
    }) +
    #force x-axis to show all years
    scale_x_discrete(limits = as.character(all_years), labels = c("'12", "'14", "'16", "'18", "'19")) +
    #y-axis with all numbers but only labels at 1, 5, 10, 15, 20
    scale_y_continuous(
      breaks = 1:20,  #ticks for all values
      labels = c(rep("", 4), "5", rep("", 4), "10", rep("", 4), "15", rep("", 4), "20"),  #labels only at key points
      limits = c(0.5, 20.5),  #adjusted to center the tiles on integer positions
      expand = c(0, 0)  #remove padding
    ) +
    #add athlete name as title
    ggtitle(athlete) +
    #no axis labels
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 6) +
    theme(
      text = element_text(family = font),
      plot.title = element_text(size = 6),
      panel.spacing = unit(1, "lines"),
      axis.ticks.y = element_line(color = "black"),
      axis.ticks.x = element_line(color = "black"),
      axis.text.y = element_text(color = "black"),
      axis.text.x = element_text(color = "black"),
      legend.position = "none",
      panel.grid = element_blank()
    ) +
    #ensure the tiles are centered on integer positions
    coord_cartesian(ylim = c(0.5, 20.5))
  
  plot_list[[athlete]] <- p
}

#create a legend plot for the bottom right slot
legend_plot <- ggplot(data.frame(gender = c("Men", "Women"), y = c(1, 2)), aes(x = 1, y = y)) +
  geom_text(aes(label = gender), color = c("blue", "red"), family = font, fontface = "bold", size = 2) +
  scale_y_continuous(limits = c(-1, 4)) + 
  theme_void(base_size = 6)

#add the legend to the plot list
plot_list[["legend"]] <- legend_plot

#manually choose who gets shown
plot_list <- list(plot_list$`Reza Alipour`, 
                  #plot_list$`Danyil Boldyrev`,
                  #plot_list$`Patrycja Chudziak`,
                  #plot_list$`Marcin Dzieński`,
                  plot_list$`Anouck Jaubert`, 
                  plot_list$`Iuliia Kaplina`, 
                  #plot_list$`Stanislav Kokorin`
                  plot_list$`Mariia Krasavina`, 
                  #plot_list$`Aleksandra Mirosław`
                  plot_list$`Zhong Qixin`,
                  #plot_list$`Dmitrii Timofeev`,
                  plot_list$legend)

#combine all plots into a grid with 3 columns
athletes_plot <- wrap_plots(plot_list, ncol = 3)

#display the combined plot
athletes_plot

#prepare data for heatmap visualization
heatmap_data <- data %>%
  #identify which holds are skipped in each attempt
  mutate(
    hold_1_skipped = ifelse(hold_1 == 0, 1, 0),
    hold_2_skipped = ifelse(hold_2 == 0, 1, 0),
    hold_3_skipped = ifelse(hold_3 == 0, 1, 0),
    hold_4_skipped = ifelse(hold_4 == 0, 1, 0),
    hold_5_skipped = ifelse(hold_5 == 0, 1, 0),
    hold_6_skipped = ifelse(hold_6 == 0, 1, 0),
    hold_7_skipped = ifelse(hold_7 == 0, 1, 0),
    hold_8_skipped = ifelse(hold_8 == 0, 1, 0),
    hold_9_skipped = ifelse(hold_9 == 0, 1, 0),
    hold_10_skipped = ifelse(hold_10 == 0, 1, 0),
    hold_11_skipped = ifelse(hold_11 == 0, 1, 0),
    hold_12_skipped = ifelse(hold_12 == 0, 1, 0),
    hold_13_skipped = ifelse(hold_13 == 0, 1, 0),
    hold_14_skipped = ifelse(hold_14 == 0, 1, 0),
    hold_15_skipped = ifelse(hold_15 == 0, 1, 0),
    hold_16_skipped = ifelse(hold_16 == 0, 1, 0),
    hold_17_skipped = ifelse(hold_17 == 0, 1, 0),
    hold_18_skipped = ifelse(hold_18 == 0, 1, 0),
    hold_19_skipped = ifelse(hold_19 == 0, 1, 0),
    hold_20_skipped = ifelse(hold_20 == 0, 1, 0)
  )

#calculate the proportion of climbers skipping each hold by year
hold_skip_proportions <- heatmap_data %>%
  group_by(year) %>%
  summarize(
    across(ends_with("_skipped"), mean, na.rm = TRUE)
  ) %>%
  pivot_longer(
    cols = ends_with("_skipped"),
    names_to = "hold",
    values_to = "skip_proportion"
  ) %>%
  mutate(
    hold_number = as.numeric(str_extract(hold, "\\d+"))
  )

#create the heatmap with your original color scheme and styling matching individual athlete plots
heatmap_plot <- ggplot(hold_skip_proportions, aes(y = hold_number, x = factor(year, levels = sort(unique(year), decreasing = FALSE)), fill = skip_proportion)) +
  geom_tile(color = "white", size = 0.1) +
  geom_text(aes(label = ifelse(skip_proportion > 0, scales::percent(skip_proportion, accuracy = 1), "")), 
            color = ifelse(hold_skip_proportions$skip_proportion > 0.5, "black", "white"), size = 1.25) +
  scale_fill_gradient(low = "black", high = "white", name = "Skipped", labels = scales::percent_format()) +
  scale_y_continuous(breaks = 1:20, labels = 1:20, limits = c(0.5, 20.5), expand = c(0, 0)) +
  labs(y = "Hold", x = "") +
  ggtitle("") + 
  theme_minimal(base_size = 6) +
  theme(
    text = element_text(family = font),
    axis.ticks.y = element_line(color = "black"),
    axis.ticks.x = element_line(color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.text.x = element_text(color = "black"),
    legend.position = "right",
    panel.grid = element_blank()
  ) +
  coord_cartesian(ylim = c(0.5, 20.5))

#set theme for cowplot
theme_set(theme_cowplot(font_family = font)) 

#combined and save panels
png("analysis/data_and_output/08_ifsc_world_championship/route_maps.png", units = "in", width = 4.5, height = 2, res = 1000)
cowplot::plot_grid(athletes_plot, ggplot(NULL) + theme_void(), heatmap_plot, labels = c("A", "B", ""), rel_widths = c(1, 0.04, 0.7), nrow = 1)
dev.off()

#also save as eps
ggsave(
  "analysis/data_and_output/08_ifsc_world_championship/route_maps.eps",
  cowplot::plot_grid(athletes_plot, ggplot(NULL) + theme_void(), heatmap_plot, labels = c("A", "B", ""), rel_widths = c(1, 0.04, 0.7), nrow = 1),
  width = 4.5, height = 2
)
