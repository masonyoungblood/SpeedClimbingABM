
#this script converts the raw data on climbing times from ifsc events, from
#https://github.com/lfigil/ifsc-data collected by lfigil on february 1, 2024,
#into the rdata files used as input to the simulation

#we manually checked whether the data include all world cups and championships from 2007 to 2019
#only missing the 2007 world championships in aviles, which we will add manually: https://ifsc.results.info/#/event/481/

#load files, temporarily removing the manually added entry
files <- list.files("/Users/masonyoungblood/Documents/Work/Fall 2023/BayesFlow Exploring/climbing_times/source_lfigil_ifsc-data")
missing_file <- files[grep("missing", files)]
files <- files[-grep("missing", files)]

#combine into single data frame
data <- do.call(rbind, lapply(files, function(x){
  read.csv(paste0("/Users/masonyoungblood/Documents/Work/Fall 2023/BayesFlow Exploring/climbing_times/source_lfigil_ifsc-data/", x))
}))

#restructure so each row is a single time and remove NA values
data <- data.frame(athlete = gsub(" ", "_", tolower(rep(paste0(data$Name, " ", data$X), 2))), 
                   gender = tolower(substr(rep(data$gender, 2), 1, 1)), 
                   time = c(data$Qualification, data$Final), 
                   event = tolower(gsub("\\(|\\)", "", rep(data$event_name, 2))), 
                   year = rep(as.numeric(gsub(".*20", "20", data$date)), 2))

#load missing file and combine with main data
missing_file <- read.csv(paste0("/Users/masonyoungblood/Documents/Work/Fall 2023/BayesFlow Exploring/climbing_times/source_lfigil_ifsc-data/", missing_file))
missing_file <- data.frame(athlete = gsub(" ", "_", tolower(rep(paste0(missing_file$first, " ", missing_file$last), 2))), 
                           gender = tolower(rep(missing_file$gender, 2)), 
                           time = c(missing_file$qualification, missing_file$final), 
                           event = "aviles_esp_2007", 
                           year = 2007)
data <- rbind(missing_file, data)

#reformat times as numeric, which converts failures and false starts into NA, and remove NA values
data$time <- as.numeric(data$time)
data <- data[-which(is.na(data$time)), ]

#check for impossible times (lower than the world record in 2019)
#qinghai in 2009 has incorrect data (ranks logged as decimals in impossible times) so drop entire event from data
#one entry from tarnov in 2009 has an impossible time of 2 seconds, which will be dropped from the data
#data[which(data$time < 5.48), ]
data <- data[-which(data$event == "qinghai_chn_2009"), ]
data <- data[-which(data$time == 2), ]

#there are two alex johnsons, recode their names with their gender to disambuiguate them
#https://en.wikipedia.org/wiki/Alex_Johnson_(climber)
#https://www.reddit.com/r/climbing/comments/sl8rm1/the_other_alex_johnson_showing_socksgang_can_send/
#unique(data$athlete)[which(lengths(lapply(unique(data$athlete), function(x){unique(data$gender[which(data$athlete == x)])})) == 2)]
data$athlete[which(data$athlete == "alex_johnson" & data$gender == "m")] <- "alex_johnson_m"
data$athlete[which(data$athlete == "alex_johnson" & data$gender == "w")] <- "alex_johnson_w"

#save all climbing times
all_climbing_times <- data
save(all_climbing_times, file = "/Users/masonyoungblood/Documents/Work/Fall 2023/BayesFlow Exploring/climbing_times/all_climbing_times.RData")

#convert to only include the best time from each climber in each year
best_climbing_times <- do.call(rbind, lapply(unique(data$year), function(x){
  #store only data from that year
  temp <- data[which(data$year == x), ]
  
  #iterate through the athletes from that year, compiling their best time in that year
  do.call(rbind, lapply(unique(temp$athlete), function(y){
    data.frame(athlete = y, gender = temp$gender[match(y, temp$athlete)], time = min(temp$time[which(temp$athlete == y)]), year = x)
  }))
}))

#save best climbing times
save(best_climbing_times, file = "/Users/masonyoungblood/Documents/Work/Fall 2023/BayesFlow Exploring/climbing_times/best_climbing_times.RData")
