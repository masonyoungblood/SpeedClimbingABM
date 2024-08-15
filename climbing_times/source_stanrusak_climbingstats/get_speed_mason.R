#load library and data
library(jsonlite)
data <- fromJSON("/Users/masonyoungblood/Documents/Work/Fall 2023/BayesFlow Exploring/climbing_times/source_stanrusak_climbingstats/full_data.json")

#set target years
years <- 2007:2019

# #confirm that naming is consistent across entire dataset
# unique(unlist(lapply(years, function(g){
#   #get events in that year
#   naming_data <- data[[which(names(data) == g)]]$events
#   
#   #return the categories
#   return(unique(unlist(lapply(1:length(naming_data), function(i){naming_data[[i]]$categories}))))
# })))

#get all climbing times between target years
#note that NA times either come from before times were recorded (early 2007 and before, obvious bc entire event)
#or denote that the climber fell (the rest of 2007 and beyond)
all_climbing_times <- do.call(rbind, lapply(years, function(g){
  #get index in data that corresponds to that year
  h <- which(names(data) == g)
  
  #store temporary version of data that only includes events from that year
  year_data <- data[[h]]$events
  
  #for each event in that year
  year_results <- lapply(1:length(year_data), function(i){
    #which events were speed climbing
    has_speed <- grep("SPEED", names(year_data[[i]]$results))
    
    #if it includes speed climbing
    if(length(has_speed) > 0){
      #then for each event that is speed climbing
      event_results <- do.call(rbind, lapply(has_speed, function(j){
        #get a data frame with the athlete ids and times
        gender_results <- do.call(rbind, lapply(1:nrow(year_data[[i]]$results[[j]]), function(x){
          #retrieve athlete id
          athlete_id <- as.numeric(year_data[[i]]$results[[j]]$athlete_id[x])
          
          #if the time is not empty, then return it
          if(length(as.numeric(year_data[[i]]$results[[j]]$rounds[[x]]$score)) == 0){
            time <- NA
          } else{
            time <- as.numeric(year_data[[i]]$results[[j]]$rounds[[x]]$score)
          }
          
          #combine into data frame
          data.frame(athlete_id = athlete_id, time = time)
        }))
        
        #if the athletes were men, add that as gender
        if(grepl(" Men", names(year_data[[i]]$results)[j])){gender_results$gender <- "M"}
        
        #if the athletes were women, add that as gender
        if(grepl(" Women", names(year_data[[i]]$results)[j])){gender_results$gender <- "W"}
        
        #add event
        gender_results$event <- year_data[[i]]$name
        
        #return data frame
        return(gender_results)
      }))
      
      #return data frame
      return(event_results)
    } else{
      #otherwise return NA
      return(NA)
    }
  })
  
  #combine all events into single object
  year_results <- do.call(rbind, year_results[-which(is.na(year_results))])
  
  #add year
  year_results$year <- g
  
  #return data frame
  return(year_results)
}))

#remove NA values (failure and before times were recorded)
all_climbing_times <- all_climbing_times[-which(is.na(all_climbing_times$time)), ]

#save data as r object
save(all_climbing_times, file = "/Users/masonyoungblood/Documents/Work/Fall 2023/BayesFlow Exploring/climbing_times/source_stanrusak_climbingstats/all_climbing_times.RData")
