library(dplyr)
library(tidyr)
library(lubridate)

# calls for time series df, identified start minute, and identified end minute: see example below
extract_intervals <- function(df, start_minute, end_minute, timezone = "America/New_York") {
  
  # Extract start times
  start <- df %>% 
    filter(minute(Time) == start_minute, second(Time) == 0) %>% 
    select(Time) %>% 
    mutate(date = date(Time),
           hour = hour(Time),
           minute = minute(Time),
           sec = second(Time))
  
  # Extract end times
  end <- df %>% 
    filter(minute(Time) == end_minute, second(Time) == 59) %>% 
    select(Time) %>% 
    mutate(date = date(Time),
           hour = hour(Time),
           minute = minute(Time),
           sec = second(Time))
  
  # Combine start and end
  startend <- rbind(start, end)
  
  # Remove duplicate entries
  startend <- startend %>% 
    select(Time, minute) %>% 
    unique()
  
  # Add row numbers
  startend$row <- seq_len(nrow(startend))
  
  # Pivot wider
  startend <- pivot_wider(data = startend,
                          names_from = minute,
                          values_from = Time)
  
  # Remove row column
  startend$row <- NULL
  
  # Extract and clean start and end columns
  startend_start <- na.omit(startend[[as.character(start_minute)]])
  startend_end <- na.omit(startend[[as.character(end_minute)]])
  
  # Combine start and end columns
  startend <- cbind(startend_start, startend_end)
  
  # Convert to data frame
  startend <- as_data_frame(startend) 
  
  # Rename columns
  colnames(startend) <- c("Start", "End")
  
  # Convert to datetime with specified timezone
  startend <- startend %>% 
    mutate(Start = as_datetime(Start, tz = timezone),
           End = as_datetime(End, tz = timezone))
  
  return(startend)
}

# Example usage:
# Replace 'your_data_frame' with the actual name of your data frame, and input the 'start minute' and 'end minute' for desired interval lengths
# result <- extract_intervals(your_data_frame, 56, 59)

result_fluxbot <- 
