merge_by_timestamp <- function(df1, df2) {
  merged_df <- merge(df1, df2, by = "timestamp", all = TRUE)
  return(merged_df)
}

merge_lists_by_timestamp <- function(list_of_df_lists, matching_string) {
  merged_results <- list() # list to store final merged dataframes
  
  # Merge dataframes within each list
  for (df_list in list_of_df_lists) {
    # Filter dataframes by matching string in file names
    matching_dfs <- Filter(function(df) grepl(matching_string, names(df)), df_list)
    
    # If there are multiple matching dataframes, merge them by timestamp
    if (length(matching_dfs) > 1) {
      merged_df <- Reduce(merge_by_timestamp, matching_dfs)
      merged_results[[names(matching_dfs[[1]])]] <- merged_df
    } else if (length(matching_dfs) == 1) {  # If there's only one matching dataframe, no need to merge
      merged_results[[names(matching_dfs[[1]])]] <- matching_dfs[[1]]
    }
  }
  
  return(merged_results)
}

# Example usage:
# Assuming list_of_df_lists is a list containing multiple lists of dataframes
# list_of_df_lists <- list(
#   list(data.frame(timestamp = seq.POSIXt(as.POSIXct("2024-01-01"), as.POSIXct("2024-01-03"), by = "day"), data1 = c(1, 2, 3), data2 = c(10, 20, 30)), 
#        data.frame(timestamp = seq.POSIXt(as.POSIXct("2024-01-02"), as.POSIXct("2024-01-04"), by = "day"), data1 = c(4, 5, 6), data2 = c(40, 50, 60))),
#   list(data.frame(timestamp = seq.POSIXt(as.POSIXct("2024-01-03"), as.POSIXct("2024-01-05"), by = "day"), data1 = c(7, 8, 9), data2 = c(70, 80, 90)), 
#        data.frame(timestamp = seq.POSIXt(as.POSIXct("2024-01-04"), as.POSIXct("2024-01-06"), by = "day"), data1 = c(10, 11, 12), data2 = c(100, 110, 120)))
# )
# 
# matching_string <- "data"
# 
# merged_result <- merge_lists_by_timestamp(list_of_df_lists, matching_string)
