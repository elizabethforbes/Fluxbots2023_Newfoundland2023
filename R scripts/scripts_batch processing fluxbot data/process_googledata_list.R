# each fluxbot sensor will return a df with two columns: one a 2min interval of timestamp and the 
# other 2min of data inside of square brackets per cell.

# this function processes a list of dataframes and transposes the within-cell
# data into long-form. It also converts the UNIX timestamps to POSIXct-form
# timestamps for easier analysis.

library(dplyr)
library(tidyr)

process_data <- function(df_list, col1, col2) {
  for (i in seq_along(df_list)) {
    df <- df_list[[i]]
    processed_df <- df %>%
      separate_rows(c(!!rlang::sym(col1), !!rlang::sym(col2))) %>%
      select({{col1}}, {{col2}}) %>%
      rename(UNIX = {{col1}}) %>%
      mutate(
        UNIX = as.numeric(UNIX),
        timestamp = as.POSIXct(UNIX, origin = "1970-01-01"),
        {{col2}} := as.numeric({{col2}})
      )
    
    # Generate filename based on index
    filename <- paste0("output_csv/dataframe_", i, ".csv")
    
    # Export processed dataframe as CSV
    write.csv(processed_df, file = filename, row.names = FALSE)
  }
}

# Example usage:
# Assuming df_list is a list of dataframes
# df_list <- list(data.frame(col1 = c("1,2,3", "4,5", "6"), col2 = c("A,B,C", "D,E", "F")), 
#                 data.frame(col1 = c("10,20,30", "40,50", "60"), col2 = c("X,Y,Z", "U,V", "W")))
# 
# process_data(df_list, "col1", "col2")