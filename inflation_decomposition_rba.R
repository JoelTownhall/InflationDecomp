# Inflation Decomposition - RBA Method (Hybrid Approach)
# HFCE for Signals (Supply/Demand) + CPI for Weights/Prices

options(repos = c(CRAN = "http://cran.us.r-project.org"))

pkgs <- c("readabs", "readrba", "tidyverse", "vars", "zoo", "lubridate", "ggplot2", "scales")
new_pkgs <- pkgs[!(pkgs %in% installed.packages()[,"Package"])]
if(length(new_pkgs)) install.packages(new_pkgs)
invisible(lapply(pkgs, library, character.only = TRUE))

dir.create("outputs", showWarnings = FALSE)

# --- 1. Data Gathering ---

# HFCE (5206.0 Table 8)
message("Downloading HFCE data...")
hfce_raw <- read_abs("5206.0", tables = 8, check_local = FALSE)

# CPI Groups (6401.0 Table 3)
message("Downloading CPI Group data...")
cpi_groups_raw <- read_abs("6401.0", tables = 3, check_local = FALSE)

# CPI Headline (RBA)
message("Downloading Headline CPI...")
cpi_headline_rba <- read_rba(series_id = "GCPIAG") %>%
  mutate(date = floor_date(date, "quarter"))

# --- 2. CPI Processing (Groups) ---

# Define Standard CPI Groups and 2023 Weights (Approx)
# Weights source: ABS 6470.0 (2023)
cpi_weights <- tibble(
  Group = c("Food and non-alcoholic beverages", 
            "Alcohol and tobacco", 
            "Clothing and footwear", 
            "Housing", 
            "Furnishings, household equipment and services", 
            "Health", 
            "Transport", 
            "Communication", 
            "Recreation and culture", 
            "Education", 
            "Insurance and financial services"),
  Weight = c(17.15, 7.70, 3.06, 21.74, 8.75, 6.16, 11.55, 2.41, 12.55, 4.36, 4.57)
)

# Extract CPI Indices for these Groups
cpi_indices <- cpi_groups_raw %>%
  filter(grepl("Index Numbers", series, ignore.case = TRUE),
         grepl("Australia", series, ignore.case = TRUE)) %>%
  # Match Group Names
  filter(series %in% paste0("Index Numbers ;  ", cpi_weights$Group, " ;  Australia ;")) %>%
  mutate(
    Group = str_remove(series, "Index Numbers ;  "),
    Group = str_remove(Group, " ;  Australia ;")
  ) %>%
  dplyr::select(date, Group, CPI_Index = value) %>%
  arrange(Group, date) %>%
  group_by(Group) %>%
  mutate(
    CPI_Inflation_Q = (CPI_Index / lag(CPI_Index) - 1) * 100,
    CPI_Inflation_YoY = (CPI_Index / lag(CPI_Index, 4) - 1) * 100
  ) %>%
  ungroup() %>%
  left_join(cpi_weights, by = "Group")

# --- 3. HFCE Processing (Aggregation to match CPI) ---

hfce_clean <- hfce_raw %>% 
  filter(frequency == "Quarter", 
         series_type == "Seasonally Adjusted") %>%
  dplyr::select(date, series, value, series_id) 

# Map HFCE specific series to CPI Groups
# We need both CP (Price) and CVM (Qty) for:
# 1. Food -> Food
# 2. Cigarettes + Alcohol -> Alcohol
# 3. Clothing -> Clothing
# 4. Rent + Electricity -> Housing
# 5. Furnishings -> Furnishings
# 6. Health -> Health
# 7. Transport -> Transport
# 8. Communications -> Communication
# 9. Recreation + Hotels -> Recreation
# 10. Education -> Education
# 11. Insurance + Other -> Insurance

# Helper to sum series by pattern
aggregate_hfce <- function(data, pattern_list, group_name) {
  # Find CP series
  cp_series <- data %>% 
    filter(grepl("Current prices", series, ignore.case = TRUE),
           !grepl("Total", series)) %>%
    filter(grepl(paste(pattern_list, collapse="|"), series, ignore.case = TRUE))
  
  # Find CVM series
  cvm_series <- data %>% 
    filter(grepl("Chain volume measures", series, ignore.case = TRUE),
           !grepl("Total", series)) %>%
    filter(grepl(paste(pattern_list, collapse="|"), series, ignore.case = TRUE))
  
  # Sum
  cp_sum <- cp_series %>% group_by(date) %>% summarise(CP = sum(value, na.rm=TRUE))
  cvm_sum <- cvm_series %>% group_by(date) %>% summarise(CVM = sum(value, na.rm=TRUE))
  
  inner_join(cp_sum, cvm_sum, by = "date") %>% mutate(Group = group_name)
}

# Define Mappings
hfce_mapped <- bind_rows(
  aggregate_hfce(hfce_clean, c("Food"), "Food and non-alcoholic beverages"),
  aggregate_hfce(hfce_clean, c("Cigarettes", "Alcoholic"), "Alcohol and tobacco"),
  aggregate_hfce(hfce_clean, c("Clothing"), "Clothing and footwear"),
  aggregate_hfce(hfce_clean, c("Rent", "Electricity"), "Housing"),
  aggregate_hfce(hfce_clean, c("Furnishings"), "Furnishings, household equipment and services"),
  aggregate_hfce(hfce_clean, c("Health"), "Health"),
  aggregate_hfce(hfce_clean, c("Transport"), "Transport"),
  aggregate_hfce(hfce_clean, c("Communications"), "Communication"),
  aggregate_hfce(hfce_clean, c("Recreation", "Hotels"), "Recreation and culture"),
  aggregate_hfce(hfce_clean, c("Education"), "Education"),
  aggregate_hfce(hfce_clean, c("Insurance", "Other goods"), "Insurance and financial services")
)

# Calculate Implicit Metrics for VAR
hfce_signals <- hfce_mapped %>%
  mutate(
    Price_Index = (CP / CVM) * 100,
    Quantity_Index = CVM
  ) %>%
  arrange(Group, date) %>%
  group_by(Group) %>%
  mutate(
    d_ln_P = log(Price_Index) - lag(log(Price_Index)),
    d_ln_Q = log(Quantity_Index) - lag(log(Quantity_Index))
  ) %>%
  ungroup() %>%
  filter(!is.na(d_ln_P))

# --- 4. Rolling VAR (Signal Extraction) ---

run_rolling_var <- function(sub_data, window_size = 40) {
  sub_data <- sub_data %>% arrange(date)
  n <- nrow(sub_data)
  results <- list()
  if(n <= window_size) return(NULL)
  
  for(i in (window_size + 1):n) {
    window_start <- i - window_size + 1
    slice_data <- sub_data[window_start:i, ]
    current_date <- slice_data$date[window_size] 
    
    y <- slice_data %>% dplyr::select(d_ln_P, d_ln_Q) %>% as.matrix()
    if(any(is.na(y))) next
    
    tryCatch({
      var_model <- VAR(y, p = 4, type = "const")
      resids <- residuals(var_model)
      last_resid <- resids[nrow(resids), ]
      sigma <- apply(resids, 2, sd)
      t_stats <- last_resid / sigma
      
      threshold <- 0.25
      is_ambiguous <- all(abs(t_stats) < threshold)
      label <- "Ambiguous"
      if(!is_ambiguous) {
        if(sign(last_resid[1]) == sign(last_resid[2])) label <- "Demand" else label <- "Supply"
      }
      results[[length(results) + 1]] <- tibble(date = current_date, label = label)
    }, error = function(e) return(NULL))
  }
  bind_rows(results)
}

message("Running VARs on HFCE groups...")
hfce_labels <- hfce_signals %>%
  group_split(Group) %>%
  map_dfr(function(df) {
    res <- run_rolling_var(df)
    if(!is.null(res)) res$Group <- unique(df$Group)
    return(res)
  })

# --- 5. Hybrid Calculation (CPI Data + HFCE Labels) ---

final_data <- cpi_indices %>%
  inner_join(hfce_labels, by = c("date", "Group")) %>%
  mutate(
    Weighted_Inflation = CPI_Inflation_Q * Weight
  ) %>%
  group_by(date, label) %>%
  summarise(
    Weighted_Sum = sum(Weighted_Inflation, na.rm=TRUE),
    Total_Weight = sum(Weight, na.rm=TRUE),
    .groups = "drop"
  ) %>%
  # Normalize to 100% of CPI (since weights sum to ~100 but we might miss some residuals)
  group_by(date) %>%
  mutate(
    # Re-scale to ensure parts sum to the Total Modeled CPI (which tracks Headline)
    # Actually, let's just use the sum. 
    # Contribution = (Points / Total Weights) approx.
    Contribution = Weighted_Sum / sum(cpi_weights$Weight) 
  ) %>%
  ungroup()

# --- 6. Visualization Prep ---

# Pivot for Chart
plot_data_q <- final_data %>%
  dplyr::select(date, label, Contribution) %>%
  pivot_wider(names_from = label, values_from = Contribution, values_fill = 0) %>%
  mutate(Total_Q = Demand + Supply + Ambiguous)

# Calculate Year-Ended
plot_data_yoy <- plot_data_q %>%
  arrange(date) %>%
  mutate(
    Demand_YoY = rollsum(Demand, 4, fill = NA, align = "right"),
    Supply_YoY = rollsum(Supply, 4, fill = NA, align = "right"),
    Ambiguous_YoY = rollsum(Ambiguous, 4, fill = NA, align = "right"),
    Total_YoY = rollsum(Total_Q, 4, fill = NA, align = "right")
  ) %>%
  # Convert to Percentage?
  # Note: CPI_Inflation_Q is already %, Weight is absolute.
  # If Weight=17, Inf=1%, Weighted=17. Sum=100. Contrib = 17/100 = 0.17%.
  # So data is already in %.
  filter(!is.na(Total_YoY))

# Get Headline
headline <- cpi_headline_rba %>%
  mutate(CPI_YoY = (value / lag(value, 4) - 1) * 100) %>%
  dplyr::select(date, CPI_YoY)

combined <- plot_data_yoy %>%
  left_join(headline, by = "date") %>%
  filter(date >= as.Date("2010-01-01"))

# Last Value
last_obs <- combined %>% filter(!is.na(CPI_YoY)) %>% tail(1)
last_val <- round(last_obs$CPI_YoY, 1)
last_date <- last_obs$date

# Plot
plot_long <- combined %>%
  dplyr::select(date, Demand_YoY, Supply_YoY, Ambiguous_YoY) %>%
  pivot_longer(cols = -date, names_to = "Component", values_to = "Value") %>%
  mutate(Component = gsub("_YoY", "", Component))

p <- ggplot() +
  geom_col(data = plot_long, aes(x = date, y = Value, fill = Component), width = 90, alpha = 0.9) +
  geom_line(data = combined, aes(x = date, y = CPI_YoY, color = "Headline CPI"), size = 1.2) +
  scale_fill_manual(values = c("Demand" = "#FF5A5F", "Supply" = "#007A87", "Ambiguous" = "#7B7B7B")) +
  scale_color_manual(values = c("Headline CPI" = "black")) +
  annotate("rect", xmin = as.Date("2020-03-01"), xmax = as.Date("2021-12-31"), ymin = -Inf, ymax = Inf, alpha = 0.15, fill = "yellow") +
  annotate("text", x = as.Date("2020-09-01"), y = 5, label = "COVID-19", vjust = -0.5, size = 3.5, fontface = "bold") +
  annotate("text", x = last_date, y = as.numeric(last_val), label = paste0(last_val, "%"), hjust = -0.3, vjust = 0.5, fontface = "bold", size = 4) +
  labs(title = "Decomposition of Australian CPI (RBA Hybrid Method)", 
       subtitle = "Signals from HFCE, Impacts calculated on CPI Groups",
       y = "Year-ended Inflation (%)", x = NULL, caption = "Source: ABS, RBA. Fixed weights approximation.") +
  theme_minimal() + theme(legend.position = "bottom") + scale_x_date(date_breaks = "2 years", date_labels = "%Y")

ggsave("outputs/inflation_decomposition.png", plot = p, width = 10, height = 6, bg = "white")
write_csv(combined, "outputs/inflation_decomposition_data_rba.csv")
message("Done.")
