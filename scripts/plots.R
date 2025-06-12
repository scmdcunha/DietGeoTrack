# -----------------------------------------------
# Analysis of Percent Identity as a function of Distance and Time
# Author: Sara Cunha
# Date: June 3, 2025
# -----------------------------------------------

# ----------------------------
# Install and load packages
# ----------------------------
packages <- c("ggplot2", "dplyr", "lubridate", "purrr", "patchwork", "scales")

installed <- rownames(installed.packages())
for (pkg in packages) {
  if (!(pkg %in% installed)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}

# ----------------------------
# Set working directory
# ----------------------------
setwd("/home/scmdcunha/metabarcoding_pipeline")

# ----------------------------
# Load and prepare data
# ----------------------------
dados <- read.delim("final_scores_clean.csv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

dados <- dados %>%
  mutate(
    Percent.Identity = as.numeric(Percent.Identity),   # Convert Percent.Identity to numeric
    Distance_km = as.numeric(Distance_km),             # Convert Distance_km to numeric
    EventDate = as.Date(EventDate),                     # Convert EventDate to Date type
    YearsSinceEvent = as.numeric(difftime(Sys.Date(), EventDate, units = "days")) / 365.25, # Calculate years since event
    Year = year(EventDate)                              # Extract year from EventDate
  ) %>%
  filter(!is.na(Percent.Identity), !is.na(Distance_km), !is.na(YearsSinceEvent))  # Remove rows with missing values

# Create groups based on fixed year intervals
dados <- dados %>%
  mutate(
    YearGroup = cut(Year,
                    breaks = c(-Inf, 2000, 2010, 2020, Inf),
                    labels = c("Before 2000", "2000-2010", "2011-2020", "2021+"),
                    right = TRUE,
                    include.lowest = TRUE)
  ) %>%
  filter(!is.na(YearGroup))

# ----------------------------
# General GLM models
# ----------------------------
# Fit GLM models with Gamma distribution and inverse link function
glm_dist <- glm(Percent.Identity ~ Distance_km, data = dados, family = Gamma(link = "inverse"))
glm_years <- glm(Percent.Identity ~ YearsSinceEvent, data = dados, family = Gamma(link = "inverse"))

# Create sequences of values for prediction
distance_seq <- seq(min(dados$Distance_km), max(dados$Distance_km), length.out = 100)
years_seq <- seq(min(dados$YearsSinceEvent), max(dados$YearsSinceEvent), length.out = 100)

# Predict Percent.Identity based on the GLM fits
predict_dist <- predict(glm_dist, newdata = data.frame(Distance_km = distance_seq), type = "response")
predict_years <- predict(glm_years, newdata = data.frame(YearsSinceEvent = years_seq), type = "response")

# ----------------------------
# General plots
# ----------------------------
# Scatter plot Percent Identity vs Distance with GLM fit line
plot1 <- ggplot(dados, aes(x = Distance_km, y = Percent.Identity)) +
  geom_point(color = "darkblue", size = 2) +
  geom_line(data = data.frame(Distance_km = distance_seq, Percent.Identity = predict_dist),
            aes(x = Distance_km, y = Percent.Identity),
            color = "red", linetype = "dashed", linewidth = 1.2) +
  labs(title = "Percent Identity vs Distance", x = "Distance (km)", y = "Percent Identity") +
  theme_minimal()

# Scatter plot Percent Identity vs Years Since Occurrence with GLM fit line
plot2 <- ggplot(dados, aes(x = YearsSinceEvent, y = Percent.Identity)) +
  geom_point(color = "darkgreen", size = 2) +
  geom_line(data = data.frame(YearsSinceEvent = years_seq, Percent.Identity = predict_years),
            aes(x = YearsSinceEvent, y = Percent.Identity),
            color = "red", linetype = "dashed", linewidth = 1.2) +
  labs(title = "Percent Identity vs Years Since Occurrence", x = "Years Since Occurrence", y = "Percent Identity") +
  theme_minimal()

print(plot1 + plot2)

# ----------------------------
# Boxplot and histogram
# ----------------------------
# Boxplot of Percent Identity grouped by YearGroup
plot3 <- ggplot(dados, aes(x = YearGroup, y = Percent.Identity)) +
  geom_boxplot(fill = "lightblue") +
  labs(title = "Percent Identity by Occurrence Year Group", x = "Year Group", y = "Percent Identity") +
  theme_minimal()

# Histogram of occurrence years distribution
plot4 <- ggplot(dados, aes(x = Year)) +
  geom_histogram(binwidth = 5, fill = "steelblue", color = "black") +
  labs(title = "Distribution of Occurrence Years", x = "Year of Occurrence", y = "Count") +
  theme_minimal()

print(plot3 + plot4)

# ----------------------------
# GLMs for fixed YearGroup categories
# ----------------------------
year_groups <- levels(dados$YearGroup)

for (yg in year_groups) {
  dados_sub <- filter(dados, YearGroup == yg)
  cat("Group:", yg, "- Number of rows:", nrow(dados_sub), "\n")

  min_year <- min(dados_sub$Year, na.rm = TRUE)
  max_year <- max(dados_sub$Year, na.rm = TRUE)

  if (is.finite(min_year) && is.finite(max_year) && max_year > min_year) {
    seq_years <- seq(min_year, max_year, length.out = 100)
    glm_fit <- try(glm(Percent.Identity ~ Year, data = dados_sub, family = Gamma(link = "inverse")), silent = TRUE)

    if (inherits(glm_fit, "try-error")) {
      p <- ggplot(dados_sub, aes(x = Year, y = Percent.Identity)) +
        geom_point(color = "darkgreen", size = 2) +
        labs(title = paste("Percent Identity vs Year -", yg), x = "Year", y = "Percent Identity") +
        theme_minimal() +
        xlim(min_year, max_year)
    } else {
      pred <- predict(glm_fit, newdata = data.frame(Year = seq_years), type = "response")
      pred_df <- data.frame(Year = seq_years, Percent.Identity = pred)

      p <- ggplot(dados_sub, aes(x = Year, y = Percent.Identity)) +
        geom_point(color = "darkgreen", size = 2) +
        geom_line(data = pred_df, aes(x = Year, y = Percent.Identity),
                  color = "red", linetype = "dashed", linewidth = 1.2) +
        labs(title = paste("Percent Identity vs Year -", yg), x = "Year", y = "Percent Identity") +
        theme_minimal() +
        xlim(min_year, max_year)
    }

    print(p)
    readline(prompt = "Press [enter] to continue...")
  } else {
    message("Invalid year range for group ", yg)
  }
}

# ----------------------------
# Mean and median Percent Identity by year
# ----------------------------
media_filtrada <- dados %>%
  group_by(Year) %>%
  summarise(
    N = n(),
    MeanPercentID = mean(Percent.Identity, na.rm = TRUE),
    MedianPercentID = median(Percent.Identity, na.rm = TRUE)
  ) %>%
  ungroup()

# Filter years with >= 30 occurrences (reliable sample size)
media_confiavel <- media_filtrada %>% filter(N >= 30)

# Year with highest mean Percent Identity
print(media_confiavel %>% arrange(desc(MeanPercentID)) %>% slice(1))

# Top 10 years with most data points
print(media_confiavel %>% arrange(desc(N)) %>% head(10))

# Visualization of mean Percent Identity per year with number of occurrences as point size
ggplot(media_confiavel, aes(x = Year, y = MeanPercentID)) +
  geom_line(color = "darkblue") +
  geom_point(aes(size = N), color = "red", alpha = 0.6) +
  scale_size_continuous(name = "Number of Occurrences") +
  labs(title = "Mean Percent Identity per Year",
       x = "Year", y = "Mean Percent Identity") +
  theme_minimal()

# ----------------------------
# Top 20 Orders by frequency
# ----------------------------
# Check if "Order" column exists
if (!"Order" %in% names(dados)) {
  stop("The 'Order' column is not present in the final_scores_clean.csv file")
}

# Count frequency of Orders
order_freq <- dados %>%
  filter(!is.na(Order) & Order != "") %>%
  group_by(Order) %>%
  summarise(Count = n()) %>%
  ungroup() %>%
  mutate(Percentage = 100 * Count / sum(Count))

# Get top 20 Orders
top_orders <- order_freq %>%
  slice_max(order_by = Count, n = 20)

# Mark the rest as "Others"
dados_orders <- dados %>%
  filter(!is.na(Order) & Order != "") %>%
  mutate(OrderGrouped = ifelse(Order %in% top_orders$Order, Order, "Others"))

# Recount with "Others"
orders_final <- dados_orders %>%
  group_by(OrderGrouped) %>%
  summarise(Count = n()) %>%
  ungroup() %>%
  mutate(Percentage = 100 * Count / sum(Count)) %>%
  arrange(desc(Percentage))

# Bar plot of Order percentages
ggplot(orders_final, aes(x = reorder(OrderGrouped, -Percentage), y = Percentage)) +
  geom_bar(stat = "identity", fill = "cornflowerblue") +
  labs(title = "Percentage of Top 20 Orders (+ 'Others')",
       x = "Order",
       y = "Percentage of Occurrences") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ----------------------------
# GLM models by Order with plots
# ----------------------------
# Loop over top Orders to fit GLM Percent.Identity ~ Distance_km and plot
for (ord in unique(dados_orders$OrderGrouped)) {
  data_ord <- filter(dados_orders, OrderGrouped == ord)

  if (nrow(data_ord) < 10) {
    message("Skipping ", ord, " due to insufficient data")
    next
  }

  glm_ord <- try(glm(Percent.Identity ~ Distance_km, data = data_ord, family = Gamma(link = "inverse")), silent = TRUE)

  if (inherits(glm_ord, "try-error")) {
    message("GLM failed for ", ord)
    next
  }

  dist_seq <- seq(min(data_ord$Distance_km), max(data_ord$Distance_km), length.out = 100)
  pred_ord <- predict(glm_ord, newdata = data.frame(Distance_km = dist_seq), type = "response")
  pred_df <- data.frame(Distance_km = dist_seq, Percent.Identity = pred_ord)

  p <- ggplot(data_ord, aes(x = Distance_km, y = Percent.Identity)) +
    geom_point(color = "darkorange", size = 2) +
    geom_line(data = pred_df, aes(x = Distance_km, y = Percent.Identity),
              color = "blue", linetype = "dashed", linewidth = 1.2) +
    labs(title = paste("Percent Identity vs Distance - Order:", ord),
         x = "Distance (km)",
         y = "Percent Identity") +
    theme_minimal()

  print(p)
  readline(prompt = "Press [enter] to continue...")
}
