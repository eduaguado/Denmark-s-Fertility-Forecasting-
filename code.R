library(tsibble)
library(lubridate)
library(fable)
library(feasts)
library(dplyr)
library(tidyr)
library(ggplot2)
library(urca)
library(forecast)


# Different Datasets Used
# 1 raw
raw_data <- dk_births_2010_2024

raw_data <- dk_births_2010_2024 %>%
  mutate(
    Month = dmy(Month),           # convert to Date
    Month = yearmonth(Month)      # convert to yearmonth class
  ) %>%
  as_tsibble(index = Month)

# raw_data plot
ggplot(raw_data, aes(x = Month, y = Births)) +
  geom_line(color = "black") +
  labs(
    title = "Monthly Births in Denmark (2010–2025)",
    x = "Month",
    y = "Number of Births (in thousands)"
  ) +
  theme_minimal()


# Transformation = bc transformed raw_data_bc
raw_data %>%
  features(Births, features = guerrero)
# Value is -0.9
raw_data %>%
  autoplot(box_cox(Births, -0.9)) + ggtitle("Box-Cox Transformed Monthly Births in Denmark (2010–2025)")
# make a new column with the transformation
raw_data %>%
  mutate(births_bc = box_cox(Births, 0.933)) -> raw_data_bc


# diff p = 1
raw_data_diff1 <- raw_data_bc %>%
  mutate(births_diff = difference(births_bc)) %>%
  filter(!is.na(births_diff)) %>%
  select(-Births, -births_bc)


# 1. STL
raw_data %>%
  model(STL(Births)) %>%
  components() %>%
  autoplot() +
  labs(
    title = "STL Decomposition of Monthly Births in Denmark (2010-2025)",
    y = "Births (in thousands)",
    x = "Month"
  )

# Seasonality (there is a clear seasonal pattern)
raw_data %>%
  gg_season(Births) + labs(title = "Seasonal Plot of Monthly Births in Denmark (2010-2025)")


# 2. ACF-PACF (Seasonal patterns)
# ACF
raw_data%>%
  ACF(Births) %>%
  autoplot() + labs(title = "ACF")

#PACF
raw_data%>%
  PACF(Births) %>%
  autoplot() + labs(title = "PACF")


# Box pierce test for autocorrelation (RH0, there is autocorrelation)
raw_data %>%
  features(Births, box_pierce, lag=24, dof=0)
# there is autocorrelation confirmed



# 3. Stationarity for diff1 model
# ADF. !! We need to differentiate on p =1 to be stationary!!
summary(ur.df(raw_data_bc$births_bc, type = "none"))
# At all conf levels, we reject H0 -> Stationary around 0 mean
# Without difference we didn't rejected it! Need diff1 to stationary

summary(ur.df(raw_data_bc$births_bc, type = "trend"))
# We reject H0 -> Series is stationary with deterministic trend
# Without difference we also rejected it! No need of diff1 but better t-value

summary(ur.df(raw_data_bc$births_bc, type = "drift"))
# At all conf levels, we reject H0 -> Series is stationary with a drift
# Without difference we also rejected it! No need of diff1 but better t-value

# KPSS
summary(ur.kpss(raw_data_bc$births_bc, type = "mu"))
# Not Reject H0 -> Series is stationary around a mean
# Without difference we didn't rejected it! No need of diff1 but better t-value

summary(ur.kpss(raw_data_bc$births_bc, type = "tau"))
# Not Reject H0 -> Series is stationary around a trend
# Without difference we didn't rejected it! No need diff1 to stationary but better t-value



# Stationarity test 1 

# Add time variable
raw_data_bc$t <- 1:nrow(raw_data_bc)

# Now allow for changing trend over time
trend_model <- breakpoints(births_bc ~ t, data = raw_data_bc)
summary(trend_model)

# Dates of breakpoints
raw_data_bc$Month[trend_model$breakpoints]




# Stationarity test 2
raw_data_bc <- raw_data_bc %>%
  mutate(l1.Births = lag(births_bc),
         l12.Births = lag(births_bc,12))

# QLR test
qlr <- Fstats(births_bc ~ l1.Births + l12.Births, data = as.ts(raw_data_bc), from = 0.15)
test <- sctest(qlr, type = "supF")
test
# Plot it
plot(qlr, alpha = 0.05, main = "F Statistics")
lines(breakpoints(qlr))

# Dates of breakpoints
breaks <- breakpoints(births_bc ~ l1.Births + l12.Births, data = raw_data_bc)$breakpoints
raw_data_bc$Month[breaks]

# Create dummies and eliminate useless variables (t, l1.Births and l12.Births)

# Create post_campaign dummies
raw_data_bc$post_campaign <- as.numeric(raw_data_bc$Month >= yearmonth("2015 May") & raw_data_bc$Month <= yearmonth("2021 Nov"))
raw_data_bc$post_covid    <- as.numeric(raw_data_bc$Month > yearmonth("2021 Nov"))
raw_data_bc <- raw_data_bc %>%
  select(-t, -l1.Births, -l12.Births)


# MODEL SELECTION (Best sarimax_111_011)

# A) Auto-Optimal SARIMA (without structural break dummies)
births_ts <- ts(raw_data_bc$births_bc, start = c(2010, 1), frequency = 12)

sarima_auto_fit <- auto.arima(
  births_ts,
  seasonal = TRUE,
  stepwise = FALSE,
  approximation = FALSE,
  trace = TRUE
)

summary(sarima_auto_fit)
# The best model is ARIMA(1,0,1)(2,1,0)[12] -197, -181 

# B) Auto-Optimal SARIMAX (with structural breaks)
# Define the external regressors matrix with your two dummies
xreg <- as.matrix(raw_data_bc[, c("post_campaign", "post_covid")])

# Run auto.arima including the external regressors
sarimax_auto_fit <- auto.arima(
  births_ts,
  xreg = xreg,
  seasonal = TRUE,
  stepwise = FALSE,
  approximation = FALSE,
  trace = TRUE
)

# View the model summary
summary(sarimax_auto_fit)
# The best model is SARIMAX(1,1,1)(2,1,0)[12]  -197, -175


# C) Manual SARIMAX Selection
models <- raw_data_bc %>%
  model(
    sarimax_010_010 = ARIMA(births_bc ~ post_campaign + post_covid + pdq(0,1,0) + PDQ(0,1,0)),
    
    sarimax_110_010 = ARIMA(births_bc ~ post_campaign + post_covid + pdq(1,1,0) + PDQ(0,1,0)),
    
    sarimax_011_010 = ARIMA(births_bc ~ post_campaign + post_covid + pdq(0,1,1) + PDQ(0,1,0)),
    
    sarimax_111_010 = ARIMA(births_bc ~ post_campaign + post_covid + pdq(1,1,1) + PDQ(0,1,0)),
    
    sarimax_111_110 = ARIMA(births_bc ~ post_campaign + post_covid + pdq(1,1,1) + PDQ(1,1,0)),
    
    sarimax_111_011 = ARIMA(births_bc ~ post_campaign + post_covid + pdq(1,1,1) + PDQ(0,1,1)),
    
    sarimax_111_111 = ARIMA(births_bc ~ post_campaign + post_covid + pdq(1,1,1) + PDQ(1,1,1))
  )


report(models)
# Best model sarimax_111_011 -220, -201 (BEST MODEL ACCORDING TO AIC AND BIC)
# Alternative models sarimax_111_111 -218, -196

sarimax_111_011 <- Arima(
  births_ts,
  order = c(1,1,1),
  seasonal = c(0,1,1),
  xreg = xreg
)

summary(sarimax_111_011)

sarimax_111_111 <- Arima(
  births_ts,
  order = c(1,1,1),
  seasonal = c(1,1,1),
  xreg = xreg
)


summary(sarimax_111_111)

# LJUNG-BOX TEST ANALYSIS (Best sarimax_111_011)
models <- list(sarima_101_210 = sarima_auto_fit, sarimax_111_210 = sarimax_auto_fit, sarimax_111_011 = sarimax_111_011, sarimax_111_111 = sarimax_111_111)

# Run Ljung-Box test for each model and collect p-values
ljung_results <- lapply(names(models), function(name) {
  mod <- models[[name]]
  res <- residuals(mod)
  
  # Ljung-Box test: HE UTILITZAT LAG 36!!!!
  test <- Box.test(res, lag = 36, type = "Ljung-Box", fitdf = length(coef(mod)))
  
  data.frame(
    Model = name,
    Ljung_Box_Statistic = test$statistic,
    Ljung_Box_p_value = test$p.value
  )
})

# Combine into a single data frame
ljung_results_df <- bind_rows(ljung_results)

# Print the results table
print(ljung_results_df)
# Best model, the SARIMAX manually selected (sarimax_111_011)!!!


# RESIDUAL ANALYSIS (sense distribucio de residuals)
tsdisplay(residuals(models$sarima_101_210), main = "SARIMA(101)(210) - Residual Diagnostics")
tsdisplay(residuals(models$sarimax_111_210), main = "SARIMAX(111)(011) - Residual Diagnostics")
tsdisplay(residuals(models$sarimax_111_011), main = "SARIMAX(111)(011) - Residual Diagnostics")
tsdisplay(residuals(models$sarimax_111_111), main = "SARIMAX(111)(111) - Residual Diagnostics")



# RESIDUAL ANALYSIS - NOMES DISTRIBUCIÓ
res <- residuals(models[["sarimax_111_011"]])
  
# Display tsdiagnostics
tsdisplay(res, main = "SARIMAX(1,1,1)(0,1,1) - Residual Distribution"))
  
# Plot histogram
hist_df <- data.frame(res)
print(
  ggplot(hist_df, aes(x = res)) +
    geom_histogram(bins = 30, fill = "grey", color = "black") +
    labs(title = paste("SARIMAX(1,1,1)(0,1,1) - Residual Distribution"), x = "Residuals", y = "Frequency") +
    theme_minimal())

# 2025 FORECAST
library(forecast)
library(ggplot2)
library(dplyr)

# Step 1: Create future exogenous variables
future_xreg <- data.frame(
  post_campaign = rep(1, 12),  # Replace with actual future values if available
  post_covid = rep(1, 12)
) %>% 
  mutate(across(everything(), as.numeric)) %>%
  as.matrix()

# Step 2: Forecast using fitted SARIMAX model
fcast <- forecast(sarimax_111_011, xreg = future_xreg, h = 12)

# Step 3: Format forecast output
f_df <- as.data.frame(fcast) %>%
  rename(Forecast = `Point Forecast`, Lo80 = `Lo 80`, Hi80 = `Hi 80`, Lo95 = `Lo 95`, Hi95 = `Hi 95`)
f_df$Date <- as.Date(time(fcast$mean))

# Step 4: Format historical data
hist_df <- data.frame(
  Date = as.Date(time(births_ts)),
  Births = as.numeric(births_ts)
)

# Step 5: Plot forecast with legend
ggplot() +
  geom_line(data = hist_df, aes(x = Date, y = Births, color = "Observed")) +
  geom_line(data = f_df, aes(x = Date, y = Forecast, color = "Forecast"), linetype = "dashed") +
  geom_ribbon(data = f_df, aes(x = Date, ymin = Lo95, ymax = Hi95, fill = "95% PI"), alpha = 0.3) +
  geom_ribbon(data = f_df, aes(x = Date, ymin = Lo80, ymax = Hi80, fill = "80% PI"), alpha = 0.4) +
  scale_color_manual(name = "Lines", values = c("Observed" = "black", "Forecast" = "red")) +
  scale_fill_manual(name = "Prediction Intervals", values = c("95% PI" = "lightblue", "80% PI" = "blue")) +
  labs(
    title = "SARIMAX(1,1,1)(0,1,1) – Forecast of Births in Denmark (2025)",
    x = "Date", y = "Number of Births"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )


