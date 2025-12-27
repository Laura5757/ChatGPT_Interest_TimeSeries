library(readxl)
library(tidyverse)
library(dplyr)
library(tseries)
library(epitools)
library(GGally)
library(zoo)
library(xts)
library(fpp2)
library(aTSA)
library(tsibble)
library(ggthemes)
library(gridExtra)
library(RcppArmadillo)
library(TSA)
library(lubridate)
library(fastDummies)
library(fable)
library(fabletools)
library(prophet)
library(fable.prophet)
library(forecast)
library(urca)


#Set R in English
Sys.setlocale("LC_TIME", "en_US.UTF-8")

#Chat-GPT search interest index from google trends 
#data <- read.csv("chatgpt_searchinterest.csv")
data <- read_csv("D:/Maria/Desktop/KTU/SA/Project/chatgpt_searchinterest.csv")
data$Day <- as.Date(data$Day, format = "%d/%m/%Y") #as date
ts_data <- ts(data$ChatGPT_interest) #to make a ts object

# 1. Plot -----------------------------------------------------------------
autoplot(ts_data) +
  ggtitle("Chat-GPT search interest (from Google trends)") +
  xlab("Time") +
  ylab("Interest index")


# 2. Statistical properties -----------------------------------------------
#Mean
summary(ts_data)
sd(ts_data, na.rm = TRUE)
var(ts_data, na.rm = TRUE)
min(ts_data, na.rm = TRUE)
max(ts_data, na.rm = TRUE)
ts_mean <- mean(ts_data, na.rm = TRUE)
ts_mean
plot.ts(ts_data, main = "Chat-GPT search interest (from Google trends)")
abline(h = ts_mean, col="red", lwd=2)

#Auto_covariance and auto-correlation function

par(mfrow = c(1, 3)) 
acf(ts_data, type = "covariance", main = "Auto-covariance Function")
acf(ts_data, type = "correlation", main = "Auto-correlation Function")
acf(ts_data, type = "partial", main = "Partial autocorrelation function")
par(mfrow = c(1, 1)) 

#Spectral analysis
#Periodogram
par(mfrow = c(1, 3))
stats::spectrum(ts_data, log="no", main = "Periodogram")
#Spectrum
stats::spectrum(ts_data, method="ar", log="dB")
spec.ar(ts_data, order=7, log="dB")
par(mfrow = c(1, 1))

# 3. Identify if seasonality exists. --------------------------------------
#Dickey-Fuller test
aTSA::adf.test(ts_data, nlag = NULL, output = TRUE)

#Look for the optimal k from the previous results of Dickey-Fuller test
#Criteria
test_AIC <- ur.df(ts_data, type = "drift", selectlags = "AIC", lags = 40)
test_BIC <- ur.df(ts_data, type = "drift", selectlags = "BIC", lags = 40)
summary(test_AIC)
print(test_AIC@lags)
print(test_BIC@lags)


# 1. Selección de Lag con Criterio AIC (Modelo Tipo 3: con tendencia)
test_trend_AIC <- ur.df(ts_data, type = "trend", selectlags = "AIC", lags = 40)

# 2. Ejecución de la prueba final (usando el lag seleccionado, que esperamos sea 40)
# Si AIC y BIC seleccionan k=40, puedes usar la siguiente línea directamente:
test_final_trend <- ur.df(ts_data, type = "trend", lags = 40)
print(test_trend_AIC@lags)
summary(test_final_trend)


#We can't reject the null hypothesis, therefore, we can consider the time series as non-stationary

#Seasonality
ts_data_w <- ts(as.numeric(ts_data), frequency = 7) #to see weekly seasonality
ts_data_a <- ts(as.numeric(ts_data), frequency = 365.25) #to see annual seasonality

#Weekly
fit_add_w <- decompose(ts_data_w, type = "additive")
fit_mult_w <- decompose(ts_data_w, type = "multiplicative")
plot(fit_add_w)
plot(fit_mult_w)
#From this, multiplicative looks better since the remainder has a lower variation in multiplicative than in additive

#Annual
fit_add_a <- decompose(ts_data_a, type = "additive")
fit_mult_a <- decompose(ts_data_a, type = "multiplicative")
plot(fit_add_a)
plot(fit_mult_a)
#As the weekly seasonality, remainder variation is lower, around 1 than the additive, which varies around -30 to 30

#To see if residuals are autocorrelated
#residuals
res_mult_w <- fit_mult_w$random
res_mult_a <- fit_mult_a$random
#Test Ljung-Box
Box.test(res_mult_w,  lag = frequency(ts_data_w), type="Ljung-Box")
Box.test(res_mult_a, lag = frequency(ts_data_a), type="Ljung-Box")
#The null hypothesis is rejected, then, both remainders are correlated
acf(na.omit(res_mult_w), main="ACF Remainder Weekly")
acf(na.omit(res_mult_a), main="ACF Remainder Annual")
#PLotting the autocorrelation, ir can be observed better the autocorrelation

# 4. Decompose time series  -----------------------------------------------
#Then it will be use multiple seasonality
ts_msts <- msts(as.numeric(ts_data), 
                seasonal.periods = c(7, 365.25),
                start = c(year(min(data$Day)), yday(min(data$Day))))
mstl_fit <- mstl(ts_msts, lambda = "auto")
autoplot(mstl_fit) + 
  ggtitle("MSTL Decomposition (Weekly + Annual Seasonality)") +
  theme_minimal()
#From this plot it can be observed that compared with the other components, then the annualy seasonality is stronger

#Components of decomposition
mstl_components <- as.data.frame(mstl_fit)
mstl_trend <- mstl_components$Trend
mstl_seasonal_w <- mstl_components$Seasonal7
mstl_seasonal_y <- mstl_components$Seasonal365.25
mstl_remainder <- mstl_components$Remainder

#Test Ljung-Box
Box.test(mstl_remainder,  lag = min(30, length(mstl_remainder)/2), type="Ljung-Box")
acf(mstl_remainder, main="ACF MSTL Remainder")

par(mfrow = c(2, 1))
plot(mstl_seasonal_w, type = "l", col = "blue",
     main = "Weekly Seasonal Component (Period = 7)",
     ylab = "Effect", xlab = "Time")
plot(mstl_seasonal_y, type = "l", col = "red",
     main = "Annual Seasonal Component (Period = 365.25)",
     ylab = "Effect", xlab = "Time")
par(mfrow = c(1, 1))

# 5.	Model dependency ---------------------------------------------------
release_version <- as.Date(c("2022-11-30", "2023-03-14", "2024-05-13", "2024-07-18", 
                             "2024-09-12", "2024-12-05", "2025-01-31", "2025-02-27", 
                             "2025-03-15", "2025-03-16", "2025-06-10", "2025-08-07", 
                             "2025-11-12"))

released_dates <- as.numeric(data$Day %in% release_version)

#Create dummies for days
day_week <- weekdays(data$Day) 
days_df <- data.frame(Day = day_week)
levels_in_order <- c("Sunday", "Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday")
days_df$Day <- factor(days_df$Day, levels = levels_in_order) 
week_d <- dummy_cols(days_df, select_columns = "Day", remove_first_dummy = TRUE)
week_d <- week_d[, -1] 

#Put altogether the x variables to use
x_var <- cbind(as.matrix(week_d), as.numeric(released_dates))

ts_data_tbl <- data %>%
  mutate(Day = as.Date(Day)) %>%
  bind_cols(x_var) %>% 
  as_tsibble(index = Day) %>%
  select(Day, interest = ChatGPT_interest,everything()) 

#Testing the release effect on time series
model_with_release <- lm(interest ~ released_dates + Day_Monday + Day_Tuesday + Day_Wednesday + 
                           Day_Thursday + Day_Friday + Day_Saturday, data = ts_data_tbl)
summary(model_with_release)
model_without_release <- lm(interest ~ Day_Monday + Day_Tuesday + Day_Wednesday + 
                              Day_Thursday + Day_Friday + Day_Saturday, data = ts_data_tbl)
summary(model_without_release)
anova(model_without_release, model_with_release)
#The p-value indicates that the released date does not improve the model, therefore, the variable will be discarded

# 7. Forecasting ----------------------------------------------------------
#Annual and weekly seasonality
ts_data <- msts(data$ChatGPT_interest, seasonal.periods = c(7, 365.25))
dummy_vars <- ts_data_tbl %>% select(Day_Monday, Day_Tuesday, Day_Wednesday, Day_Thursday, Day_Friday, Day_Saturday)
#Determine train/test split 
test_periods <- 21  #3 weeks
train_size <- length(ts_data) - test_periods
train <- window(ts_data, end = time(ts_data)[train_size])
test <- window(ts_data, start = time(ts_data)[train_size + 1])

#Split dummy variables
train_dummies <- as.matrix(dummy_vars[1:train_size, -7])
test_dummies <- as.matrix(dummy_vars[(train_size + 1):(train_size + test_periods), -7])

length(train)
length(test)

#Weekly seasonality sets
ts_data_f <- ts(data$ChatGPT_interest, frequency = 7)
train_ts <- window(ts_data_f, end = time(ts_data_f)[train_size])
test_ts <- window(ts_data_f, start = time(ts_data_f)[train_size + 1])


#7.1 Methods -----------------------------------------------------------------
h <- test_periods  

#a. Seasonal Naive Forecast
fit_snaive <- snaive(train, h = h)
forecast_snaive <- forecast(fit_snaive, h = h)
summary(fit_snaive)

#b. Random Walk with Drift
fit_rwdrift <- rwf(train, h = h, drift = TRUE)
forecast_rwdrift <- forecast(fit_rwdrift, h = h)
summary(fit_rwdrift)

#c. Forecasting using Decomposition
fit_stlf <- stlf(train, h = h, method = "ets", s.window = "periodic")
forecast_stlf <- forecast(fit_stlf, h = h)
summary(fit_stlf)

#d. Exponential Smoothing
fit_ets <- ets(train_ts)
forecast_ets <- forecast(fit_ets, h = h)
summary(fit_ets)

#e. ARIMA/SARIMA Model
fit_arima <- auto.arima(train_ts, seasonal = TRUE, stepwise = FALSE, 
                        approximation = FALSE)
forecast_arima <- forecast(fit_arima, h = h)
print(fit_arima)
summary(fit_arima)

#f. ARIMAX with Weekday Dummy Variables
fit_arimax <- auto.arima(train_ts, xreg = train_dummies, seasonal = TRUE,
                         stepwise = FALSE, approximation = FALSE)
forecast_arimax <- forecast(fit_arimax, h = h, xreg = test_dummies)
print(fit_arimax)
summary(fit_arimax)

#g. ARIMA with Fourier Terms
K_week <- 3
fourier_train <- fourier(train_ts, K = K_week)
fourier_test <- fourier(train_ts, K = K_week, h = h)

fit_fourier <- auto.arima(train_ts, xreg = fourier_train, seasonal = FALSE,
                          stepwise = FALSE, approximation = FALSE)
forecast_fourier <- forecast(fit_fourier, h = h, xreg = fourier_test)
print(fit_fourier)
summary(fit_fourier)

#h. Dynamic Regression with trend
time_trend <- 1:length(train_ts)
time_trend_test <- (length(train_ts) + 1):(length(train_ts) + h)
fit_dynreg <- auto.arima(train_ts, xreg = cbind(train_dummies, time_trend), 
                         seasonal = TRUE, stepwise = FALSE, 
                         approximation = FALSE)
forecast_dynreg <- forecast(fit_dynreg, h = h, 
                            xreg = cbind(test_dummies, time_trend_test))
print(fit_dynreg)
summary(fit_dynreg)


# RMSE --------------------------------------------------------------------
calculate_rmse <- function(actual, predicted) {
  actual_vec <- as.numeric(actual)
  predicted_vec <- as.numeric(predicted)
  sqrt(mean((actual_vec - predicted_vec)^2))
}

#Get RSME of each method
rmse_snaive <- calculate_rmse(test, forecast_snaive$mean)
rmse_rwdrift <- calculate_rmse(test, forecast_rwdrift$mean)
rmse_stlf <- calculate_rmse(test, forecast_stlf$mean)
rmse_ets <- calculate_rmse(test_ts, forecast_ets$mean)
rmse_arima <- calculate_rmse(test, forecast_arima$mean)
rmse_arimax <- calculate_rmse(test_ts, forecast_arimax$mean)
rmse_fourier <- calculate_rmse(test_ts, forecast_fourier$mean)
rmse_dynreg <- calculate_rmse(test_ts, forecast_dynreg$mean)

#RMSE dataframe
results <- data.frame(
  Method = c("Seasonal Naive", "RW with Drift", "STL Decomposition", 
             "Exponential Smoothing", "ARIMA", "ARIMAX (Weekday Dummies)",
             "ARIMA + Fourier", 
             "Dynamic Regression"),
  RMSE = c(rmse_snaive, rmse_rwdrift, rmse_stlf, rmse_ets, rmse_arima,
           rmse_arimax, rmse_fourier, rmse_dynreg),
  Uses_Dummies = c("No", "No", "No", "No", "No", "Yes", "No", "Yes"),
  Uses_Fourier = c("No", "No", "No", "No", "No", "No", "Yes", "No")
)

results <- results %>% arrange(RMSE)
results


# LJUNG BOX TEST ----------------------------------------------------------
ljungbox_pvalue <- function(fit_model) {
  check <- checkresiduals(fit_model, plot = FALSE)
  #pvalue test
  return(check$p.value[1]) 
}

#Get p value of each model
p_snaive <- ljungbox_pvalue(fit_snaive)
p_rwdrift <- ljungbox_pvalue(fit_rwdrift)
p_stlf <- ljungbox_pvalue(fit_stlf)
p_ets <- ljungbox_pvalue(fit_ets)
p_arima <- ljungbox_pvalue(fit_arima)
p_arimax <- ljungbox_pvalue(fit_arimax)
p_fourier <- ljungbox_pvalue(fit_fourier)
p_dynreg <- ljungbox_pvalue(fit_dynreg)

#pvalue dataframe
results_p <- data.frame(
  Method = c("Seasonal Naive", "RW with Drift", "STL Decomposition", 
             "Exponential Smoothing", "ARIMA", "ARIMAX (Weekday Dummies)",
             "ARIMA + Fourier", "Dynamic Regression"),
  LjungBox_pvalue = c(p_snaive, p_rwdrift, p_stlf, p_ets, p_arima,
                      p_arimax, p_fourier, p_dynreg), # <-- NUEVA COLUMNA
  Uses_Dummies = c("No", "No", "No", "No", "No", "Yes", "No", "Yes"),
  Uses_Fourier = c("No", "No", "No", "No", "No", "No", "Yes", "No")
)

results_p <- results_p %>% arrange(LjungBox_pvalue)
results_p

# Forecasting methods plot ------------------------------------------------
par(mfrow = c(3, 3), mar = c(4, 4, 2, 1))

plot(forecast_snaive, main = "Seasonal Naive", 
     xlab = "Time", ylab = "Interest")
lines(test, col = "red", lwd = 2)

plot(forecast_rwdrift, main = "RW with Drift", 
     xlab = "Time", ylab = "Interest")
lines(test, col = "red", lwd = 2)

plot(forecast_stlf, main = "STL + ETS", 
     xlab = "Time", ylab = "Interest")
lines(test, col = "red", lwd = 2)

plot(forecast_ets, main = "ETS", 
     xlab = "Time", ylab = "Interest")
lines(test_ts, col = "red", lwd = 2)

plot(forecast_arima, main = "ARIMA", 
     xlab = "Time", ylab = "Interest")
lines(test_ts, col = "red", lwd = 2)

plot(forecast_arimax, main = "ARIMAX (Dummies)", 
     xlab = "Time", ylab = "Interest")
lines(test_ts, col = "red", lwd = 2)

plot(forecast_fourier, main = "ARIMA + Fourier", 
     xlab = "Time", ylab = "Interest")
lines(test_ts, col = "red", lwd = 2)


plot(forecast_dynreg, main = "Dynamic Regression", 
     xlab = "Time", ylab = "Interest")
lines(test_ts, col = "red", lwd = 2)

par(mfrow = c(1, 1))

#Comparison 5 methods
#Plot annual
top_methods_a <- results$Method[c(6,7,8)]
plot(test, type = "l", lwd = 3, col = "black", 
     main = "Annual scale forecasts vs Actual", 
     xlab = "Time", ylab = "Interest")
colors <- c("red", "blue", "green", "orange")
forecasts_list <- list(forecast_snaive, forecast_stlf, forecast_rwdrift)
names(forecasts_list) <- results$Method[c(6,7,8)]

for(i in 1:5) {
  method_name <- top_methods_a[i]
  lines(forecasts_list[[method_name]]$mean, col = colors[i], lwd = 2)
}

legend("topleft", 
       legend = c("Actual", top_methods_a),
       col = c("black", colors),
       lty = 1, lwd = 2, cex = 0.7)


#Plot weekly
top_methods_w <- results$Method[c(1,2,3,4,5)]
plot(test_ts, type = "l", lwd = 3, col = "black", 
     main = "Weekly scale forecasts vs Actual", 
     xlab = "Time", ylab = "Interest")

colors <- c("red", "blue", "green", "purple","orange")
forecasts_list <- list( forecast_ets, forecast_arimax,
                        forecast_fourier, forecast_dynreg, forecast_arima)
names(forecasts_list) <- results$Method[c(1,2,3,4,5)]

for(i in 1:5) {
  method_name <- top_methods_w[i]
  lines(forecasts_list[[method_name]]$mean, col = colors[i], lwd = 2)
}

legend("topleft", 
       legend = c("Actual", top_methods_w),
       col = c("black", colors),
       lty = 1, lwd = 2, cex = 0.7)

# 6. Forecast next 3 periods with final one -------------------------------
checkresiduals(forecast_arima)
full_data_ts <- ts(data$ChatGPT_interest, frequency = 7)
h_final <- 21 #three weeks
#Fit with all the data
fit_arima_final <- auto.arima(
  full_data_ts, 
  seasonal = TRUE, 
  stepwise = FALSE, 
  approximation = FALSE
)
print(fit_arima_final)
#Forecast
forecast_final <- forecast(fit_arima_final, h = h_final)
print(forecast_final)

#Plot forecast
plot(
  forecast_final, 
  main = "Forecast ARIMA next 3 periods (weeks)",
  xlab = "Time", 
  ylab = "Interest"
)
summary(forecast_final)
checkresiduals(forecast_final)


# 8. Own forecasting method ---------------------------------------------------------------
#ARIMA Model with weekday residual correction 
#Using the ARIMA which was the best model
fit_arima_main <- auto.arima(train_ts, seasonal = TRUE, stepwise = FALSE, approximation = FALSE)
forecast_arima_main <- forecast::forecast(fit_arima_main, h = h)

#Get the residuals and capture the behaviour of its
arima_fitted <- fitted(fit_arima_main)
arima_residuals <- residuals(fit_arima_main)

train_days <- weekdays(data$Day[1:train_size])
residual_df <- data.frame(
  residual = as.numeric(arima_residuals),
  day = train_days
)

#Effect of each day
day_effects <- residual_df %>%
  group_by(day) %>%
  summarise(
    mean_residual = mean(residual, na.rm = TRUE),
    sd_residual = sd(residual, na.rm = TRUE)
  ) %>%
  arrange(factor(day, levels = c("Sunday", "Monday", "Tuesday", 
                                 "Wednesday", "Thursday", "Friday", "Saturday")))
print(day_effects)

#Correct the forecast 
forecast_own <- numeric(h)
test_days <- weekdays(data$Day[(train_size + 1):(train_size + h)])
damping_factor <- 0.5 #how much trend is preserved

for (i in 1:h) {
  main_forecast <- forecast_arima_main$mean[i]
  current_day <- test_days[i]
  #Effect of the day
  day_correction <- day_effects$mean_residual[day_effects$day == current_day]
  forecast_own[i] <- main_forecast + (damping_factor * day_correction)
}

#Plot test forecast
ylim_range <- range(c(test_ts, forecast_arima_main$mean, forecast_own))

plot(test_ts, 
     type = "l", 
     col = "black", 
     lwd = 3,
     main = "ARIMA and ARIMA + Correction vs Actual",
     xlab = "Time", 
     ylab = "Interest",
     ylim = ylim_range
)
lines(forecast_arima_main$mean, 
      col = "blue", 
      lwd = 2, 
      lty = 2) 
forecast_own_ts <- ts(forecast_own, 
                      start = start(forecast_arima_main$mean), 
                      frequency = frequency(forecast_arima_main$mean))

lines(forecast_own_ts, 
      col = "red", 
      lwd = 2)  

legend("topleft",
       legend = c("Actual", "ARIMA", "ARIMA + Correction"),
       col = c("black", "blue", "red"),
       lty = c(1, 2, 1), 
       lwd = c(3, 2, 2),
       bg = "white")

grid()

# 9. Measure precision of own method --------------------------------------
#Get RMSE
rmse_main <- calculate_rmse(test_ts, forecast_arima_main$mean)
rmse_own <- calculate_rmse(test_ts, forecast_own)
rmse_own

#Get the residuals of the model for Ljung Box test
own_fitted <- numeric(train_size)
train_days_all <- weekdays(data$Day[1:train_size])

for (i in 1:train_size) {
  main_fit <- arima_fitted[i]
  current_day <- train_days_all[i]
  day_correction <- day_effects$mean_residual[day_effects$day == current_day]
  own_fitted[i] <- main_fit + (damping_factor * day_correction)
}

own_residuals <- train_ts - own_fitted
lb_test_main <- Box.test(forecast_arima_main$residuals, 
                        lag = min(20, length(own_residuals)/5), 
                        type = "Ljung-Box")
lb_test_own <- Box.test(own_residuals, 
                        lag = min(20, length(own_residuals)/5), 
                        type = "Ljung-Box")
lb_test_main
lb_test_own

# Forecasting next 3 weeks with own method ------------------------------------------------
#Fit own model with full dataset
full_data_ts <- ts(data$ChatGPT_interest, frequency = 7)
n_full <- length(full_data_ts)
fit_arima_full_final <- auto.arima(full_data_ts, seasonal = TRUE, 
                                   stepwise = FALSE, approximation = FALSE)
forecast_arima_full_final <- forecast::forecast(fit_arima_full_final, h = 21)
print(fit_arima_full_final)

#Get the residuals and weekdays effects from full dataset
arima_fitted_full_final <- fitted(fit_arima_full_final)
arima_residuals_full_final <- residuals(fit_arima_full_final)
full_days <- weekdays(data$Day)
residual_df_full_final <- data.frame(
  residual = as.numeric(arima_residuals_full_final),
  day = full_days
)

#Calculate day effects from full dataset
day_effects_full_final <- residual_df_full_final %>%
  group_by(day) %>%
  summarise(
    mean_residual = mean(residual, na.rm = TRUE),
    sd_residual = sd(residual, na.rm = TRUE)
  ) %>%
  arrange(factor(day, levels = c("Sunday", "Monday", "Tuesday", 
                                 "Wednesday", "Thursday", "Friday", "Saturday")))
print(day_effects_full_final)

#Get the date and weekdays
last_date <- max(data$Day)
future_dates <- seq(last_date + 1, by = "day", length.out = 21)
future_days <- weekdays(future_dates)

#Damping factor (same as before) for the future forecast
forecast_full_final_own <- numeric(21)
damping_factor <- 0.5
for (i in 1:21) {
  base_forecast <- forecast_arima_full_final$mean[i]
  current_day <- future_days[i]
  #Get the correction
  day_correction <- day_effects_full_final$mean_residual[day_effects_full_final$day == current_day]
  #Apply correction (damping)
  forecast_full_final_own[i] <- base_forecast + (damping_factor * day_correction)
}

#Summary
forecast_summary <- data.frame(
  Period = 1:21,
  Date = as.character(future_dates),
  Weekday = future_days,
  Base_ARIMA = round(as.numeric(forecast_arima_full_final$mean), 2),
  Day_Correction = round(damping_factor * 
                           sapply(future_days, function(d) {
                             day_effects_full_final$mean_residual[day_effects_full_final$day == d]
                           }), 2),
  Full_Final_Forecast = round(forecast_full_final_own, 2),
  Lower_80 = round(as.numeric(forecast_arima_full_final$lower[, 1]), 2),
  Upper_80 = round(as.numeric(forecast_arima_full_final$upper[, 1]), 2),
  Lower_95 = round(as.numeric(forecast_arima_full_final$lower[, 2]), 2),
  Upper_95 = round(as.numeric(forecast_arima_full_final$upper[, 2]), 2)
)
print(forecast_summary)

#Plot forecast
n_full <- length(full_data_ts)
plot_window <- 500
start_idx <- n_full - plot_window + 1

plot(start_idx:n_full, tail(full_data_ts, plot_window), 
     type = "l", lwd = 2, col = "black",
     xlim = c(start_idx, n_full + 21),
     ylim = range(c(tail(full_data_ts, plot_window), 
                    forecast_full_final_own, 
                    forecast_arima_full_final$upper[,2])),
     main = "ARIMA + Correction forecast next 3 periods (21 days)",
     xlab = "Time", 
     ylab = "ChatGPT Search Interest")
#95% CI
polygon(x = c((n_full + 1):(n_full + 21),
              rev((n_full + 1):(n_full + 21))),
        y = c(forecast_arima_full_final$lower[,2], 
              rev(forecast_arima_full_final$upper[,2])),
        col = rgb(0, 0, 1, 0.15), border = NA)
#80% CI
polygon(x = c((n_full + 1):(n_full + 21),
              rev((n_full + 1):(n_full + 21))),
        y = c(forecast_arima_full_final$lower[,1], 
              rev(forecast_arima_full_final$upper[,1])),
        col = rgb(0, 0, 1, 0.25), border = NA)
#Own method forecast
lines(x = (n_full + 1):(n_full + 21),
      y = forecast_full_final_own, 
      col = "red", lwd = 3)

points(x = (n_full + 1):(n_full + 21),
       y = forecast_full_final_own, 
       col = "red", pch = 19, cex = 1.8)

legend("topleft",
       legend = c("Historical Data", 
                  "ARIMA + Correction)", 
                  "80% Confidence Interval",
                  "95% Confidence Interval"),
       col = c("black", "red", 
               rgb(0, 0, 1, 0.25), rgb(0, 0, 1, 0.15)),
       lty = c(1, 1, 1, 1), 
       lwd = c(2, 3, 10, 10),
       pch = c(NA, 19, 15, 15),
       bg = "white")

grid()

#Get metrics
fit_arima_full_final$aic
fit_arima_full_final$bic
Box.test(fit_arima_full_final$residuals, 
         lag = min(20, length(fit_arima_full_final)/5), 
         type = "Ljung-Box")
training_accuracy <- accuracy(fit_arima_full_final)
rmse_training <- training_accuracy[1, "RMSE"]
rmse_training


# 10.	Estimate forecasting accuracy dependency on horizon for the own method proposed.--------
#Rolling origin, 20 forecast for horizon 1, 3, 7, 14, 21, 28
horizons <- c(1, 3, 7, 14, 21, 28)  #Days 
damping_factor <- 0.5

#Using the full time series
full_data_ts <- ts(data$ChatGPT_interest, frequency = 7)
n_full <- length(full_data_ts)

min_train_size <- 365  #Training size of a year
max_horizon <- max(horizons)
n_origins <- 20 

#Get the origin points for each rolling
origin_points <- seq(from = min_train_size, 
                     to = n_full - max_horizon, 
                     length.out = n_origins)
origin_points <- round(origin_points)

#Dataframe to save the results
rmse_by_horizon <- data.frame(
  Horizon = integer(),
  RMSE = numeric(),
  Origin = integer()
)

#Loop to forecast with n horizon in each origin
for (origin_idx in 1:length(origin_points)) {
  origin <- origin_points[origin_idx]
  train_data <- window(full_data_ts, end = time(full_data_ts)[origin]) #train set
  
  #Fit and correct the own method
  fit_arima <- auto.arima(train_data, seasonal = TRUE, 
                          stepwise = FALSE, approximation = FALSE)
  arima_residuals <- residuals(fit_arima)
  train_days <- weekdays(data$Day[1:origin])
  
  #Save the day effects
  day_effects <- data.frame(
    residual = as.numeric(arima_residuals), 
    day = train_days
  ) %>%
    group_by(day) %>%
    summarise(mean_residual = mean(residual, na.rm = TRUE))
  
  #Forecast fo the respective horizon
  for (h in horizons) {
    #To see if there is enough data
    if (origin + h > n_full) next
    
    #Actual values
    actual_values <- full_data_ts[(origin + 1):(origin + h)]
    #Forecast
    forecast_base <- forecast(fit_arima, h = h)
    
    #Forecast + correction
    forecast_days <- weekdays(data$Day[(origin + 1):(origin + h)])
    forecast_corrected <- numeric(h)
    for (i in 1:h) {
      base_pred <- forecast_base$mean[i]
      current_day <- forecast_days[i]
      day_correction <- day_effects$mean_residual[day_effects$day == current_day]
      forecast_corrected[i] <- base_pred + (damping_factor * day_correction)
    }
    
    #Get RMSE of the forecast
    rmse_val <- sqrt(mean((actual_values - forecast_corrected)^2))
    rmse_by_horizon <- rbind(rmse_by_horizon, 
                             data.frame(Horizon = h, 
                                        RMSE = rmse_val,
                                        Origin = origin))
  }
}

#Compute the mean RMSE of each horizon
rmse_summary <- rmse_by_horizon %>%
  group_by(Horizon) %>%
  summarise(
    Mean_RMSE = mean(RMSE),
    SD_RMSE = sd(RMSE),
    Min_RMSE = min(RMSE),
    Max_RMSE = max(RMSE),
    N = n()
  ) %>%
  arrange(Horizon)
print(rmse_summary)

#Plot mean RMSE vs horizon
ggplot(rmse_summary, aes(x = Horizon, y = Mean_RMSE)) +
  geom_line(color = "darkred", linewidth = 1.2) +
  geom_point(color = "red", size = 4) +
  geom_errorbar(aes(ymin = Mean_RMSE - SD_RMSE, 
                    ymax = Mean_RMSE + SD_RMSE),
                width = 0.5, color = "gray50", linewidth = 0.8) +
  labs(
    title = "Mean RMSE forecast vs Horizon",
    subtitle = "ARIMA + Correction method ",
    x = "Forecast Horizon (days)",
    y = "Mean RMSE"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(size = 12)
  ) +
  scale_x_continuous(breaks = horizons)

#Linear regression for RMSE and horizon to see its dependency
lm_rmse <- lm(RMSE ~ Horizon, data = rmse_by_horizon)
print(summary(lm_rmse))

