# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Install the required packages if necessary:
list_of_packages = c("ggplot2", "zoo", "dplyr", "lubridate", "forecast", "prophet",
                     "xts", "smooth", "tsintermittent", "tscount")
new_packages = list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)){
    cat("INSTALLING THE FOLLOWING PACKAGES:\n", 
        paste(list_of_packages, collapse = "\n"), "\n", sep = "")
    install.packages(new_packages)
}
    

library(ggplot2)

source("fire_season_analysis.R")

cat("FIRE SEASONS SEVERITY\n---------------------\n")
cat("Example for the time series in figure 1 (India) in the paper.\n\n")
devAskNewPage(ask = TRUE)
par(ask=TRUE)

# ---------------------------------------------------------------------------------
cat("1) Load and plot time series\n")
# Loading data
ts_example = read.csv(file = "data_exemple.csv", colClasses=c("Date", "numeric"))
# Plotting 
show(ggplot(ts_example, aes(x=time, y=fc)) + 
    geom_line() + theme_bw() + 
    ggtitle("Daily fire counts", "raw data"))

# ---------------------------------------------------------------------------------
cat("2) Plotting the periods with or without fire\n")
lp = largest_periods_without_fire(ts_example$fc, num_periods = 15, consider_more_than_days = 2, sma_k = 7)
fs_nofire = do.call(rbind,lapply(1:nrow(lp), function(i) data.frame(time=ts_example$time[lp[i,"start"]:lp[i,"end"]], fc=0)))
show(ggplot(ts_example, aes(x=time, y=fc)) +
    geom_line(size=0.1, color="#d25865") +
    geom_point(data = fs_nofire, mapping = aes(x=time, y=fc), size=0.1, shape=46, color="#86aad2") +
    theme_bw() +
    ggtitle("Daily fire counts", "Red: periods of fire. Blue: periods without fire"))

# ---------------------------------------------------------------------------------
cat("3) Obtain the monthly-accumulated fire counts time series\n")
# Defining the length of the fire season in months
FIRE_SEASON_TIME_WINDOW_MONTHS = 7
# Get the month with the highest occurance of fire
maxmonth = month_with_most_frequent_peaks_climatology(ts_example)
# Extract the fire seasons centered in the month with the highest occurance of fire
fire_seasons = estimate_fire_season_max_month_centered(ts = ts_example, 
                                                month_break = maxmonth, 
                                                month_window_length = FIRE_SEASON_TIME_WINDOW_MONTHS)
# Remove the first and last season to avoid measuring broken seasons
fire_seasons = fire_seasons[c(-1, -length(fire_seasons))]
# Group by month
fire_seasons = do.call(rbind, lapply(fire_seasons, group_by_month))
# Plotting
show(ggplot(fire_seasons, aes(x=time, y=fc)) + 
    geom_line() + geom_point() + theme_bw() + 
    ggtitle("Fire seasons", "Monthly fire counts"))

# ---------------------------------------------------------------------------------
cat("4) Training the forecasting methods (may take while to run)\n")

ADD_CONST = 1
BOXCOX_LAMBDA = "auto"

# Select the first 10 seasons as the training set
ts_train = head(fire_seasons, 10 * FIRE_SEASON_TIME_WINDOW_MONTHS)
# Adding a small constant value to avoid zeros
ts_train$fc = ts_train$fc + ADD_CONST
# Transforming it into a ts object
ts_train_ts = ts(ts_train$fc, frequency=FIRE_SEASON_TIME_WINDOW_MONTHS) 
# Select the lasting 3 seasons as the test set
ts_test = tail(fire_seasons, 3 * FIRE_SEASON_TIME_WINDOW_MONTHS)
# Adding a small constant value to avoid zeros
ts_test$fc = ts_test$fc + ADD_CONST
# Transforming it into a ts object
ts_test_ts = ts(ts_test$fc, frequency=FIRE_SEASON_TIME_WINDOW_MONTHS, start = end(ts_train_ts)[1] + 1)

# List to keep the forecasting results
models_predictions = list()

# 1. Seasonal naive (snaive):
prediction_snaive = snaive(ts_train_ts, h= 3 * FIRE_SEASON_TIME_WINDOW_MONTHS, lambda = BOXCOX_LAMBDA)
prediction_snaive = pmax(prediction_snaive$mean - ADD_CONST, 0)
prediction_snaive = data.frame(time=ts_test$time, fc=as.numeric(prediction_snaive))
models_predictions[["snaive"]] = prediction_snaive

# 2. ARIMA:
fit_arima = auto.arima(ts_train_ts, lambda = BOXCOX_LAMBDA)
prediction_arima = forecast(fit_arima, h = 3 * FIRE_SEASON_TIME_WINDOW_MONTHS)
prediction_arima = data.frame(time=ts_test$time, 
                              fc=round(pmax(as.numeric(prediction_arima$mean - ADD_CONST),0)))
models_predictions[["ARIMA"]] = prediction_arima

# 3. ETS:
fit_ets = ets(ts_train_ts, lambda = BOXCOX_LAMBDA)
prediction_ets = forecast(fit_ets, 3 * FIRE_SEASON_TIME_WINDOW_MONTHS)
prediction_ets = data.frame(time=ts_test$time, 
                            fc=round(pmax(as.numeric(prediction_ets$mean - ADD_CONST),0)))
models_predictions[["ETS"]] = prediction_ets

# 4. STLF
prediction_stlf = stlf(ts_train_ts, h = 3 * FIRE_SEASON_TIME_WINDOW_MONTHS, lambda = BOXCOX_LAMBDA)
prediction_stlf = data.frame(time=ts_test$time, 
                            fc=round(pmax(as.numeric(prediction_stlf$mean - ADD_CONST),0)))
models_predictions[["STLF"]] = prediction_stlf

# 5. TBATS
fit_tbats = tbats(ts_train_ts)
prediction_tbats = forecast(fit_tbats, 3 * FIRE_SEASON_TIME_WINDOW_MONTHS)
prediction_tbats = data.frame(time=ts_test$time, 
                             fc=round(pmax(as.numeric(prediction_tbats$mean - ADD_CONST),0)))
models_predictions[["TBATS"]] = prediction_tbats

# 6. TSGLM
fit_tsglm = tsglm(ts_train_ts, model = list(past_obs = c(1, FIRE_SEASON_TIME_WINDOW_MONTHS), 
                                      past_mean = FIRE_SEASON_TIME_WINDOW_MONTHS))
prediction_tsglm = predict(fit_tsglm, n.ahead =  3 * FIRE_SEASON_TIME_WINDOW_MONTHS)
prediction_tsglm = data.frame(time=ts_test$time, 
                              fc=round(pmax(as.numeric(prediction_tsglm$median - ADD_CONST),0)))
models_predictions[["TSGLM"]] = prediction_tsglm

# 7. NEURAL NETWORKS
fit_ann = nnetar(ts_train_ts, lambda = BOXCOX_LAMBDA)
prediction_ann = forecast(fit_ann, 3 * FIRE_SEASON_TIME_WINDOW_MONTHS)
prediction_ann = data.frame(time=ts_test$time, 
                            fc=round(pmax(as.numeric(prediction_ann$mean - ADD_CONST),0)))
models_predictions[["neural_net"]] = prediction_ann

# 8. Prophet:
names(ts_train) = names(ts_test) = c("ds", "y")
prediction_prophet = forecast_facebook(ts_train = ts_train, ts_test = ts_test,
                                       add_const = 0, boxcox_lambda = "auto",
                                       return_accuracy = F)
prediction_prophet = data.frame(time=ts_test$ds, 
                                fc=round(pmax(prediction_prophet - ADD_CONST,0)))
models_predictions[["prophet"]] = prediction_prophet

# Constructing a table with the results
results = do.call(rbind, lapply(1:length(models_predictions), function(id){
    model = models_predictions[[id]]
    fire_seasons$data = "raw";
    model$data = "forecast" 
    result = rbind(fire_seasons, model)
    result$method = names(models_predictions)[id]
    result
}))
results$data = factor(results$data, levels = sort(unique(results$data), decreasing = T), ordered = T)
results$method = factor(results$method, levels = names(models_predictions), ordered = T)

# plotting
show(ggplot(results, aes(x=time, y=fc, color=data)) +
    geom_line(size=0.5, alpha=0.65) +
    geom_point(size=0.8,  alpha=0.65) +
    scale_color_manual(values = c("gray10", "#E41A1C")) +
    facet_grid(method~.) +
    ylab("FC") + 
    theme(legend.position="top", axis.title.x = element_blank()))

devAskNewPage(ask = FALSE)
par(ask=FALSE)