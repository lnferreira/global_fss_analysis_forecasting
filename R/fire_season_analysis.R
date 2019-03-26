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


# =================================================================================
# Method to find the months with the largest occurance of fire by year
# This method groups by months, add and resturn the max
# =================================================================================
max_fire_count_month <- function(ts, min_firedays_in_year) {
    ts$year = year(ts$time)
    ts$month = month(ts$time)
    fc_sum_year = ts %>% group_by(year) %>% tally(fc)
    if (any(fc_sum_year$n < min_firedays_in_year))
        return(NaN)
    fc_sum_month = ts %>% group_by(year, month) %>% tally(fc)
    years =  sort(unique(fc_sum_month$year))
    max_months = sapply(years, function(y) which.max(fc_sum_month[fc_sum_month$year == y,]$n))
    names(max_months) = years
    max_months
}

# =================================================================================
# Method to find the most frquent month with peak of fire counts
# This method groups by months, add and returns the max
# =================================================================================
month_with_most_frequent_peaks_climatology <- function(ts) {
    ts$year = year(ts$time)
    ts$month = month(ts$time)
    # ts$has_fire = ifelse(ts$fc > 0, 1, 0)
    # fc_sum_month = ts %>% group_by(month) %>% tally(has_fire)
    # fc_sum_month[which.max(fc_sum_month$n),]$month
    months_peak = ts %>% group_by(year, month) %>% tally(fc) %>% 
        filter(n == max(n)) %>% pull(month)
    as.numeric(names(sort(table(months_peak), T))[1])
}

# =================================================================================
# Returns the periods without fire.
#   ts: time series
#   num_periods: desired number of periods
#   consider_more_than_days: min number of days without fire to consider a period 
#                            without fire
#   sort_result: sort the results according to starting time
#   sma_k: rolling means smooth applied to the time series
#   warn: show warnings
# =================================================================================
largest_periods_without_fire <- function(ts, num_periods=1, consider_more_than_days=1, 
                                         sort_result=TRUE, sma_k=NA, warn=TRUE){
    ts = as.numeric(ts)
    if (!is.na(sma_k))
        ts = floor(rollmean(ts, sma_k))
    if (length(ts) == 0)
        return(0)
    current_zeros = 0
    all_zeros_windows_lengths = c()
    all_periods = list()
    current_period = c()
    consecutive_days_with_fire = 0
    for (i in 1:length(ts)){
        if (ts[i] == 0) {
            current_zeros = current_zeros + 1
            current_period = c(current_period, i)
            consecutive_days_with_fire = 0
        }
        if (ts[i] != 0) {
            consecutive_days_with_fire = consecutive_days_with_fire + 1
            if (consecutive_days_with_fire >= consider_more_than_days & current_zeros > 0){
                all_zeros_windows_lengths = c(all_zeros_windows_lengths, current_zeros)
                all_periods[[length(all_periods)+1]] = current_period
                current_zeros = 0
                current_period = c()
            }
        }
    }
    if (current_zeros > 0){ # add the last period
        all_zeros_windows_lengths = c(all_zeros_windows_lengths, current_zeros)
        all_periods[[length(all_periods)+1]] = current_period
    } 
    if (is.infinite(num_periods))
        num_periods = length(num_periods)
    result = all_zeros_windows_lengths
    if (length(result) > 1 & sort_result){
        if (num_periods > length(result)){
            if (warn)
                warning(sprintf("Desired number of periods: %d, max periods found: %d\n", num_periods, length(result)))
            num_periods = length(result)
        }
        longest_indices = order(all_zeros_windows_lengths, decreasing = T)[1:num_periods]
        result = result[longest_indices]
        all_periods = all_periods[longest_indices]
    }
    if (length(result) > 0){
        all_periods = as.data.frame(do.call(rbind, lapply(all_periods, range)))
        names(all_periods) = c("start", "end")
        all_periods$length = apply(all_periods, 1, diff) + 1
        result = all_periods
        if (sort_result)
            result = result[order(result$start),]
    }
    result
}

# =================================================================================
# Estimate the fire seasons from a time series of fire counts.
#   ts: time series
#   num_seasons: desired number of seasons
#   return_time: return time indexes or time series values
#   consider_more_than_days: min number of days without fire to consider a period 
#                            without fire
#   min_function_remove: if the number of found seasons is higher than the 
#                        desired number of seasons (num_seasons), this function is 
#                        applied to decreasingly sort the seasons and keep those 
#                        with the highest values.
#   sma_k: rolling mean smooth applied to the time series
#   warn: show warnings
# =================================================================================
estimate_fire_season <- function(ts, num_seasons, return_time=TRUE, consider_more_than_days=1, 
                                 min_function_remove=length, sma_k=NA, warn=TRUE) {
    lpwf = largest_periods_without_fire(ts, num_seasons, consider_more_than_days, sma_k = sma_k, warn = warn)
    if (!is.null(lpwf)){
        for (i in 1:nrow(lpwf)){
            ts[lpwf[i,]$start:lpwf[i,]$end] = NaN
        }
    }
    seasons = list()
    current_season = c()
    in_season = FALSE
    for (i in 1:length(ts)) {
        if (!is.nan(ts[i])){
            in_season = TRUE
            current_season = c(current_season, ifelse(return_time, i, ts[i]))
        } else {
            if (in_season){
                seasons[[length(seasons)+1]] = current_season
                current_season = c()
                in_season= FALSE
            }
        }
    }
    if (length(current_season)) { # Add the last season
        seasons[[length(seasons)+1]] = current_season
    }
    if (length(seasons) > num_seasons) {
        all_indexes = 1:length(seasons)
        keep = order(sapply(seasons, min_function_remove), decreasing = T)[1:num_seasons]
        seasons = seasons[all_indexes %in% keep]
    }
    if (warn & length(seasons) != num_seasons)
        warning(sprintf("Desired number of seasons: %d, returned: %d\n", num_seasons, length(seasons)))
    seasons
}

# =================================================================================
# 
# =================================================================================
estimate_fire_season_max_month_centered  <- function(ts, month_break, month_window_length=7) {
    years = unique(year(ts$time))
    fire_seasons = list()
    for (y in years){
        mb = ISOdate(y, month_break, 1, hour = 0, min = 0, sec = 0)
        start = mb %m-% months((month_window_length-1)/2)
        end = mb %m+% months(((month_window_length-1)/2) + 1) - days(1)
        int = interval(start, end)
        fire_seasons[[length(fire_seasons)+1]] = ts[ts$time %within% int,]
    }
    fire_seasons
}

# =================================================================================
# This function makes all fire seasons the same length by removing day from the
# beginning and end
# =================================================================================
equalize_lengths_fire_seasons <- function(fire_seasons, remove_first_last_season=T){
    if (remove_first_last_season){
        fire_seasons = lapply(fire_seasons, function(x) x[!(x$season %in% range(x$season)),])
        for (i in 1:length(fire_seasons))
            fire_seasons[[i]]$season = fire_seasons[[i]]$season - 1
    }
    fss_lengths = lapply(fire_seasons, function(fs) (fs %>% group_by(season) %>% tally())$n)
    min_lenght = sort(unique(sapply(fss_lengths, unique)))[1]
    seasons = sort(as.numeric(unique(do.call(rbind, lapply(fire_seasons, function(fs) unique(fs$season))))))
    for (i in 1:length(fire_seasons)) {
        fs = fire_seasons[[i]]
        season_length = length(which(fs$season == 1))
        if (season_length != min_lenght) {
            remove_begin = floor((season_length - min_lenght) / 2)    
            remove_end = season_length - min_lenght - remove_begin
            fs = do.call(rbind, lapply(seasons, function(s){
                f = fs[fs$season == s,]
                f[(1+remove_begin):(nrow(f)-remove_end),]
            }))
            sl = (fs %>% group_by(season) %>% tally())$n
            if (any(sl != min_lenght))
                stop("Different length than expected (id=", i, ")\n")
            fire_seasons[[i]] = fs
        } 
    }
    fire_seasons
}

# =================================================================================
# Remove outliers from a vector x
# =================================================================================
remove_outliers <- function(x, return_type=c("remove", "id", "na")) {
    ids = which(x %in% boxplot.stats(x)$out)
    if (return_type[1] == "remove" & length(ids) > 0){
        return(x[-(ids)])
    } else if (return_type[1] == "na" & length(ids) > 0){
        x[ids] = NA
        return(x)
    }
    ids
}