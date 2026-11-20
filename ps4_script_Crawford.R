pacman::p_load(tidyverse, fpp3, BVAR, feasts)

#Importing data set from file folder
liquor_raw <- read_csv("Data/MRTSSM4453USN.csv",
                       na = ".",
                       col_types = cols(
                         observation_date = col_date(format = "%Y-%m-%d"),
                         MRTSSM4453USN = col_double()))

# Creating time series object indexed monthly and renaming columns
liquor <- liquor_raw |> 
  transmute(
    Month = yearmonth(observation_date),
    sales = MRTSSM4453USN
  ) |> 
  as_tsibble( index = Month) |> 
  filter(Month <= yearmonth("2019 Dec") )

# checking for gaps in time series 
has_gaps(liquor)

#graphing sales with autoplot
autoplot(liquor, sales) +
  labs(
    title = "Monthly Retail Sales of liquor in the United States (Millions of Dollars)",
    x = "Month",
    y = "Sales (Millions of Dollars)")

#Graphing seasonal monthly sales by year with gg_season
gg_season(liquor, log(sales)) +
  labs(
    title = "Monthly Retail Sales of liquor in the United States (Log-Scale)",
    x = "Month",
    y = "Sales (Log-Scale)"
  )

#Fitting an ARIMA model and ETS model to the log of sales
liquor_fit <- liquor |> 
  model(arima = ARIMA(log(sales)),
        ets = ETS(log(sales)))

#Assessing specification of each model
liquor_fit |> select(arima) |>  report()
liquor_fit |> select(ets) |> report()

# Checking residuals of ARIMA model
liquor_fit |> select(arima) |> gg_tsresiduals()

#Ljung-Box test of residuals of ARIMA model
liquor_fit |> 
  select(arima) |> 
  augment() |> 
  features(.innov, ljung_box, lag = 24, dof = 6)

#Horse Race
liquor_cv <- liquor |> 
  stretch_tsibble(.init = 120, .step = 1) |> 
  filter(.id < max(.id))

liquor_cv_fit <- liquor_cv |> 
  model(snaive = SNAIVE(sales),
        ets = ETS(log(sales)),
        arima = ARIMA(log(sales)))

saveRDS(liquor_cv_fit, "results/liquor_cv_fit.rds")

# Producing one-month ahead forecasts
liquor_cv_fc <- liquor_cv_fit |> 
  forecast(h = 1)

#Calculating accuracy of forecasts
liquor_cv_fc |> 
  accuracy(liquor, measure = list(RMSE = RMSE, MAE = MAE)) |> 
  knitr::kable(digits = 2)

#Producing first and last forecast
range(liquor_cv_fc$Month)

#creating tibble of forecast errors for each model
liquor_errs <- liquor_cv_fc |> 
  as_tibble() |> 
  left_join(as.tibble(liquor) |> select(Month, actual = sales),
            by = "Month") |> 
  mutate(e = actual - .mean) |> 
  select(Month, .model, e) |> 
  pivot_wider(names_from = .model, values_from = e) 

#Diebold-Mariano test comparing ARIMA to ETS
forecast::dm.test(liquor_errs$arima, liquor_errs$ets, h =1, power = 2)

#Diebold-Mariano test comparing ARIMA to SNAIVE
forecast::dm.test(liquor_errs$arima, liquor_errs$snaive, h = 1, power = 2)

# Model Combination of ARIMA and ETS
liquor_combo_fc <- liquor_cv_fit |> 
  mutate(combo = (arima + ets) / 2 ) |> 
  forecast(h = 1)

#Testing accuracy of model combinations to their individual models
liquor_combo_fc |> 
  accuracy(liquor, measures = list(RMSE = RMSE, MAE = MAE)) |> 
  knitr::kable(digits = 2)

#Computing correlation of errors between ARIMA and ETS
cor(liquor_errs$arima, liquor_errs$ets)
           