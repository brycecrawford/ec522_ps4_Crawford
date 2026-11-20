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

#Importing quarterly data set and slicing off the first two rows
qd_raw <- read_csv("data/fredqd_2026-09.csv", show_col_types = FALSE) |> 
  slice(-(1:2))
        
#Building the growth series
growth <- function(x) 100 * difference(log(x))

qd <- qd_raw |> 
  mutate( date = as.Date(sasdate, "%m/%d/%Y")) |> 
  transmute( date, 
             gdp    = growth(GDPC1),
             cons   = growth(PCECC96),
             iv     = growth(GPDIC1),
             emp    = growth(PAYEMS),
             infl   = growth(GDPCTPI),
             unrate = difference(UNRATE),
             ffr    = difference(FEDFUNDS),
             gs10   = difference(GS10)) |> 
  filter(date >= as.Date("1960-03-01"), date <= as.Date("2019-12-01"))

#Creating matrix for BVAR
y <- qd |> 
  select(-date) |> 
  as.matrix()

nrow(y)
head(y)

#Setting up TS-CV
N <- nrow(y)
origins <- 120:(N - 1)
actual <- y[121:N, ]

# First and last forecast dates
range(qd$date[121:N])

# Number of forecasts per variable
length(origins) * ncol(y)

#Estimating VAR
var_ols_fc <- function(z,p) {
  e <-  embed(z,p +1)
  M <- ncol(z)
  B <- coef(lm(e[, 1:M] ~ e[, -(1:M)]))
  x <- c(1, t(z[nrow(z):(nrow(z) - p + 1), ]))
  x %*% B
}

fc_var1 <-  NULL
fc_var4 <- NULL
for (T in origins) {
  fc_var1 <- rbind(fc_var1, var_ols_fc(y[1:T, ],1))
  fc_var4 <- rbind(fc_var4, var_ols_fc(y[1:T, ],4))
}

#Computing RMSE of each model for each variable
rmse_by_var <- function(f) sqrt(colMeans((actual - f)^2))

rmse_by_var(fc_var1)
rmse_by_var(fc_var4)

# Estimating BVAR(4) with a Minnesota Prior
bvar_cv <- function(lambda) {
  pr <- bv_priors(hyper = "lambda",
                  mn = bv_mn(lambda = bv_lambda(mode = lambda, sd = 0.4,
                                                min = 0.999 * lambda,
                                                max = 1.001 * lambda),
                             b = 0))
  fc <- NULL
for (T in origins) {
  fit_T <- bvar(y[1:T, ], lags = 4, priors = pr,
                n_draw = 2000, n_burn = 1000, verbose = FALSE)
  x <- c(1, t(y[T:(T - 3), ]))
  fc <- rbind(fc, x %*% coef(fit_T))
}
  fc
}

set.seed(422)
fc_bvar_03 <- bvar_cv(0.3)
saveRDS(fc_bvar_03, "results/fc_bvar_03.rds")

# Estimating BVAR(4) with tighter lambda
set.seed(422)
fc_bvar_05 <- bvar_cv(0.05)
saveRDS(fc_bvar_05, "results/fc_bvar_05.rds")

#VAR and BVAR Results
var_results <- tibble(
  Variable = colnames(y),
  `VAR(1)` = rmse_by_var(fc_var1),
  `VAR(4)` = rmse_by_var(fc_var4),
  `BVAR(4) lambda = 0.3` = rmse_by_var(fc_bvar_03),
  `BVAR(4) lambda = 0.05` = rmse_by_var(fc_bvar_05))

saveRDS(var_results, "results/var_results.rds")

var_results |> knitr::kable(digits = 3)

e_var4 <- actual[, "gdp"] - fc_var4[, 1]
e_bvar <- actual[, "gdp"] - fc_bvar_03[, 1]

forecast::dm.test(e_var4, e_bvar, h = 1, power = 2)

