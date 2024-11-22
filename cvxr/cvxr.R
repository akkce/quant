# Kevin Elliott
# 2024-11-21
#
# A simple example of an ETF portfolio optimization via CVXR and publicly 
# available returns data from Yahoo Finance. 
#
# The portfolio is re-balanced once a quarter on the 3rd Friday using the prior
# four quarters of returns data to build the risk model. The objective function
# uses each ETF's historical cumulative returns as predicted alpha and attempts
# to maximize it, subject to a penalty for residual variance essentially
# making this a momentum strategy. Hard constraints are the portfolio needs to
# be fully invested, long-only and have a resulting volatility equal to or less 
# than the benchmark (a standard 60/40 portfolio). Additionally after the first
# replace period a turnover constraint of 20% or less is in place.
#
# The resulting portfolio weights and associated statistics (volatility, 
# active risk, turnover and tracking error to the benchmark) are saved in the
# pre-allocated results matrix where they are subsequently used to compare 
# performance an holdings to the benchmark.

# DISCLAIMER: This is academic example and demonstration of R and CVXR's ability
# to be used for portfolio optimization. This is not a real world portfolio 
# model and should NOT be used for any actual investing. 

# load libraries
library("tidyverse")
library("lubridate")
library("tidyquant")
library("PerformanceAnalytics")
library("yfR")
library("CVXR")

# set date range for portfolio
sdt <- as.Date("2010-01-01")
edt <- Sys.Date()

# get OHLC pricing data via YF
#   we decrement our starting date so we can ensure we have price return
#   calculations starting on the first market day of 2010
etfs <- c("SPY", "MDY", "IJR", "AGG", "BIL", "GLD")
ohlc_df <- yf_get(
  tickers = etfs,
  first_date = sdt - 5,
  last_date = edt
)



# portfolio re-balance dates:
#   get month and weekdays and then for every quarter
#   get the 3rd Friday and following Monday
rebal_dates_v <- ohlc_df %>%
  select(dt = ref_date) %>%
  filter(dt >= sdt) %>% 
  mutate(
    yyyy = format(dt, "%Y"),
    month = format(dt, "%b"),
    dow = format(dt, "%a")) %>%
  group_by(yyyy, month, dow) %>%
  mutate(n = 1:n()) %>%
  filter(str_detect(month, "Mar|Jun|Sep|Dec") & str_detect(dow, "Mon|Fri") & n == 3) %>%
  pull(dt)

# handle the special cases of the first and last re-balance periods
rebal_dates_v <- c(sdt, rebal_dates_v, max(ohlc_df$ref_date))

# off-set dates to get rebalance date starting and ending pairs
rebal_dates <- data.frame(
  rebal.sdt = rebal_dates_v[seq(1, length(rebal_dates_v), 2)],
  rebal.edt = rebal_dates_v[seq(2, length(rebal_dates_v), 2)]
)

# calculate days in each rebalance period
rebal_dates$rebal.days <- rebal_dates$rebal.edt - rebal_dates$rebal.sdt

# create corresponding sample period dates for each rebalance date pair
#   out-of-sample period is 1 year (4 quarters)
rebal_dates$sample.sdt <- lag(rebal_dates$rebal.sdt, 4)
rebal_dates$sample.edt <- lag(rebal_dates$rebal.edt, 1)
rebal_dates$sample.days <- rebal_dates$sample.edt - rebal_dates$sample.sdt

# remove first re-balance date w/o a sample
rebal_dates <- tail(rebal_dates, -4)

# indexer variable
rebal_dates$n <- 1:nrow(rebal_dates)

# re-order columns to make easier to read
rebal_dates <- rebal_dates[,
  c("n", "sample.sdt", "sample.edt", "sample.days", "rebal.sdt", "rebal.edt", "rebal.days")]

print("Sample and rebalance date pairs:")
print(rebal_dates)



# create benchmark 60/40 portfolio
bench_wgts_v <- vector("numeric", length(etfs))
names(bench_wgts_v) <- etfs

bench_wgts_v[names(bench_wgts_v) == "SPY"] <- 0.42
bench_wgts_v[names(bench_wgts_v) == "MDY"] <- 0.095
bench_wgts_v[names(bench_wgts_v) == "IJR"] <- 0.085
bench_wgts_v[names(bench_wgts_v) == "AGG"] <- 0.4

print("Benchmark portfolio:")
print(bench_wgts_v)



# Global Variable: create matrix to hold resulting portfolio weights, and benchmark portfolio characteristics
#
# NOTE: it is much more performant to pre-allocate an object and index into it iterativly versus 
# "growing" the object over each iteration. See Patrick Burn's R Inferno Chapter 2.

results_mtx <<- matrix(0.0, nrow = max(rebal_dates$n), ncol = (length(etfs)+5))
colnames(results_mtx) <- c(etfs, "bench.vol", "port.vol", "ar", "to", "te")
rownames(results_mtx) <- rebal_dates$n

# Global Variable: create matrix to hold per-ticker cumulative returns in sample
cr_mtx <<- matrix(0.0, nrow = max(rebal_dates$n), ncol = length(etfs))
colnames(cr_mtx) <- etfs 
rownames(cr_mtx) <- rebal_dates$n

print(cr_mtx[1:3, ])
print(results_mtx[1:3,])



# iterate through each re-balance period building a risk model,
# optimizing the portfolio and saving relevant characteristics into
# the results matrix

for(i in 1:max(rebal_dates$n)) {
 
  print(rebal_dates[i, ]$n)
  
  # subset out our return data and take it from long to wide
  # to coerce NAs out of any implicit missing data
  ohlc_w_df <- ohlc_df %>%
    filter(ref_date >= rebal_dates[i, ]$sample.sdt &
      ref_date < rebal_dates[i, ]$sample.edt) %>%
    select(dt = ref_date, ticker, r = ret_closing_prices) %>%
    arrange(dt, match(ticker, etfs)) %>%
    pivot_wider(names_from = "ticker", values_from = "r")
  
  # bail if we don't have complete cases 
  if(any(apply(ohlc_w_df[, -1], 2, function(x) { sum(is.na(x)) }) > 0)) {
    stop("Missing return data")
  }
 
  # calculate per-ticker cumulative return during sample period 
  cr_v <- apply(ohlc_w_df[, -1], 2, function(R) { prod(1 + R) - 1})
  cr_mtx[i, ] <- cr_v
  
  # create covariance matrix and annualize it
  Sigma <- as.matrix(ohlc_w_df[, -1])
  Sigma <- apply(Sigma, 2, function(x) { x - mean(x)})
  Sigma <- cov(Sigma)
  Sigma <- Sigma * 252
 
  # calculate the benchmark's volatility to be used as constraint 
  bench_vol_i <- ohlc_df %>%
    filter(ref_date >= rebal_dates[i, ]$sample.sdt &
      ref_date < rebal_dates[i, ]$sample.edt) %>%
    select(dt = ref_date, ticker, r = ret_closing_prices) %>%
    arrange(dt, match(ticker, etfs)) %>%
  tq_portfolio(., ticker, r, bench_wgts_v, col_rename = "r") %>%
  summarise(bench.volatility = sqrt(252 * var(r))) %>%
  pull()
  
  # check if our covariance matrix is positive definite
  #   a requirement for convex optimization problems, bail if not
  
  if(!matrixcalc::is.positive.definite(Sigma)){
    stop("Covariance matrix is not positive definite")
  }

  # no turnover constraint on initial optimization, setup up CVXR problem 
  # differently depending on iteration number
  if(i == 1) {
    
    # portfolio weights
    x <- Variable(ncol(Sigma))
 
    # objective function:
    #   maximize expected returns expressed via historical performance
    #   penalized by residual variance
    # constraints:
    #   fully invested - weights must sum to 100%
    #   long only - no weights less than 0.0
    #   volatility less than or equal to benchmark portfolio
    
    prob <- Problem(
      Minimize(
        -(t(cr_v) %*% x) +
          quad_form(x, Sigma)
       ),
      constraints = list(
        sum(x) == 1.0,
        x >= 0.0,
        p_norm((chol(Sigma) %*% x)) <= bench_vol_i 
      )
    )
    
  } else {
    
    # portfolio weights
    x <- Variable(ncol(Sigma))
    
    # turnover: portfolio weights - prior portfolio weights
    dx <- x - results_mtx[(i - 1), 1:length(etfs)]
 
    # objective function:
    #   maximize expected returns expressed via historical performance
    #   penalized by residual variance
    # constraints:
    #   fully invested - weights must sum to 100%
    #   long only - no weights less than 0.0
    #   volatility less than or equal to benchmark portfolio
    #   turnover: less than or equal to 20%
     
    prob <- Problem(
      Minimize(
        -(t(cr_v) %*% x) + 
          quad_form(x, Sigma)
       ),
      constraints = list(
        sum(x) == 1.0,
        x >= 0.0,
        p_norm((chol(Sigma) %*% x)) <= bench_vol_i,
        p_norm(dx, 1) <= 0.20
      )
    )
    
  }

  # attempt to solve problem and extract weights
  solution <- solve(prob)
  
  # if solution is infeasible just keep last re-balance periods weights
  # otherwise get weight vector from CVXR solution
  if(solution$status == "infeasible") {
    
    opt_w <- results_mtx[i-1, etfs]
    
  } else {
    
    opt_w <- as.vector(solution$getValue(x))
    names(opt_w) <- colnames(Sigma)
    
  }
 
  # set tiny weights to 0.0 and re-scale
  opt_w[opt_w < 0.0001] <- 0.0
  opt_w <- opt_w / sum(opt_w)
 
  opt_w
  cr_v[rev(order(cr_v))]
  
  # volatility
  vol_i <- as.numeric(sqrt(t(opt_w) %*% Sigma %*% opt_w))
 
  # active risk 
  ar_i <- sum(abs(bench_wgts_v - opt_w))
  
  # turnover
  if(i == 1) {
    to_i <- sum(abs(opt_w - rep(0, length(etfs))))
  } else {
    to_i <- sum(abs(opt_w - results_mtx[(i - 1), 1:length(etfs)]))
  }
  
  # tracking error
  te_i <- sqrt(t((opt_w - bench_wgts_v)) %*% Sigma %*% (opt_w - bench_wgts_v))
    
  # save results
  results_mtx[i, ] <- c(
    opt_w, bench_vol_i, vol_i, ar_i, to_i, te_i)  
  print(results_mtx[1:i, ])

  # clean up (should be unnecessary due to scope)
  rm(x)  
  rm(prob)
  
  
}

View(results_mtx)




# get an aggregated set of portfolio and benchmark returns
aggreg_rets <- map_dfr(rebal_dates$n, function(n) {
  
  print(n)
  
  port_rets_n <- ohlc_df %>%
    filter(ref_date >= rebal_dates[n, ]$rebal.sdt &
      ref_date < rebal_dates[n, ]$rebal.edt) %>%
    select(dt = ref_date, ticker, r = ret_closing_prices) %>%
    tq_portfolio(ticker, r, results_mtx[n, 1:length(etfs)], col_rename = "port.returns")
  
  bench_rets_n <- ohlc_df %>%
    filter(ref_date >= rebal_dates[n, ]$rebal.sdt &
      ref_date < rebal_dates[n, ]$rebal.edt) %>%
    select(dt = ref_date, ticker, r = ret_closing_prices) %>%
    tq_portfolio(ticker, r, bench_wgts_v, col_rename = "bench.returns")
  
  as.data.frame(inner_join(port_rets_n, bench_rets_n, by = "dt"))
  
})

# charts are nice
aggreg_rets_xts <- xts(aggreg_rets[, 2:3], order.by = aggreg_rets[, 1])
chart.CumReturns(aggreg_rets_xts, legend.loc = "bottom")
chart.RelativePerformance(aggreg_rets_xts[, 1], aggreg_rets_xts[, 2])

# glue the portfolio and dates together for extra nice charts
portfolio <- cbind(rebal_dates[, c("rebal.sdt", "rebal.edt")], data.frame(results_mtx))
View(cbind(rebal_dates[, c("rebal.sdt", "rebal.edt")], data.frame(cr_mtx)))

portfolio[, c("rebal.sdt", etfs)] %>%
  pivot_longer(where(is.numeric), names_to = "ticker", values_to = "wgt") %>%
  mutate(date = factor(rebal.sdt)) %>%
  ggplot(aes(x = date, y = wgt, fill = ticker, group = ticker)) +
    geom_area(alpha = 0.8, size=.5, colour = "white") + 
    theme_tq() +
    scale_fill_tq() + 
    labs(title = "Weights") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ylab("Weight") +
    xlab("Date")

# a static portfolio vs. our optimized one  
static_portfolio_w <- c(0.31, 0.16, 0.12, 0.20, 0.10, 0.11)  
sum(static_portfolio_w)
names(static_portfolio_w) <- etfs
print(static_portfolio_w)

# mean ETF weights over whole time frame
apply(portfolio[, etfs], 2, mean)

# look at the static portfolio's performance
ohlc_df %>%
  filter(ref_date >= sdt &
    ref_date < edt) %>%
  select(dt = ref_date, ticker, r = ret_closing_prices) %>%
  tq_portfolio(ticker, r, static_portfolio_w, col_rename = "static.port.returns") %>%
  summarise(cr = Return.cumulative(static.port.returns)) %>%
  as.data.frame()

# look at our portfolio and benchmark's performance
aggreg_rets %>% 
  pivot_longer(where(is.numeric), names_to = "name", values_to = "ret") %>%
  group_by(name) %>%
  summarise(cr = Return.cumulative(ret)) %>%
  as.data.frame()
